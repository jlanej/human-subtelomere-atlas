#!/usr/bin/env python3
"""Per-form verdicts for single origin and relative age (inputs: results/origin_forms.tsv, form_freq_by_superpop.tsv,
ape_junction_hits_diagnosis.tsv). Criteria (fixed before inspection of individual forms):
  clustered       : nn_carrier_frac >= 0.9 and p_nn < 0.01 (carriers are each other's nearest haplotypes)
  diagnostic      : >= 1 SNV at alt frequency >= 0.9 in carriers and <= 0.1 in natives, or the smallest UPGMA clade holding
                    95% of carriers is >= 80% carriers
  single_origin   : 'consistent' if clustered and diagnostic and pi_carrier < pi_native (carriers less diverse than the donor
                    arm's native alleles); 'clustered, carriers as diverse as natives' if clustered but pi_carrier >= pi_native
                    (old single origin with subsequent conversion, or several origins of similar sequence); 'not supported'
                    otherwise. 'clustered, low diversity, no carrier-fixed SNVs' = clustered and pi_carrier < pi_native but not diagnostic
                    (carriers form a low-diversity lineage nested within native variation).
Relative age: pi_ratio = pi_carrier / pi_native (mu-independent). Neutral reference: expected age of a neutral allele at
  frequency p in units of 2Ne generations, -2p ln p / (1 - p) (Kimura & Ohta 1973, doi:10.1093/genetics/75.1.199), compared to
  the mean pairwise coalescence time of the sample (2Ne): a carrier/native ratio well below this reference indicates a
  younger lineage than expected for a neutral allele at that frequency (or natives that are older than average).
Clock age: rho (mean SNVs per carrier from the carrier consensus per bp) / mu, mu = 1.25e-8 per bp per generation,
  29 years per generation; star genealogy assumed. Upper-biased by paralogous gene conversion and alignment noise.
Robust clock age: median carrier distance to the consensus / mu (age_kyr_median_mu).
Out-of-Africa pattern: pi among AFR carriers > pi among non-AFR carriers (sign test across forms).
Outputs: results/origin_verdicts.tsv, results/origin_summary.json
"""
import os, json, numpy as np, pandas as pd
from scipy.stats import binomtest, wilcoxon
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
R = pd.read_csv(os.path.join(RES, 'origin_forms.tsv'), sep='\t')
F = pd.read_csv(os.path.join(RES, 'form_freq_by_superpop.tsv'), sep='\t')
F = F[F.superpopulation == 'ALL'].set_index('form')
R['freq'] = R.form.map(F.freq)
ok = R.status == 'ok'
R['clustered'] = ok & (R.nn_carrier_frac >= 0.9) & (R.p_nn < 0.01)
R['clade95_purity'] = (R.clade95_size - R.clade95_noncarriers) / R.clade95_size
R['diagnostic'] = ok & ((R.diag_snv_carrier >= 1) | (R.clade95_purity >= 0.8))
def verdict(r):
    if r.status != 'ok':
        return 'not analysed: ' + r.status
    if r.clustered and r.diagnostic and r.pi_carrier < r.pi_native:
        return 'consistent with single origin'
    if r.clustered and r.pi_carrier < r.pi_native:
        return 'clustered, low diversity, no carrier-fixed SNVs'
    if r.clustered:
        return 'clustered, carriers as diverse as natives'
    return 'not supported'
R['single_origin'] = R.apply(verdict, axis=1)
p = R.freq.clip(1e-6, 1 - 1e-6)
R['neutral_age_ref_2Ne'] = (-2 * p * np.log(p) / (1 - p)).round(3)
R['younger_than_neutral_ref'] = ok & (R.pi_ratio_carrier_native < R.neutral_age_ref_2Ne)
R['afr_carrier_pi_higher'] = np.where(R.pi_carrier_AFR.notna() & R.pi_carrier_nonAFR.notna(), R.pi_carrier_AFR > R.pi_carrier_nonAFR, np.nan)
# nearest non-carrier haplotype to the carrier consensus (results/origin_nearest_to_consensus.tsv)
N = pd.read_csv(os.path.join(RES, 'origin_nearest_to_consensus.tsv'), sep='\t')
nn1 = N.sort_values('dist_to_carrier_consensus_per_kb').groupby('form').head(1).set_index('form')
R['nearest_noncarrier_group'] = R.form.map(nn1.group)
R['nearest_noncarrier_form'] = R.form.map(nn1.member_form.fillna(''))
R['nearest_noncarrier_superpop'] = R.form.map(nn1.superpopulation)
R['nearest_noncarrier_dist_per_kb'] = R.form.map(nn1.dist_to_carrier_consensus_per_kb)
R['carrier_median_dist_per_kb'] = R.form.map(nn1.carrier_median_dist_to_consensus)
# donor lineage: 'nested' = a non-carrier is closer to the carrier consensus than the median carrier (the source lineage
# survives among non-carriers); 'distinct' = all sampled non-carriers are farther than 3x the median carrier distance
R['donor_lineage'] = np.where(R.nearest_noncarrier_dist_per_kb <= R.carrier_median_dist_per_kb, 'nested in sampled non-carriers',
                      np.where(R.nearest_noncarrier_dist_per_kb > 3 * R.carrier_median_dist_per_kb, 'distinct (no close non-carrier)', 'intermediate'))
R.loc[~(R.status == 'ok'), 'donor_lineage'] = np.nan
# robust clock age from the MEDIAN carrier distance to the consensus (less sensitive to a minority of divergent carriers,
# e.g. converted or mis-assigned tips): T = median / mu generations
MU, GEN = 1.25e-8, 29
R['age_kyr_median_mu'] = (R.carrier_median_dist_per_kb / 1000 / MU * GEN / 1000).round(1)
R.to_csv(os.path.join(RES, 'origin_verdicts.tsv'), sep='\t', index=False)
O = R[R.status == 'ok']
both = O.dropna(subset=['pi_carrier_AFR', 'pi_carrier_nonAFR'])
k = int((both.pi_carrier_AFR > both.pi_carrier_nonAFR).sum())
D = pd.read_csv(os.path.join(RES, 'ape_junction_hits_diagnosis.tsv'), sep='\t')
P = pd.read_csv(os.path.join(RES, 'ape_junction_probes.tsv'), sep='\t')
cons = O[O.single_origin == 'consistent with single origin']
summ = dict(
    forms_selected=len(R), forms_analysed=int(ok.sum()), not_analysed=R[~ok].status.value_counts().to_dict(),
    clustered=int(O.clustered.sum()), single_origin_counts=O.single_origin.value_counts().to_dict(),
    strictly_monophyletic=int((O.clade_all_noncarriers == 0).sum()),
    clade_all_le1_noncarrier=O[O.clade_all_noncarriers <= 1].form.tolist(),
    median_nn_carrier_frac=float(O.nn_carrier_frac.median()), median_nn_expected=float(O.nn_perm_expected.median()),
    median_pi_ratio=float(O.pi_ratio_carrier_native.median()),
    pi_ratio_lt1=int((O.pi_ratio_carrier_native < 1).sum()),
    younger_than_neutral_ref=int(O.younger_than_neutral_ref.sum()),
    consistent_forms=cons.form.tolist(),
    donor_lineage_counts=O.donor_lineage.value_counts().to_dict() if 'donor_lineage' in O else None,
    nearest_group_counts=O.nearest_noncarrier_group.value_counts().to_dict(),
    consistent_median_age_kyr_mu=float(cons.age_kyr_mu.median()) if len(cons) else None,
    consistent_median_age_kyr_median_mu=float(cons.age_kyr_median_mu.median()) if len(cons) else None,
    forms_age_median_mu_lt_300kyr=O[O.age_kyr_median_mu < 300].form.tolist(),
    afr_vs_nonafr_carrier_pi=dict(forms=len(both), afr_higher=k, p_sign=binomtest(k, len(both), 0.5).pvalue,
                                  p_wilcoxon=float(wilcoxon(both.pi_carrier_AFR, both.pi_carrier_nonAFR).pvalue)),
    mutation_model=dict(mu_per_bp_per_gen=1.25e-8, years_per_gen=29, estimator='rho (mean distance to carrier consensus), star genealogy'),
    ape=dict(junction_clusters_tested=int(P.junction_cluster.nunique()), species=int(P.species.nunique()), tests=len(P),
             carrier_probe_spanning_hits=len(D), clusters_with_spanning_hit=int(D.junction_cluster.nunique()),
             diagnosis=D.diagnosis.value_counts().to_dict(),
             hits_where_carrier_probe_spans_chm13_own_arm=int(D.carrier_probe_spans_chm13A.sum()),
             candidate_shared_junctions=int((D.diagnosis == 'candidate_shared_junction').sum())),
)
headline = ['20q>3q', '3q>11p', '6p>11p', '8p>1p', '16p>9q', '1p>11p']
summ['headline'] = R[R.form.isin(headline)][['form', 'carriers_used', 'natives_used', 'interval_bp', 'pi_carrier', 'pi_native',
    'pi_ratio_carrier_native', 'pi_ratio_ci95', 'nn_carrier_frac', 'p_nn', 'clade_all_noncarriers', 'clade95_noncarriers',
    'diag_snv_carrier', 'age_kyr_mu', 'age_kyr_ci95', 'age_kyr_median_mu', 'neutral_age_ref_2Ne', 'pi_carrier_AFR', 'pi_carrier_nonAFR', 'single_origin',
    'nearest_noncarrier_group', 'nearest_noncarrier_form', 'nearest_noncarrier_superpop', 'nearest_noncarrier_dist_per_kb',
    'carrier_median_dist_per_kb', 'donor_lineage']].to_dict('records')
json.dump(summ, open(os.path.join(RES, 'origin_summary.json'), 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in summ.items() if k != 'headline'}, indent=1, default=float))
print(pd.DataFrame(summ['headline']).to_string())
