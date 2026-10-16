#!/usr/bin/env python3
"""Collect the atlas-data values quoted in literature/NOVELTY.md and literature/REVIEW.md.

Inputs : ATLAS_INPUT_DIR (environment variable, e.g. from the git-ignored work/config.local.sh): the analysis
         home with results/subtel_*.tsv; and the final outputs of the atlas tracks in this repository
         (analyses/telomere_cis, analyses/junctions, analyses/allele_ages, data/catalogue revision 2,
         data/validation/tan2026_summary.json, and, if present, tests/validation summary tables).
Output : literature/evidence_numbers.json  {key: {value, source}} plus '_sources' (sha256 and mtime of every file read).
Sources are given relative to ATLAS_INPUT_DIR ('input:') or to the repository root ('repo:').
Analysis-track values (telomere_cis, junctions, allele_ages) were computed on the earlier per-window form calls;
catalogue values are from catalogue v0.1 revision 2 (paralog-family labels, nativity classes).
"""
import csv, datetime, hashlib, json, os, sys
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IN = os.environ.get("ATLAS_INPUT_DIR")
if not IN:
    sys.exit("set ATLAS_INPUT_DIR (source work/config.local.sh) to the analysis home")
E = {}
def path(src):
    pref, rel = src.split(":", 1)
    return os.path.join(IN if pref == "input" else REPO, rel)
def tsv(src):
    with open(path(src)) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))
def js(src):
    return json.load(open(path(src)))
def put(k, v, src):
    E[k] = {"value": v, "source": src}
def pct(a, b, nd=1):
    return round(100.0 * a / b, nd)

# --- upstream placement, validation and pedigree (analysis home) ---
S = "input:results/subtel_foreign_tips.tsv"
tips = tsv(S); nf = sum(r["t20_state"] == "foreign" for r in tips)
put("named_tips", len(tips), S); put("foreign_tips_t20", nf, S); put("foreign_pct_t20", pct(nf, len(tips)), S)
put("n_forms", len(tsv("input:results/subtel_foreign_forms.tsv")), "input:results/subtel_foreign_forms.tsv")
fs = [r for r in tips if r["t20_state"] == "foreign"]
s = [r for r in fs if r["form"] == "20q>3q"]
put("form_20q3q_tips", len(s), S); put("form_20q3q_samples", len({r["sample"] for r in s}), S)
jt = [r for r in tips if r.get("junction")]
put("junction_tips", len(jt), S)
put("junction_tips_shared", sum(1 for r in jt if r["junction_cluster_samples"] and float(r["junction_cluster_samples"]) >= 2), S)
S = "input:results/subtel_verkko_concordance.tsv"
vf = [r for r in tsv(S) if r["t20_state"] == "foreign"]
put("verkko_samples", len({r["sample"] for r in tsv(S)}), S)
put("verkko_foreign_tips", len(vf), S)
put("verkko_full_length", sum(1 for r in vf if r["aligned_bp"] and float(r["aligned_bp"]) / float(r["query_bp"]) >= 0.95 and float(r["ident"]) >= 0.995), S)
put("verkko_junction_concordant", sum(r["junction_concordant"] == "True" for r in vf), S)
put("verkko_junction_tested", sum(r["junction_concordant"] in ("True", "False") for r in vf), S)
S = "input:results/subtel_hg002_junction_reads.tsv"
cnt = {}
for r in tsv(S):
    if r["dataset"] == "ont_ul_dorado" and r["crosses_10kb"] == "True" and r["conflict"] == "False" and float(r["mapq"]) >= 20:
        cnt[r["q"]] = cnt.get(r["q"], 0) + 1
put("hg002_switches", len(cnt), S); put("hg002_ont_cross_min", min(cnt.values()), S); put("hg002_ont_cross_max", max(cnt.values()), S)
S = "input:results/subtel_pedigree_transmission_summary.tsv"
ped = [r for r in tsv(S) if r["window"] == "100k"]
put("ped_foreign", sum(int(r["n_foreign"]) for r in ped), S); put("ped_foreign_intact", sum(int(r["n_foreign_intact"]) for r in ped), S)
put("ped_transmissions", len(ped), S)

# --- allele states by window (telomere_cis fixed windows) ---
S = "repo:analyses/telomere_cis/results/fixed_window_states.tsv"
fw = tsv(S)
for w in ("S5", "S10", "S20"):
    put(f"foreign_pct_{w}", pct(sum(r[f"{w}_state"] == "foreign" for r in fw), len(fw)), S)
put("fixed_window_tips", len(fw), S)

# --- catalogue (v0.1 revision 2) ---
S = "repo:data/catalogue/alleles_structural.tsv"
al = tsv(S); put("cat_alleles", len(al), S)
for t in "ABC":
    put(f"cat_tier_{t}", sum(r["tier"] == t for r in al), S)
for c in ("native", "within_family", "foreign"):
    put(f"cat_alleles_{c}", sum(r["nativity_class"] == c for r in al), S)
put("cat_haplogroups", len(tsv("repo:data/catalogue/haplogroups.tsv")), "repo:data/catalogue/haplogroups.tsv")
S = "repo:data/catalogue/paralog_families.tsv"
pf = tsv(S); put("cat_families", len(pf), S)
put("cat_family_list", "; ".join(f"{r['family']} ({r['members'].replace(',', ', ')})" for r in pf), S)
S = "repo:data/catalogue/ends.tsv"
en = tsv(S); r2 = [r for r in en if r["set"] == "HPRC_R2"]
put("cat_ends_r2", len(r2), S); put("cat_ends_chm13", sum(r["set"] == "CHM13" for r in en), S)
put("cat_ends_hg002", sum(r["set"] == "HG002_Q100" for r in en), S)
for c in ("native", "within_family", "foreign"):
    put(f"cat_r2_{c}", sum(r["nativity_class"] == c for r in r2), S)
put("cat_r2_nonnative_pct", pct(sum(r["nativity_class"] in ("within_family", "foreign") for r in r2), len(r2)), S)
put("cat_r2_foreign_pct", pct(sum(r["nativity_class"] == "foreign" for r in r2), len(r2)), S)
S = "repo:data/catalogue/QC.md"
for line in open(path(S)):
    c_ = [x.strip() for x in line.strip().strip("|").split("|")]
    if len(c_) >= 6 and c_[1] == "non-native" and c_[0] in ("-500", "500"):
        put("cat_shift_nonnative_" + ("minus500" if c_[0] == "-500" else "plus500"), round(100 * float(c_[3].split()[0]), 1), S + " (section 9 table, allele column)")
        put("cat_shift_n_nonnative", int(c_[2]), S + " (section 9 table)")

# --- telomere_cis ---
S = "repo:analyses/telomere_cis/results/baseline_stats.json"
b = js(S); put("tl_karimian_rho", round(b["mean_rel_tl"]["spearman_rho"], 2), S); put("tl_karimian_arms", b["mean_rel_tl"]["n"], S)
put("tl_karimian_p", f"{b['mean_rel_tl']['spearman_p']:.1e}", S)
S = "repo:analyses/telomere_cis/results/paired_tests_t20.tsv"
t20 = next(r for r in tsv(S) if r["test"] == "pooled_all_het_pairs")
put("tl_t20_d", round(float(t20["mean_d"])), S); put("tl_t20_p", round(float(t20["wilcoxon_p"]), 3), S)
S = "repo:analyses/telomere_cis/results/window_robustness.tsv"
wr = {r["test"]: r for r in tsv(S)}
for w in ("S20", "S10", "S5", "W515"):
    put(f"tl_{w}_d", round(float(wr[w]["mean_d"])), S); put(f"tl_{w}_p", float(f"{float(wr[w]['wilcoxon_p']):.2g}"), S)
ext = [wr[f"{w} | boundary-robust TL_ext"] for w in ("S20", "S10", "S5", "W515")]
put("tl_ext_min", round(min(float(r["mean_d"]) for r in ext)), S); put("tl_ext_max", round(max(float(r["mean_d"]) for r in ext)), S)
put("tl_ext_min_p", round(min(float(r["wilcoxon_p"]) for r in ext), 2), S)
put("tl_S20_pairs", int(wr["S20"]["n"]), S); put("tl_S20_mde", round(float(wr["S20"]["mde80"])), S)
S = "repo:analyses/telomere_cis/results/window_robustness.json"
w = js(S)
put("tl_forms_tested", w["forms_tl_ext"]["n_forms"], S); put("tl_forms_sig_ext", w["forms_tl_ext"]["n_bh_q_lt_0_05"], S)
d = w["donor_likeness_tl_ext_S20_pairs"]
put("tl_donor_slope_ext", round(d["slope"], 2), S); put("tl_donor_se_ext", round(d["se"], 2), S); put("tl_donor_p_ext", round(d["p"], 2), S)
put("tl_its_shift_per100", round(-w["boundary_adjusted"]["S5"]["beta_dits_per_100bp"]), S)
put("tl_S5_minus_W515", round(w["window_heterogeneity_bootstrap"]["S5-W515"]["diff"]), S)
put("tl_S5_minus_W515_p", w["window_heterogeneity_bootstrap"]["S5-W515"]["p_boot"], S)
S = "repo:analyses/telomere_cis/results/hg002_stats.json"
h = js(S)["het_arms"]["S20|ont_ul_dorado"]
put("tl_hg002_het_arms", h["n_arms"], S); put("tl_hg002_foreign_longer", h["n_foreign_longer"], S); put("tl_hg002_sign_p", h["sign_test_p"], S)

# --- junctions ---
S = "repo:analyses/junctions/results/reciprocal_summary.json"
rs = js(S); pc = rs["pooled"]
put("recip_carriers_assessed", pc["carriers"]["assessed"], S); put("recip_carriers_reciprocal", pc["carriers"]["reciprocal"], S)
put("recip_pct_carriers", round(100 * pc["frac"]["carriers"]["reciprocal"], 2), S)
put("recip_pct_noncarriers", round(100 * pc["frac"]["noncarriers"]["reciprocal"], 2), S)
put("recip_chi2_p", round(pc["chi2_donor_status"]["p"], 2), S); put("recip_forms", rs["forms"], S)
put("recip_complement_matched", rs["complement_clusters_matched"], S); put("recip_complement_tested", rs["complement_clusters_tested"], S)
S = "repo:analyses/junctions/results/junction_summary.json"
j = js(S); o = j["mh_classes_obs"]; jc = j["junction_class"]
put("jx_clusters", j["analysed"], S)
put("jx_gt25", o.get("homology_26_200", 0) + o.get("homology_gt200", 0), S)
put("jx_null_pct_gt25", round(100 * (j["mh_classes_null_frac"].get("homology_26_200", 0) + j["mh_classes_null_frac"].get("homology_gt200", 0)), 3), S)
put("jx_spanning", jc["homology_spanning"], S); put("jx_spanning_pct", pct(jc["homology_spanning"], j["analysed"], 0), S)
put("jx_nonhom", j["analysed"] - jc["homology_spanning"], S)
put("jx_blunt", jc["blunt"], S); put("jx_micro", jc["microhomology"], S); put("jx_shorthom", jc["short_homology"], S)
put("jx_far", j["far_from_telomere"], S)
S = "repo:analyses/junctions/results/junction_null_summary.tsv"
nul = {(r["subset"], r["feature"]): r for r in tsv(S)}
m = nul[("all", "mh_parental")]
put("jx_mh_obs_mean", round(float(m["obs_mean"]), 1), S); put("jx_mh_null_mean", round(float(m["null_mean"]), 2), S); put("jx_mh_perm_p", float(m["perm_p_greater"]), S)
a = nul[("all", "its_run_c_any")]; f = nul[("far", "its_run_c_any")]
put("its_all_pct", round(100 * float(a["obs_mean"]), 1), S); put("its_all_null_pct", round(100 * float(a["null_mean"]), 2), S)
put("its_far_pct", round(100 * float(f["obs_mean"]), 2), S); put("its_far_null_pct", round(100 * float(f["null_mean"]), 2), S)
put("its_far_p", round(float(f["perm_p_greater"]), 2), S)

# --- allele_ages ---
S = "repo:analyses/allele_ages/results/pop_summary.json"
p = js(S)
put("pop_forms_tested", p["forms_tested"], S); put("pop_forms_q05", p["forms_q05_perm"], S)
put("pop_forms_all5", p["forms_in_all_5_superpops"], S); put("pop_median_fst", round(p["median_fst_anova"], 3), S)
put("pop_big_clusters", p["junction_clusters_ge10_samples"], S); put("pop_big_clusters_ooa", p["junction_clusters_ge10_shared_AFRcont_OOA"], S)
S = "repo:analyses/allele_ages/results/origin_summary.json"
og = js(S)
put("age_forms_analysed", og["forms_analysed"], S); put("age_forms_clustered_crit", og["clustered"], S)
put("age_single_origin", og["single_origin_counts"]["consistent with single origin"], S); put("age_monophyletic", og["strictly_monophyletic"], S)
put("age_forms_selected", og["forms_selected"], S); put("age_forms_untestable", og["forms_selected"] - og["forms_analysed"], S)
sc_ = og["single_origin_counts"]
put("age_clustered_diverse", sc_["clustered, carriers as diverse as natives"], S)
put("age_clustered_lowdiv", sc_["clustered, low diversity, no carrier-fixed SNVs"], S)
put("age_not_supported", sc_["not supported"], S)
put("age_afr_higher", og["afr_vs_nonafr_carrier_pi"]["afr_higher"], S); put("age_afr_forms", og["afr_vs_nonafr_carrier_pi"]["forms"], S)
put("age_afr_p", round(og["afr_vs_nonafr_carrier_pi"]["p_sign"], 3), S)
put("age_single_origin_mean_age_kyr", round(og["consistent_median_age_kyr_mu"]), S)
ape = og["ape"]
put("ape_clusters", ape["junction_clusters_tested"], S); put("ape_shared", ape["candidate_shared_junctions"], S)
put("ape_hits", ape["carrier_probe_spanning_hits"], S); put("ape_hit_clusters", ape["clusters_with_spanning_hit"], S)
put("ape_hits_duplicon", ape["diagnosis"].get("shared_duplicon", 0), S)
S = "repo:analyses/allele_ages/results/origin_verdicts.tsv"
ov = [r for r in tsv(S) if r["status"] == "ok"]
put("age_nn_sig", sum(float(r["p_nn"]) <= 0.001 for r in ov), S)
cl = [r for r in ov if r["clustered"] == "True"]
put("age_clustered", len(cl), S)
put("age_clustered_02_15", sum(200 <= float(r["age_kyr_median_mu"]) <= 1500 for r in cl), S)
for r in ov:
    if r["form"] in ("20q>3q", "3q>11p", "6p>11p"):
        put("age_" + r["form"].replace(">", ""), round(float(r["age_kyr_median_mu"])), S)

# --- cross-check against Tan et al. 2026 sequences (data/validation) ---
S = "repo:data/validation/tan2026_summary.json"
t = js(S); sh = t["same_haplotype"]; pr = sh["primary"]; ic = t["independent_cohorts"]
put("tan_primary_pairs", pr["n_pairs"], S); put("tan_samples", sh["n_samples"], S)
put("tan_sa_agree", round(100 * pr["sa_agree"], 1), S); put("tan_sa_n", pr["n_sa_compared"], S)
put("tan_s20_agree", round(100 * pr["S20_state_and_family_agree"], 1), S)
put("tan_switch_identical", round(100 * pr["switch_base_frac_proj_exact"], 1), S); put("tan_switch_n", pr["n_switches_base_both"], S)
nb = sh["by_nativity_class"]
for c in ("native", "foreign", "within_family"):
    put(f"tan_sa_agree_{c}", round(100 * nb[c]["sa_agree"], 1), S)
put("tan_cohort_records", ic["n_primary"], S); put("tan_cohort_samples", ic["n_primary_samples"], S)
fa = ic["frac_assigned_by_cohort"]
put("tan_cohort_assigned_min", round(100 * min(fa.values()), 1), S); put("tan_cohort_assigned_max", round(100 * max(fa.values()), 1), S)
rp = ic["replication_tierA_common"]
put("tan_rep_n", rp["n_alleles"], S); put("tan_rep_obs", rp["n_observed"], S); put("tan_rep_obs3", rp["n_observed_ge_3_samples"], S)
put("tan_rep_by_cohort", ", ".join(f"{k} {rp[k]['n_observed']}/{rp[k]['n_alleles']}" for k in ("HGSVC3", "CPC", "APC", "JASA")), S)
ab = ic["tierAB_observed_in_independent"]
put("tan_tierA_obs", f"{ab['A']['observed']}/{ab['A']['n']}", S); put("tan_tierB_obs", f"{ab['B']['observed']}/{ab['B']['n']}", S)
fc = ic["freq_correlation_non_native"]
put("tan_freq_r", round(fc["pearson_r"], 3), S); put("tan_freq_n", fc["n_alleles"], S)
put("tan_sig_sa", ic["n_sa_tests_q_lt_0.05"], S); put("tan_n_sa", ic["n_sa_tests"], S)
put("tan_sig_form", ic["n_form_tests_q_lt_0.05"], S); put("tan_n_form", ic["n_form_tests"], S)
nv = ic["novel"]
put("tan_novel_records", nv["n_unassigned_primary"], S); put("tan_novel_pct", round(100 * nv["frac_unassigned_primary"], 1), S)
put("tan_novel_candidates", nv["n_candidate_alleles"], S)
S = "repo:data/validation/tan2026_cohort_forms.tsv"
cf = [r for r in tsv(S) if r["hprc_reference"] != "ALL"]
put("tan_matched_tests", len(cf), S)
put("tan_matched_sig", sum(1 for r in cf if r["bh_q"] not in ("", "nan") and float(r["bh_q"]) < 0.05), S)

# --- query-tool validation (tests/validation), read if present ---
S = "repo:tests/validation/self_consistency_summary.tsv"
L = "repo:tests/validation/loso_summary.tsv"
if os.path.exists(path(S)) and os.path.exists(path(L)):
    sc = next(r for r in tsv(S) if r["stratum"] == "all" and r["metric"] == "sa_match")
    lo = next(r for r in tsv(L) if r["allele_type"] == "structural_allele" and r["stratum"] == "all")
    put("tool_sentence", f"The query tool reproduces the catalogue structural allele for {sc['n_agree']} of {sc['n']} re-queried catalogue windows "
        f"and assigns {lo['n_correct']} of {lo['n']} ends correctly in leave-one-sample-out tests (tests/validation/VALIDATION.md).", S + "; " + L)
else:
    E["tool_sentence"] = {"value": "Query-tool validation results (tests/validation/VALIDATION.md) were not available when this file was rendered.", "source": S + " (absent)"}

snap = {}
for src in sorted({v["source"].split(" ")[0] for v in E.values()}):
    pth = path(src)
    if os.path.isfile(pth):
        snap[src] = {"sha256": hashlib.sha256(open(pth, "rb").read()).hexdigest(),
                     "mtime_utc": datetime.datetime.fromtimestamp(os.path.getmtime(pth), datetime.timezone.utc).strftime("%Y-%m-%dT%H:%MZ")}
E["_sources"] = {"value": snap, "source": "sha256 and modification time of each source file at extraction"}
json.dump(E, open(os.path.join(REPO, "literature", "evidence_numbers.json"), "w"), indent=1)
print(json.dumps({k: v["value"] for k, v in E.items() if not k.startswith("_")}))
