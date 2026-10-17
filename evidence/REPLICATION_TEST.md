# Replication test of `evidence/replicate_claim.sh`

**Result: 9 of 9 runs matched exactly.** Each run covered every cited read, the claim
verdict (v1) and the definitive tier. The runs were 4 HG002 Q100 claims, 4 1000 Genomes ONT claims (one confirmed but not
definitive) and 1 strictly native control.

## Setup

* **Clean tree.** A copy of the repository without `work/`, `.git/` and `evidence/showcase/` (`rsync --exclude work/ ...`). No file
  produced by the evidence pipelines was available to the runs other than the committed tables (`evidence/claims.tsv`,
  `evidence/reads.tsv.gz`, `evidence/controls.tsv`) and the catalogue (`data/catalogue/`).
* **Fresh cache** (`SUBTEL_REPLICATE_CACHE`, outside the repository). CHM13 v2.0 came from its public URL
  (`https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz`, 981 MB + `.gzi`;
  73 s). The `.fai` is built locally because the server has none. The 48 CHM13 ends came from `subtel-atlas fetch-refs` (range
  reads, SHA-256-checked; about 30 s on the first run, then cached). Assembly windows were cut by range from the public HPRC
  Release 2 and HG002 Q100 FASTA files and checked against `ends.tsv` `window_sha256`. Reads came from the public HG002 Q100 BAMs
  and the 1000 Genomes ONT CRAMs, by region and read ID.
* **Machine.** One laptop (macOS, 10 cores), two runs in parallel with 3 threads each. Software: minimap2/mappy
  2.31-r1302, samtools 1.24, Python 3.11.
* **Selection.** The 4 HG002 claims are the 3q, 9p and 17p showcase claims plus 15q paternal (2 definitive reads, single
  dataset). The 1000 Genomes ONT claims are the 3 showcase claims plus NA18960 2q SA-2q-02, which is v1-confirmed but not
  definitive (shared anchors). The control is the 1KG-ONT strictly native control with the most cited reads.
* **Not tested.** HG003/HG004 trio claims: their ONT reads are GRCh38 CRAMs and need a GRCh38 FASTA (`--grch38`).

## Results

Times are in seconds: fetch (range requests), stage 1 (diagnostic units and switch), genome mapping (minimap2 builds the CHM13
index in memory, about 45 s of each value), placement (end windows cut from the public assembly and reads placed), wall
(whole run, warm cache).

| claim | kind | cited reads | fetch requests | published (v1 / definitive) | replicated | reads identical | fetch | stage 1 | genome | placement | wall | result |
|---|---|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---|
| `HG002#2#chr3_MATERNAL#R` | HG002 Q100 | 136 | 7 | confirmed / definitive_multi_dataset (7 conclusive) | confirmed / definitive_multi_dataset (7 conclusive) | 136/136 | 14 | 206 | 52 | 50 | 328 | **MATCH** |
| `HG002#1#chr17_PATERNAL#L` | HG002 Q100 | 107 | 7 | confirmed / definitive_multi_dataset (7 conclusive) | confirmed / definitive_multi_dataset (7 conclusive) | 107/107 | 44 | 482 | 269 | 77 | 907 | **MATCH** |
| `HG002#1#chr15_PATERNAL#R` | HG002 Q100 | 36 | 10 | confirmed / definitive (2 conclusive) | confirmed / definitive (2 conclusive) | 36/36 | 73 | 83 | 60 | 43 | 295 | **MATCH** |
| `HG002#2#chr9_MATERNAL#L` | HG002 Q100 | 190 | 10 | confirmed / definitive_multi_dataset (74 conclusive) | confirmed / definitive_multi_dataset (74 conclusive) | 190/190 | 56 | 307 | 67 | 58 | 528 | **MATCH** |
| `HG01258|20q|SA-20q-01` | 1KG-ONT | 5 | 1 | confirmed / definitive (4 conclusive) | confirmed / definitive (4 conclusive) | 5/5 | 10 | 1 | 46 | 36 | 126 | **MATCH** |
| `HG01258|17p|SA-17p-02` | 1KG-ONT | 6 | 2 | confirmed / definitive (5 conclusive) | confirmed / definitive (5 conclusive) | 6/6 | 12 | 19 | 55 | 44 | 135 | **MATCH** |
| `NA20346|2q|SA-2q-04` | 1KG-ONT | 10 | 1 | confirmed / definitive (3 conclusive) | confirmed / definitive (3 conclusive) | 10/10 | 4 | 6 | 47 | 40 | 131 | **MATCH** |
| `NA18960|2q|SA-2q-02` | 1KG-ONT | 16 | 2 | confirmed / confirmed_not_definitive (8 conclusive) | confirmed / confirmed_not_definitive (8 conclusive) | 16/16 | 6 | 6 | 46 | 43 | 106 | **MATCH** |
| `control|HG00133|20p|SA-20p-10` | control (1KG-ONT) | 20 | 1 | 0 conclusive | 0 conclusive | 20/20 | 4 | 2 | 50 |  | 106 | **MATCH** |

Machine-readable: [replication_results.tsv](replication_results.tsv). Per-read comparisons are printed by the script and saved as
`replicate_summary.json` in the run directory.

## Example (abridged output)

```
$ evidence/replicate_claim.sh 'HG01258|20q|SA-20q-01'
[criteria: evidence/criteria_v1.json a89c5eaf023ed032, packaged a89c5eaf023ed032, required a89c5eaf023ed032]
model HG01258#2#...: D=chr3q P=chr20q R=chr20q pattern=DP
fetched 5 cited reads from HG01258.t2t.cram chr20:66205256-66210255
read                                  dataset          published   replicated  conc pub/rep  switch_ci_sub (rep)
30a25f04-4069-425c-aa44-15cc0de6b642  1KG_ONT_HG01258  conclusive  conclusive  T/T           22100-31700
...
claim HG01258|20q|SA-20q-01: v1 verdict replicated confirmed (published confirmed); conclusive concordant 4 (4) ...
definitive tier: 4 definitive reads (4 molecules ...) -> definitive (published definitive); per-read tiers match: True
RESULT: MATCH  (5/5 reads identical)
```

## Notes

* The first HG002 run stopped after the genome mapping on a bookkeeping error in the script (a missing-read check read a field
  that is set only in stage 2). The script was fixed and the run repeated; the table shows the repeat. No verdict logic was
  involved.
* Wall time is dominated by stage 1 for long HG002 ultra-long reads (up to about 500 kb), and by the in-memory CHM13 index
  (about 45 s per run). `SUBTEL_CHM13_MMI` points the script to a prebuilt index (`minimap2 -x map-ont -d chm13.mmi chm13v2.0.fa.gz`,
  7.5 GB on disk) and saves that time.
* The published verdicts were computed with the same minimap2/mappy version. Other versions can move individual alignments.
