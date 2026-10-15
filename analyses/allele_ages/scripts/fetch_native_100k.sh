#!/bin/bash
# Range-read the regions listed by fetch_native_100k.py from the public HPRC S3 assemblies (samtools faidx over https;
# htslib fetches .fai/.gzi). Idempotent: assemblies with a non-empty .ok are skipped. Parallelism: ATLAS_THREADS (default 3).
set -euo pipefail
cd "$(dirname "$0")/../../../work/allele_ages/r2_100k_native"
P=${ATLAS_THREADS:-3}
cat jobs.txt | xargs -P "$P" -L 1 sh -c '[ -s $0.ok ] && exit 0; for t in 1 2 3; do samtools faidx -r $0.regions "$1" > $0.fa 2> $0.log && [ $(grep -c ">" $0.fa) -eq $(wc -l < $0.regions) ] && { date -u +%Y-%m-%dT%H:%M:%SZ > $0.ok; exit 0; }; sleep 5; done; echo FAILED $0 >> failed.txt'
echo "done: $(ls *.ok | wc -l) of $(wc -l < jobs.txt) assemblies"
