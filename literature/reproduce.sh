#!/usr/bin/env bash
# Reproduce the literature review track of the human subtelomere atlas.
# Requirements: python3 (standard library only); network access to api.openalex.org,
# api.crossref.org, eutils.ncbi.nlm.nih.gov and doi.org.
# Environment:
#   OPENALEX_API_KEY  OpenAlex API key (required by OpenAlex for every request)
#   ATLAS_INPUT_DIR   analysis home containing results/subtel_*.tsv (read from work/config.local.sh if unset)
# Raw API responses are cached in work/literature/ (git-ignored); pass --refresh to
# fetch_literature.py to re-query the live services (results may then differ from the release).
set -euo pipefail
cd "$(dirname "$0")/.."
: "${OPENALEX_API_KEY:?set OPENALEX_API_KEY}"
if [ -z "${ATLAS_INPUT_DIR:-}" ] && [ -f work/config.local.sh ]; then source work/config.local.sh; fi
: "${ATLAS_INPUT_DIR:?set ATLAS_INPUT_DIR (work/config.local.sh) to the analysis home}"
python3 literature/scripts/fetch_literature.py "$@"
python3 literature/scripts/screen_candidates.py
python3 literature/scripts/fetch_pubmed_abstracts.py
python3 literature/scripts/build_records.py
python3 literature/scripts/evidence_numbers.py > /dev/null
python3 literature/scripts/render_docs.py
python3 literature/scripts/write_readme.py
( cd literature && sha256sum -c SHA256SUMS 2>/dev/null || shasum -a 256 -c SHA256SUMS )
