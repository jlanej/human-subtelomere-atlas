#!/bin/bash
# Shared pipeline configuration. Machine-specific values belong in work/config.local.sh (git-ignored), e.g.
#   export ATLAS_INPUT_DIR=/path/to/analysis_home   # optional: original analysis home (tipdup/, subtel/, results/, ref/, hg002/, meta/)
#   export ATLAS_THREADS=4
# Without ATLAS_INPUT_DIR the pipeline runs in public mode: small derived inputs from data/inputs/ (committed snapshot) and large
# public files downloaded by pipeline/fetch_public_inputs.sh into work/public/ (see pipeline/INPUTS.md).
ATLAS_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export ATLAS_REPO
[ -f "$ATLAS_REPO/work/config.local.sh" ] && . "$ATLAS_REPO/work/config.local.sh"
if [ -n "${ATLAS_INPUT_DIR:-}" ]; then export ATLAS_INPUT_DIR; else unset ATLAS_INPUT_DIR; export ATLAS_INPUT_MODE=public; fi
export ATLAS_THREADS="${ATLAS_THREADS:-3}"
export ATLAS_S3_HTTPS="${ATLAS_S3_HTTPS:-https://human-pangenomics.s3.amazonaws.com}"
