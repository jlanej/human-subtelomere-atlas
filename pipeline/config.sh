#!/bin/bash
# Shared pipeline configuration. Machine-specific values belong in work/config.local.sh (git-ignored), e.g.
#   export ATLAS_INPUT_DIR=/path/to/analysis_home   # directory holding tipdup/, subtel/, results/, ref/, hg002/, meta/
#   export ATLAS_THREADS=4
ATLAS_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export ATLAS_REPO
[ -f "$ATLAS_REPO/work/config.local.sh" ] && . "$ATLAS_REPO/work/config.local.sh"
: "${ATLAS_INPUT_DIR:?set ATLAS_INPUT_DIR (analysis home with tipdup/, subtel/, results/) in work/config.local.sh}"
export ATLAS_THREADS="${ATLAS_THREADS:-3}"
export ATLAS_S3_HTTPS="${ATLAS_S3_HTTPS:-https://human-pangenomics.s3.amazonaws.com}"
