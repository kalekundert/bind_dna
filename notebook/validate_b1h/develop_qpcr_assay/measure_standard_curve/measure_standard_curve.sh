#!/usr/bin/env bash
set -euo pipefail

sw b1h/qpcr_assay \
  sz224 sz228 sz229 sz230 sz231 \
  -P sr1/sr2 \
  -P sr3/sr5 \
  -d 3
