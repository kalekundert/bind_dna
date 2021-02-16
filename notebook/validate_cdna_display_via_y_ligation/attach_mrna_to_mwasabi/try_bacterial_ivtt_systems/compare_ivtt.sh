#!/usr/bin/env bash
set -euo pipefail

sw zap |
sw purex f89 -t '30 min' -v 6 -r |
sw nebex f89 -t '30 min' 
