#!/usr/bin/env bash
set -euo pipefail

sw ivtt p107 |
sw measure_protein_conc "PURExpress product" -v 5 -s 4 -x 0.1


