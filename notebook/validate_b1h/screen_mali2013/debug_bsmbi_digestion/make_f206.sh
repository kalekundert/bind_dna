#!/usr/bin/env bash
set -euo pipefail

# Simultaneously make and optimize the PCR annealing temperature for f206.

sw pcr -u f206 -g 15 -n 9 -v 16 |

sw page_purify f206 -D -P tbe/purify -v 16
