#!/usr/bin/env bash
set -euo pipefail

codon_harmony --input linkers.prot --output linkers.fa --max-relax=0.3
