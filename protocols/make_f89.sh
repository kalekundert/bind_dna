#!/usr/bin/env bash
set -euo pipefail

# Once I make the improvements I have in mind for the `sw make` command, this 
# will be simply `sw make f42 f85 f89`.  Until then, though...
sw make f42 |
sw phenol |
sw ethanol -v 2 -b 'nuclease-free water' |
sw make f85 |
sw cdna/make f85 o129 -l f89
