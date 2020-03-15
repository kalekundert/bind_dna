#!/usr/bin/env bash
set -euo pipefail

GG="BsaI,BbsI"
klab_reverse_translate his.prot -d his.fa -R $GG
klab_reverse_translate ura.prot -d ura.fa -R $GG
