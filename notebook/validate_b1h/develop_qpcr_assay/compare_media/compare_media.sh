#!/usr/bin/env bash
set -euo pipefail

sw step "Dilute each template 10x."

sw pcr \
  -p ssoadv \
  -d 2 \
  --skip-thermocycler \
  sz224+RT/IPTG=0/Carb=1,sr1,sr2 \
  sz224+RT/IPTG=10/Carb=1,sr1,sr2 \
  sz224+RT/IPTG=100/Carb=1,sr1,sr2 \
  sz224+RT/IPTG=10/Carb=0.5,sr1,sr2 \
  sz228+RT/IPTG=0/Carb=1,sr1,sr2 \
  sz228+RT/IPTG=10/Carb=1,sr1,sr2 \
  sz228+RT/IPTG=100/Carb=1,sr1,sr2 \
  sz228+RT/IPTG=10/Carb=0.5,sr1,sr2 \
  \
  sz224-RT/IPTG=0/Carb=1,sr1,sr2 \
  sz224-RT/IPTG=10/Carb=1,sr1,sr2 \
  sz224-RT/IPTG=100/Carb=1,sr1,sr2 \
  sz224-RT/IPTG=10/Carb=0.5,sr1,sr2 \
  sz228-RT/IPTG=0/Carb=1,sr1,sr2 \
  sz228-RT/IPTG=10/Carb=1,sr1,sr2 \
  sz228-RT/IPTG=100/Carb=1,sr1,sr2 \
  sz228-RT/IPTG=10/Carb=0.5,sr1,sr2 \
  |

sw pcr \
  -p ssoadv \
  -d 2 \
  --skip-thermocycler \
  sz224+RT/IPTG=0/Carb=1,sr3,sr5 \
  sz224+RT/IPTG=10/Carb=1,sr3,sr5 \
  sz224+RT/IPTG=100/Carb=1,sr3,sr5 \
  sz224+RT/IPTG=10/Carb=0.5,sr3,sr5 \
  sz228+RT/IPTG=0/Carb=1,sr3,sr5 \
  sz228+RT/IPTG=10/Carb=1,sr3,sr5 \
  sz228+RT/IPTG=100/Carb=1,sr3,sr5 \
  sz228+RT/IPTG=10/Carb=0.5,sr3,sr5 \
  \
  sz224-RT/IPTG=0/Carb=1,sr3,sr5 \
  sz224-RT/IPTG=10/Carb=1,sr3,sr5 \
  sz224-RT/IPTG=100/Carb=1,sr3,sr5 \
  sz224-RT/IPTG=10/Carb=0.5,sr3,sr5 \
  sz228-RT/IPTG=0/Carb=1,sr3,sr5 \
  sz228-RT/IPTG=10/Carb=1,sr3,sr5 \
  sz228-RT/IPTG=100/Carb=1,sr3,sr5 \
  sz228-RT/IPTG=10/Carb=0.5,sr3,sr5 \
