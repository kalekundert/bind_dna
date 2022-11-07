#!/usr/bin/env bash
set -euo pipefail

# - In :expt:`165`, my Cq values (for the +RT samples) were in the range 11-15.
# - I probably want to try going all the way to the limit of detection, e.g. 35 
#   or so.
# - Each increase of 1 Cq unit corresponds to a 2x decrease in template.  So an 
#   increase of 20 Cq units corresponds to a 2**20 = 1e6x dilution.
# - So I should do 6 10x dilution steps.
sw serial 1 / 10 -n 6 -v 20µL -m cDNA -d water |

sw pcr \
  -p ssoadv \
  -d 2 \
  --skip-thermocycler \
  sz224+RT/1e0,s1,s2 \
  sz224+RT/1e1,s1,s2 \
  sz224+RT/1e2,s1,s2 \
  sz224+RT/1e3,s1,s2 \
  sz224+RT/1e4,s1,s2 \
  sz224+RT/1e5,s1,s2 \
  \
  sz224-RT/1e0,s1,s2 \
  sz224-RT/1e1,s1,s2 \
  sz224-RT/1e2,s1,s2 \
  sz224-RT/1e3,s1,s2 \
  sz224-RT/1e4,s1,s2 \
  sz224-RT/1e5,s1,s2 \
  \
  sz228+RT/1e0,s1,s2 \
  sz228+RT/1e1,s1,s2 \
  sz228+RT/1e2,s1,s2 \
  sz228+RT/1e3,s1,s2 \
  sz228+RT/1e4,s1,s2 \
  sz228+RT/1e5,s1,s2 \
  \
  sz228-RT/1e0,s1,s2 \
  sz228-RT/1e1,s1,s2 \
  sz228-RT/1e2,s1,s2 \
  sz228-RT/1e3,s1,s2 \
  sz228-RT/1e4,s1,s2 \
  sz228-RT/1e5,s1,s2 \
  |

sw pcr \
  -p ssoadv \
  -d 2 \
  sz224+RT/1e0,s3,s5 \
  sz224+RT/1e1,s3,s5 \
  sz224+RT/1e2,s3,s5 \
  sz224+RT/1e3,s3,s5 \
  sz224+RT/1e4,s3,s5 \
  sz224+RT/1e5,s3,s5 \
  \
  sz224-RT/1e0,s3,s5 \
  sz224-RT/1e1,s3,s5 \
  sz224-RT/1e2,s3,s5 \
  sz224-RT/1e3,s3,s5 \
  sz224-RT/1e4,s3,s5 \
  sz224-RT/1e5,s3,s5 \
  \
  sz228+RT/1e0,s3,s5 \
  sz228+RT/1e1,s3,s5 \
  sz228+RT/1e2,s3,s5 \
  sz228+RT/1e3,s3,s5 \
  sz228+RT/1e4,s3,s5 \
  sz228+RT/1e5,s3,s5 \
  \
  sz228-RT/1e0,s3,s5 \
  sz228-RT/1e1,s3,s5 \
  sz228-RT/1e2,s3,s5 \
  sz228-RT/1e3,s3,s5 \
  sz228-RT/1e4,s3,s5 \
  sz228-RT/1e5,s3,s5 \
