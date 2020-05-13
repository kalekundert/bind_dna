#!/usr/bin/env bash
set -euo pipefail

stepwise anneal $@ -c1.25 -C10
