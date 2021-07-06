#!/usr/bin/env bash
set -euo pipefail
stepwise stain sybr-green-ii/page -I | stepwise laser blue red
