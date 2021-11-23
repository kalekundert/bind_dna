#!/usr/bin/env bash
set -euo pipefail

sw digest p2 HaeII |

sw ligate 'HaeII fragments' -c 100 -l 2686 |

sw e_gel
