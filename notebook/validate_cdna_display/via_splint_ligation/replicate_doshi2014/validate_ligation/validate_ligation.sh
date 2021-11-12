#!/usr/bin/env bash
set -euo pipefail

sw zap |

sw cond validate_ligation.xlsx |

sw cdna/splint f145 o278 -d |

sw gel urea/ligate/doshi2014 7
