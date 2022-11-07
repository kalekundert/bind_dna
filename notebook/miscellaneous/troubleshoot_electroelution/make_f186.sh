#!/usr/bin/env bash
set -euo pipefail

sw digest p236 NotI-HF,HindIII-HF -n 2 |

sw step "Label the product: f186" |

sw gel agarose/1.5 f186 -l 10 -b TAE -s 'stain gelgreen' |

sw note "The G-CAPSULE protocol notes the TBE can interfere with ligation, so TAE is preferred." "Run a gel" |

sw step "Recover DNA by electroelution.~Follow G-CAPSULE manufacturer's protocol.~Expected band is 2.4 kb."
