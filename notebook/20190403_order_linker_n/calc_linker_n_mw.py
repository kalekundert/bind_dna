#!/usr/bin/env python3

branch = {
        'h': 12+5+7,
        'c': 6+5+5,
        'n': 3,
        'o': 6 + 1,  # backbone is oxidized to PO₄ before release.
        'p': 1,
}
dmt = {
        'h': 6 + 4 + 4 + 5,
        'c': 21,
        'n': 0,
        'o': 2,
        'p': 0,
}
lev = {
        'h': 7,
        'c': 5,
        'n': 0,
        'o': 2,
        'p': 0,
}
nipr = {
        'h': 14,
        'c': 6,
        'n': 1,
        'o': 0,
        'p': 0,
}
cnet = {
        'h': 4,
        'c': 3,
        'n': 1,
        'o': 0,
        'p': 0,
}

branch_prot = {}
for atoms in [branch,dmt, lev,nipr,cnet]:
    for k, v in atoms.items():
        branch_prot.setdefault(k, 0)
        branch_prot[k] += v

atomic_masses = {
        'h':  1.008,
        'c': 12.011,
        'n': 14.007,
        'o': 15.999,
        'p': 30.974,
}

def mw(atoms):
    return sum(atoms[k] * atomic_masses[k] for k in atoms)

print(f"{mw(branch):.3f} Da")
print(f"Δ = {402.36-mw(branch):.3f} Da")
#print(f"{mw(branch_prot):.3f} Da")

# I think the extra 1 Da is because they're assuming the levulinyl will be 
# replaced by an H, not so another nucleotide.
