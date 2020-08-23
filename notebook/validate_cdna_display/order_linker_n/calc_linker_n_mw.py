#!/usr/bin/env python3

# Glen Research provides "M.W." and "F.W." for each phosphoramidite.  "M.W." is 
# the molecular weight of the protected phosphoramidite, and I can reproduce 
# this number myself.  I'm not sure exactly what "F.W." is.  It seems like it 
# should be the weight of the phosphoramidite once it's been deprotected and 
# incorporated into a chain, but this never quite adds up:
# 
# - For most phosphoramidites, the "F.W." is 1 hydrogen mass larger than the 
#   deprotected/incorporated mass that I calculate.
#
# - For puromycin, the "F.W." is 63 Da greater than my deprotected/ 
#   incorporated mass, and I can't think of any way to account for the 
#   difference.
#
# These discrepancies were bothering me, because it seemed most likely that I 
# was overlooking something.  But I'm pretty sure now that I'm right.  Here is 
# the CAS page for deoxyadenosine monophosphate (653-63-4):
#
# https://www.commonchemistry.org/ChemicalDetail.aspx?ref=653=63=4
#
# The formula of this molecule is C10 H14 N5 O6 P, which corresponds to a MW of 
# 331.21 Da.  The "F.W." of dA from Glen Research is 313.21.  The difference is 
# 18 Da, which is the mass of H₂O.  However, if you look at the chemical 
# structure on the above page, you'll see that 3 hydrogens (2 on the 5' PO₄H₂ 
# and 1 on the 3' OH) and 1 oxygen need to be removed to form the repeating 
# monomer unit.  This means that the correct mass for an deprotected/ 
# incorporated deoxyadenosine monomer is 312.20 Da, which is what I calculated.

# I just noticed that the BEX quote specifies both nmol and µg.  They used 
# roughly the same MW as I calculated, which is reassuring.

# CPG: "Controlled pore glass", a solid-phase support for synthesis, used at 3' 
# end.

fw = {
        'dA': 313.21,
        'dT': 304.20,
        'dC': 289.18,
        'dG': 329.21,
        '5-Me-dC': None,
}
atomic_masses = {
        'h':  1.008,
        'c': 12.011,
        'n': 14.007,
        'o': 15.999,
        'p': 30.974,
        'f': 18.998,
        'cl': 35.45,
}

def mw(*atoms):
    all_atoms = sum_atoms(*atoms)
    return sum(all_atoms[k] * atomic_masses[k] for k in all_atoms)

def sum_atoms(*counts):
    all_atoms = {}
    for atoms in counts:
        for k, v in atoms.items():
            all_atoms.setdefault(k, 0)
            all_atoms[k] += v
    return all_atoms

def formula(*counts):
    all_atoms = sum_atoms(*counts)
    atom_counts = [
            f'{k.upper()}{v}'
            for k, v in all_atoms.items()
    ]
    return ' '.join(atom_counts)

def atoms_from_seq(seq):
    atoms = [term_5_oh]
    atoms_from_nt = {
            'A': da,
            'C': dc,
            'T': dt,
            'G': dg,
    }

    for nt in seq:
        atoms.append(atoms_from_nt[nt])

    atoms.append(term_3)
    return sum_atoms(*atoms)


# Protecting groups

dmt = {
        'h': 5+7+7,
        'c': 21,
        'o': 2,
}
mmt = {
        'h': 5+5+7,
        'c': 20,
        'o': 1,
}
nipr = {
        'h': 14,
        'c': 6,
        'n': 1,
        'o': -1,  # N(iPr)₂ is replaced by O when backbone is oxidized after
        'p': 0,   # synthesis is complete.
}
cnet = {
        'h': 4,
        'c': 3,
        'n': 1,
        'o': 0,
        'p': 0,
}
lev = {
        'h': 7,
        'c': 5,
        'n': 0,
        'o': 2,
        'p': 0,
}
bz = {
        # Benzoyl
        'c': 7,
        'h': 5 - 1,
        'o': 1,
}
ibu = {
        # isobutyryl
        'c': 4,
        'h': 7 - 1,
        'o': 1,
}
cocf3 = {
        'c': 2,
        'o': 1,
        'f': 3,
        'h': -1,
}
succinyl = {
        'c': 4,
        'o': 2,
        'h': 4,
}

# dA (Glen 10-1000)

da = {
        'h': 4+7,
        'c': 5+5,
        'n': 5,
        'o': 5,
        'p': 1,
}
da_protect = sum_atoms(da, dmt, nipr, cnet, bz)

# dC (Glen 10-1010)
dc = {
        'h': 11,
        'c': 9,
        'n': 3,
        'o': 6,
        'p': 1,
}
dc_protect = sum_atoms(dc, dmt, nipr, cnet, bz)

# dG (Glen 10-1020)

dg = {
        'h': 4+7,
        'c': 10,
        'n': 5,
        'o': 6,
        'p': 1,
}
dg_protect = sum_atoms(dg, dmt, nipr, cnet, ibu)

# dT (Glen 10-1030)

dt = {
        'h': 5+7,
        'c': 10,
        'n': 2,
        'o': 7,
        'p': 1,
}
dt_protect = sum_atoms(dt, dmt, nipr, cnet)

# 5-Me-dC (Glen 10-1018)

branch = {
        'h': 12+5+7,
        'c': 6+5+5,
        'n': 3,
        'o': 6 + 1,
        'p': 1,
}
branch_protect = sum_atoms(branch, dmt, nipr, cnet, lev)

# Spacer18 (Glen 10-1918)
spacer18 = {
        'c': 12,
        'o': 9,
        'h': 24,
        'p': 1,
}
spacer18_protect = sum_atoms(spacer18, dmt, nipr, cnet)

# iCy5 (Glen 10-5915)
cy5 = {
        'c': 3+10+5+10+3,
        'n': 2,
        'o': 4,
        'h': 6+10+5+10+6,
        'p': 1,
}
cy5_protect = sum_atoms(cy5, mmt, nipr, cnet, {'cl': 1})

# Puromycin (Glen 10-4040, 10-4140)

puromycin = {
        'c': 7+5+10,
        'h': 8+7+13,
        'n': 7,  
        'o': 5,
}
puromycin_protect = sum_atoms(puromycin, dmt, cocf3)

# 5' phosphate
term_5_po4 = {
        'p': 1,
        'o': 3,
        'h': 1,
}
term_5_oh = {
        'h': 1,
}

# 3' OH
term_3 = {
        'p': -1,
        'o': -2,
        'h': 1,
}

# Linker-N
#                            123456789 12345678           123456789 12345678
# Main chain: 5′-(Phosphate)-CCCCCCCGCCGCCCCCCG-(5-Me-dC)-AAAAAAAAAAAAAAAAAA-(Spacer18)-(Spacer18)-(Spacer18)-(Cy5)-(Spacer18)-CC-(Puromycin)-3′
# Side chain: 5′-CCTTG-3′
linker_n = [
        term_5_po4,

        dc,
        dc,
        dc,
        dc,
        dc,
        dc,
        dc,
         dg,
        dc,
        dc,
         dg,
        dc,
        dc,
        dc,
        dc,
        dc,
        dc,
         dg,

        branch,

        da,
        da,
        da,
        da,
        da,
        da,
        da,
        da,
        da,
        da,
        da,
        da,
        da,
        da,
        da,
        da,
        da,
        da,

        spacer18,
        spacer18,
        spacer18,
        cy5,
        spacer18,
        dc,
        dc,
        puromycin,

        dc, 
        dc, 
        dt,
        dt,
        dg,
        term_3,
]
linker_n = sum_atoms(*linker_n)

# Debugging

wiki_da = {
        'c': 10,
        'h': 14,
        'n': 5,
        'o': 6,
        'p': 1,
}
h2o = {
        'h': -2,
        'o': -1,
}
a = sum_atoms(term_5_oh, da, term_3)
aa = sum_atoms(term_5_oh, da, da, term_3)
print('A')
print(formula(a))
print(f"{mw(a):.2f} Da")
print()
print('AA')
print(formula(aa))
print(f"{mw(aa):.2f} Da")
print()

print('Linker-N')
print(formula(linker_n))
print(f"{mw(linker_n):.2f} Da")
