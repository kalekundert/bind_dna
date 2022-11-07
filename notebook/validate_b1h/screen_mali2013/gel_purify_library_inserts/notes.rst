**************************
Gel purify library inserts
**************************

After amplifying and digesting the library oligos, I think it might make sense 
to gel purify them.  I would need to use PAGE for this, since the oligos are 
very small (and not much bigger than the digested ends).

I noticed that the PAGE gel purification protocols from the `QIAEX II manual`_ 
and the `Novex TBE gel FAQ page`_ have clear similarities.  They also have 
similarities to the PAGE purification protocol that I got from Fitzy:

- "Diffusion/elution" buffers:

  - Novex/QIAEX II:

    - 500 mM ammonium acetate
    - 10 mM magnesium acetate/sulfate (the counter-ion is the only difference)
    - 1 mM EDTA (pH 8.0)
    - 0.1% SDS

  - Fitzy:

    - 10 mM tris, pH 7.5
    - 500 mM NaCl
    - 1 mM EDTA
    - 0.1% SDS

  - The presence/absence of magnesium seems like the biggest difference between 
    these buffers.

- Diffusion/elution protocols:

  - Novex:

    - Crush the gel
    - Incubate at 37°C for 3-4h (for fragments <500 bp)
    - Wash polyacrylamide pellet once

  - QIAEX II:

    - Incubate at 50°C for 30 min.

  - Fitzy:

    - Crush the gel
    - Incubate at 55°C overnight with 800 rpm mixing.

- Remove residual polyacrylamide:

  - Novex/QIAEX II:

    - Use empty spin column
    - I assume something like Biorad 7326204

  - Fitzy:

    - Use cellulose-acetate Spin-X column

  - I like the Spin-X columns more:

    - Smaller pore size, which should give best removal.
    - Cellulose acetate known to have low DNA binding.

  - Probably doesn't really matter, though.

- DNA-recovery protocol:

  - Novex: ethanol precipitation

  - QIAEX II: itself, obviously

  - Fitzy: spin-column

