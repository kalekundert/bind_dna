************************
Optimize linker ligation
************************

Although the annealing and ligation reactions seem to work with the poly-A 
linker (o129), only about half of the mRNA is modified: :expt:`17`.  I think 
it's probably possible to do better than that, since [Naimudden2016]_ seem to 
be closer to 90% (although with the caveat that they might be confusing FITC 
and SYBR green).

Considerations
==============
After comparing a number of published T4 RNA ligase protocols (see below), I 
think there are several things to think about trying:

- The 10 min incubation recommended by [Naimudden2016]_ is an outlier, although 
  I understand the desire to be gentle with RNA.  Still, I should try longer 
  incubation times.

  This is from the NEB FAQ for `T4 RNA ligase I 
  <https://international.neb.com/faqs/2018/01/30/what-is-the-optimal-reaction-temperature-and-time-for-t4-rna-ligase-i>`__:

    In general, increased reaction time and lowered reaction temperature yield 
    more complete ligation reactions. A typical ligation reaction should be 
    carried out at 25°C for up to 2 hours. For longer oligos, 2 hours of 
    incubation at 25°C followed by overnight incubation at 16°C may improve 
    yield.

- The protocols really differ in the amount of enzyme used:

  - KBK: 4 U/pmol
  - [Naimudden2016]_: 0.4 U/pmol
  - [Kubo2020]_: Units not specified
  - Takara: 2.5–5.0 U/pmol
  - NEB: 0.5 U/pmol
  - [Kitamura2002]_: 5 U/pmol

  That said, I'm already on the high end, so I don't think that adding more 
  ligase is likely to help much.

- All of the non-cDNA-display protocols have 25% PEG, which is known to 
  dramatically improve ligase activity [Bauer2017]_.

- NEB recommends using a 2–10x excess of linker.  That's in line with what I 
  was thinking about trying, anyways.

- Most of the protocols call for adding ATP separately from the buffer.  It's 
  possible that this is important to keep the ATP active, but I do a good job 
  taking care of my T4 DNA ligase buffer, so I don't think I need to worry 
  about this.

My plan (in order):

- Try +/- PEG, since I'm pretty sure that will help.

- Do a timecourse, because I think longer will help, but I also want to keep 
  the incubation as short as possible.  Do every timepoint at 25°C and 16°C, 
  since those are the consensus temperatures.  Time points:

  - 10 min (what I'm currently doing)
  - 2h (NEB recommendation)
  - 16h (most common recommendation)

  How to quench?  I've got several options:

  - Incubate at 65°C.  Don't know where this idea came from though, so it might 
    not actually work.
  - Spin column (as recommended by NEB).  This might get rid of nligated 
    linker, though, which might make it harder to quantify yield.
  - Add EDTA (as recommended by Takara).
  - Formamide loading buffer.  This is assuming that formamide would 
    effectively denature the ligase, which may or may not be the case.
  - Freeze

  I could do a 0h control to make sure the quenching works.  I could also just 
  run three gels...

- Titrate the linker.


:expt:`17`
~~~~~~~~~~
- Setup the ligation reaction:

  - 5 pmol mRNA
  - 5 pmol linker
  - 0.1x PBS (left over from annealing reaction)
  - 1x T4 DNA ligase buffer
    - 50 mM Tris-HCl, pH 7.5
    - 10 mM MgCl₂
    - 10 mM DTT
    - 1 mM ATP
  - 0.01% BSA
  - 20 U T4 RNA ligase
  - water to 40 µL

- Incubate at 25°C for 10 min, then 65°C for 10 min

I'm not sure where exactly this protocol came from.  Some thoughts:

- 10x diluted PBS: my own optimizations of the annealing reaction.

- 20 U T4 RNA ligase: Probably [Naimudden2016]_, even though I'm using 10x less 
  mRNA/linker (I probably assumed that more enzyme wouldn't hurt, and this is 
  already about as little as I can pipet).

- BSA: Probably because Takara included BSA with the enzyme, although I'm using 
  a different concentration that their online protocol suggests.

- 10 min incubation time at 25°C: [Naimudden2016]_.  I don't know where the 
  65°C incubation/heat denaturation came from, though.

[Naimudden2016]_
~~~~~~~~~~~~~~~~
- Setup the ligation reaction:

  - 50 pmol mRNA
  - 50 pmol linker-N
  - 3 U T4 PNK (NEB)
  - 20 U T4 RNA ligase (Takara)
  - Unspecified total volume, but assuming 50 µL based on [Kubo2020]_.
  - Unspecified buffer (possibly no buffer)

- Incubate at 25°C for 10, 20, or 40 min.

[Kubo2020]_
~~~~~~~~~~~
- Setup ligation reaction:
  - 40 pmol mRNA
  - 40 pmol linker
  - 1x T4 DNA ligase buffer
  - 2 µL T4 PNK (unspecified concentration)
  - 2 µL T4 RNA ligase (unspecified concentration)
  - water to 40 µL

- Incubate at 25°C for 15 min (or at 16°C for 2h)

Takara --- T4 RNA Ligase (2050A)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
`Original protocol 
<https://www.takarabio.com/assets/documents/User%20Manual/2050A_DS.v1901Da.pdf>`__.

- Setup the ligation reaction:

  - 1-2 µg ssRNA (8-16 pmol f11)
  - 1x T4 RNA Ligase buffer (same as NEB T4 DNA ligase buffer)
    - 50 mM Tris-HCl, pH 7.5
    - 10 mM MgCl₂
    - 10 mM DTT
    - 1 mM ATP
  - 0.006% BSA
  - 25% PEG 6000
  - 40-50 U T4 RNA ligase
  - water to 50 µL

- Incubate at 5-16°C for 16-18h.

- Stop the reaction by adding 2 µL 500 mM EDTA.

NEB --- T4 RNA Ligase I (M0204)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
`Original protocol 
<https://www.neb.com/protocols/2018/10/17/protocol=ligation=of=an=oligo=to=the=3=end=of=rna=using=t4=rna=ligase=1m0204>`__.  
I'm using T4 RNA ligase from Takara, but this protocol may still be relevant.

- Setup the ligation reaction:

  - 20 pmol RNA
  - 40-200 pmol DNA or RNA oligo
  - 1x T4 RNA Ligase Reaction Buffer (T4 DNA ligase buffer w/o ATP)
    - 50 mM Tris-HCl, pH 7.5
    - 10 mM MgCl₂
    - 1 mM DTT
  - 1 mM ATP (this makes the buffer equivalent to T4 DNA ligase buffer)
  - 10% DMSO (optional)
  - 15-25% PEG 8000
  - 20 U Murine RNase inhibitor (M0314, optional)
  - 10 U T4 RNA Ligase 1
  - water to 20 µL

- Incubate at 25°C for 2 hours or at 16°C for 16 hours

- Stop the reaction with a spin-column cleanup.

[Kitamura2002]_
~~~~~~~~~~~~~~~
- Setup the ligation reaction:

  - 10 pmol 5' end
  - 10 pmol 3' end
  - 50 mM Tris, pH 8.0
  - 10 mM MgCl₂
  - 0.1 mM ATP
  - 0.001% BSA (10 mg/L)
  - 1 mM hexamminecobalt (III) chloride
  - 25% PEG 6000
  - 50 U T4 RNA ligase
  - water to 10 µL

- Incubate at 25°C for 16h.
