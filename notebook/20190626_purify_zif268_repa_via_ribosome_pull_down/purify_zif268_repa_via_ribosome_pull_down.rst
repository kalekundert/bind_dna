*****************************************
Purify Zif268-repA via ribosome pull-down
*****************************************

See :expt:`20190626_purify_zif268_repa_via_ribosome_pull_down` for why we'd 
like to remove the ribosomes from the IVTT reaction.  Another way to achieve 
this is to use the RiboMinus kit.  The kit contains LNAs with SD sequences, 
which the ribosomes will bind.  These LNAs are attached to magnetic beads, 
which can be pulled down.  Although these kits are relatively expensive, it'd 
be good to at least know if they work.

One nice thing about using magnetic beads to bind the ribosomes is that I could 
add magnetic Ni-NTA beads at the same time, to remove the high-MW PURExpress 
components (e.g. T7 polymerase) for no extra effort.

Considerations
==============

Ribosome concentration
----------------------

- The molecular weight of the E. coli ribosome is 2700 kDa `(Bionumbers 
  100118)`__.  From this I can calculate that 12 pmol is 32.4 μg.  
  
  As a sanity check, this corresponds with the product page for NEB P0763 (the 
  ribosome used in the PURExpress kit), which says that 13.3 μM = 33.3 mg/mL.  
  From this I can calculate that 12 pmol is 30.0 μg.  For future calculations, 
  I think the above value from Bionumbers is more accurate.

  __ https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=10&id=100118

- I can't find anything describing the concentration of ribosomes in the 
  PURExpress kit I've been using.  But in several places---including the 
  `Product Notes`__ page for the PURExpress Δ Ribosome kit (E3313S)---NEB 
  recommends using 60 pmol ribosomes for a 25 μL reaction.  This is 10x more 
  than [Shimizu2001]_ (12 pmol for a 50 μL reaction).

  __ https://international.neb.com/products/e3313-purexpress-delta-ribosome-kit#Product%20Information_Product%20Notes

  Assuming that the regular PURExpress kit uses the same concentration of 
  ribosomes that NEB recommends for the Δ Ribosome kit, my 10 μL PURExpress 
  reactions have 24 pmol (65 μg, 2.4 µM) of ribosomes.

  .. update:: 2020/04/27

     FAQ #7 in NEB's :download:`PURExpress FAQ <purexpress_faq.pdf>` states 
     that the concentration of ribosomes in a standard PURExpress reaction is 2 
     µM ± 20%.  This is pretty close to what I calculated.

- The RiboMinus kit calls for 2-10 μg of total RNA.  Since ~80% of total RNA is 
  rRNA, thats ~2-8 μg of rRNA.  I expect my reactions to have ~8x that much 
  rRNA, so I should expect to need significantly more beads than the protocol 
  calls for.

Magnetic bead value
-------------------
.. datatable:: magnetic_bead_prices.xlsx

The Pierce beads seem to be the best value, especially if I end up buying in 
bulk.  I am trusting the companies self-reported capacity numbers, which I 
don't really trust.  The beads will probably have lower capacity for ribosomes 
than for IgG, but IgG was the largest molecule for which capacity seems to be 
widely reported.


Methods
=======

RiboMinus
---------
.. protocol:: 20190701_purexpress.txt

   - I can use DHFR as a control.  I'd expect to lose it when I do a spin to 
     concentrate, but I should see it in the eluate of the magnetic beads.  
     
   - The lanes I'd run:

      - crude
      - bead eluate
      - 100K flow-through
      - 100K retentate

   - The genes I'd express:

      - DHFR
      - 11

      - Not 11 - ORI.  This construct has given confusing results with the reverse 
        His purification, and I don't think I'd learn anything from repeating it 
        here.

.. protocol::

   Follow RiboMinus protocol with following exceptions:

   - Add 100 μL Ni-NTA magnetic beads (5% solution) along with RiboMinus beads 
     before washing.

   - Mix magetic beads every 3 min while binding ribosomes.

   - Purify protein using 100K spin filters.

A precipitate formed immediately when I added the loading buffer to my samples 
when preparing the SDS-PAGE gel.  Kettner thought that this could be due to the 
presence potassium in my buffers.  Basically, the potassium salt of SDS is much 
less soluble than the sodium salt.  So if there's too much potassium in the 
buffer, the SDS precipitates.
Note that PBS actually includes potassium (2.7 mM KCl, 1.8 mM KH₂PO₄), but 
significantly more sodium (137 mM NaCl, 10 mM Na₂PO₄), which is probably why 
I've been able to run gels with PBS.

I wasn't able to find the composition of hybridization buffer B10 from the 
RiboMinus kit, but it does contain guanidine thiocyanate.  Guanidine (i.e.  
guanidinium) also seems to be able to precipitate SDS is the same manner as 
potassium.  I found a good description of the use of GTC in RNA purification 
protocols from 
[Farrell2010]_:

    Guanidine thiocyanate (GTC) is a stronger protein denaturant than 
    guanidine hydrochloride and is the denaturant of choice for the 
    preparation of RNA from sources enriched in RNase activity, especially 
    pancreatic tissue (Chirgwin et al., 1979). It is routinely used at a 
    working concentration of 4M.

I suspect that GTC is the problem.  Especially if it is 4M, which would 
explain why the ~100x spin-filter dilution wasn't enough to get rid of it.  
Even ignoring its role in precipitating SDS, strong denaturants are also 
incompatible with CIS-display (cDNA display, which is covalent, would work if 
the proteins refold correctly).  So if I want to continue using this protocol 
for CIS-display, I'll need to stop using buffer B10.  That might not be 
possible though---the LNA probes might not be able to reach their binding 
sites without the denaturant.

If I just want to know if the ribosome purification worked, I can just repeat 
the experiment and wash the final retentate several times.  If I want to see 
all the intermediate steps on a gel, I can see 4 options:

 - Use desalting columns.
 - Do drop dialysis.
 - Run a native gel.
 - Don't use buffer B10.
   
Actually, after a brief look, both the Thermo and Biorad desalting columns have 
a 40 kDa MW cutoff.  And I kinda know already that a native gel will be smeary.  
It would probably take several rounds of drop dialysis to get rid of a 4M 
solute, and each round would be pretty tedious.  The pros and cons of using a 
different buffer were discussed above.

.. figure:: 20190702_ribosome_pulldown.svg

- The bands are faint because (presumably) a lot of the protein got caught in 
  the precipitate.  It's also hard to draw conclusions from the absence of a 
  band, because that protein could just be more affected by the precipitate.

- Both DHFR and Zif268-repA are present in the crude reactions.  DHFR can be 
  seen in the bead eluate, but Zif268-repA is cannot.

- Disconcertingly, the ribosomes seem to be eluted from the beads and retained 
  through all the filtering steps.  This is with the GTC buffer (that I want to 
  cut out) and a probable excess of beads (that I want to use fewer of).  If I 
  couldn't even remove the ribosomes in these conditions, it doesn't bode well 
  for this protocol moving forward.  

  .. update:: 2019/07/10

      The beads are not in excess, see Considerations_ above.  This is likely 
      why some ribosomes were retained.

  It is interesting to be that the GTC treatment didn't seem to disassemble to 
  ribosomes, as they were still retained by the 100K spin filter.

Shine-Dalgarno oligos
---------------------
DNA oligos with Shine-Dalgarno (SD) sequences bind to the ribosome with ~30 nM 
affinity [Damian2009]_.  So I might be able to make my own ribosome pulldown 
protocol by simply ordering 5'-biotin-modified oligos and streptavidin-coated 
magnetic beads.  The advantage of this approach is that by targeting intact 
ribosomes, it should not require either denaturing solvents or elevated 
temperatures.

- I want to address the following points:
  
   - Which of the three oligos I ordered works the best?

   - How much beads/oligos should I uses?

   - How many batch purifications should I do?

- I also want to use as little PURExpress as possible.  I expect that doing ~3 
  purifications with a small excess of beads/oligos will work the best.  I 
  think my plan is to first test all the oligos in large excess, then to go 
  from there.
  
Oligos:

- My oligos are 100 μM, i.e. 100 pmol/μL.  If I have 24 pmol ribosomes in my 
  reactions, I'll need at least 0.24 μL of each oligo.

- I probably don't want a super-huge excess of oligos, because unbound oligos 
  could compete with the bound oligos for spots on the beads (even though the 
  unbound ones could fit in a lot of spots that the bound ones couldn't).  An 
  excess will help dive the ribosome binding reaction to completion, though.  I 
  might also need an excess to help out-compete the DNA added to the reaction.

- For my first experiment, I'll use 2.4 μL (10x excess).

Beads:

- 1 μL of the Pierce beads have a `capacity`__ of:
  
   - 551 ng (3.6 pmol) of IgG, a 150 kDa protein.
     
   - 22.6 ng (35 pmol) of biotinylated fluorescein, a 644.71 Da small molecule.

  __ https://assets.thermofisher.com/TFS-Assets/LSG/figures/streptavidin-magnetic-beads.jpg-650.jpg

- To roughly predict the capacity of the beads for intact ribosomes, I'll make 
  the following assumptions:

   - Binding capacity is proportional to the surface area occluded by the 
     target.

   - Surface area is proportional to volume**(2/3).

   - Volume is proportional to mass.

  From this, I calculate that IgG occludes ~40x the surface area of 
  biotinylated fluorescein.  This roughly corresponds (e.g. same order of 
  magnitude) to the 10x difference in bead capacity for these two targets.

  If this relationship holds, an intact ribosome would have ~7x the surface 
  area of IgG.  If this corresponds to a 7x decrease in bead capacity, 1 μL of 
  beads could bind 0.5 pmol intact ribosome.

- My 10 μL PURExpress reactions have 24 pmol of ribosome, so I would expect to 
  need about 50 μL of beads per reaction.

- For my first experiment, I'll use 100 μL.  The calculations above are pretty 
  approximate, so this may or may not be enough to get rid of all the 
  ribosomes, but hopefully it'll be enough to see a difference.

.. protocol:: 20190719_purexpress.txt

   - Setup the IVTT reactions without template DNA.  The template may interfere 
     with oligo binding, and for now I just want to know if this idea could 
     work in the most ideal circumstances.

   - Incubate at 37°C for 5 min (just to warm everything up).

   - Add 2.4 μL 100 μM oligos.

   - Incubate at 37°C for 1h.

   - Wash 50 μL beads in TBST.

   - Dilute ribosomes+oligos to 30 μL with TBST.

   - Add diluted ribosomes to washed beads.

   - Mix at RT for 1h

   - Keep supernatant

   - Run E-gel
      
      - Very smeary. I think the salt or tween is messing with the gel.

      - I tried running a 10x dilution, but the bands were very faint.

   - Nanodrop

.. datatable:: nanodrop.xlsx

   RNA concentrations as measured by nanodrop in "duplex RNA" mode.

- The negative control is probably lower than everything else because it didn't 
  get as much master mix.  That was a real flaw in how I set up the experiment.  
  I should've made excess master mix (rather than making just enough and using 
  whatever is leftover as the negative control) because comparisons with the 
  negative control are the whole point of this experiment.

- Regardless, I can still say that none of the oligos seemed to deplete the 
  ribosomes at all.

- Would be nice to visualize the bead retentate, but I'm not sure how to do 
  this reliably.

Results
=======
I've concluded that purifying the reaction by pulling down the ribosomes is a 
dead end.  The established pulldown methods are too harsh, this gentle method 
doesn't give any indication of working, and both approaches struggle with the 
sheer quantity of ribosomes in the PURExpress reactions.

