************************************
Express mWasabi-repA with NEBexpress
************************************

- T7 instead of tac.  

  But I could use tac.  The RNAP is added exogenously, and I've already ordered 
  wildtype RNAP.  But I'll have to figure out how much to add.

- Should give higher expression with less template DNA.

- Is an extract; so should have rho factor.  But not clear how rho factor 
  interacts with T7 RNAP.

- First test expression, then make all the templates and see how well 
  Cis-display is working.


Results
=======

2020/03/05
----------

.. protocol:: 20200305_nebex_page_page.txt

.. figure:: 20200305_nebexpress_f17_f22.svg

- I think the DNA template is being degraded:

  - Neither of the NEBExpress lanes have a Cy5 band.

  - There isn't really a band corresponding to the fusion product in the B 
    lanes, but there are lower MW bands that are not present in any of the 
    other reactions.  This could be a consequence of the DNA template being 
    partly degraded before a significant amount of mRNA is made.

  The obvious thing to try doing about this is using the GamS nuclease 
  inhibitor that I didn't have today.
  
- It may also make sense to try using plasmid DNA, which should be more 
  resistant to degradation than linear DNA.  The issue would be finding a way 
  to see the unlabeled DNA, but I have some options:
  
  - I might be able to use GelRed.  GelRed is better excited by the UV lamp on 
    the Gel Doc (300 nm) than GFP is.  GelRed also emits at higher wavelengths 
    than GFP does, although I don't know what filter (if any) the Gel Doc uses.  
    Likely there would still be some cross-talk, but with controls I might be 
    able to see what's going on.

  - I could express the same mWasabi-repA fusions from two plasmids of 
    significantly different size (e.g. 3 kb vs 6 kb).  If the fusion is 
    successfully binding DNA, I should see the mWasabi bands shifted by 
    different amounts in the two reactions.

   - I could purify the protein.DNA complex, digest the protein, then run a 
     gel.  The only way I should see DNA is if it is bound by the protein being 
     purified.

- Unlike with PURExpress, it doesn't seem like any of the expressed 
  mWasabi-repA is getting stuck in the transcription/translation machinery.  Of 
  course, it also doesn't really seem like any full-length mWasabi-repA is 
  being expressed, so I'll need to get better expression (maybe by using 
  plasmid DNA) before reading too much into this.

- I should try purifying again.  My inability to purify mWasabi-repA from 
  PURExpress is a big reason why I suspect that it's getting stuck in the 
  transcription/translation machinery, so maybe I should use purification as a 
  proxy for "well-behaved".

  Note that I did this with Promega S-30 and got nothing.  I'll want to see 
  expression of full-length mWasabi-repA before trying this.
  
- This gel looks almost identical to the one I made using the Promega S30 
  extract, see :expt:`20200108_express_mwasabi_repa_in_s30_lysate`.  This is 
  kinda reassuring, and probably just means that it doesn't matter too much 
  which extract I use.

.. todo::

   Find a way to express full-length protein using NEBExpress/Promega S30 
   without losing the DNA.
   
   - Add GamS.
   - Stop the reaction sooner.
   - Use plasmid DNA.

