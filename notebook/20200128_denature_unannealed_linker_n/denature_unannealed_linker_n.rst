****************************
Denature unannealed linker-N
****************************

I noticed in :expt:`20200127_eliminate_mrna_degradation` that a significant 
amount of linker (o93) appeared to bind the mRNA even without any of the 
ligation reagents or incubation steps.  This suggests that the two oligos 
anneal readily and are not fully denatured by the treatment prior to loading 
the gel, which consists of:

- urea loading buffer
- incubate at 70°C for 3 min

I need a way to fully denature this control, otherwise it will be impossible to 
use gel electrophoresis to monitor the ligation reaction.  Here I'm 
experimenting with ways of achieving this.

Results
=======

Formamide, 10 min incubation --- 2020/01/28
-------------------------------------------
I'm trying a formamide loading buffer (NEB B0363) and a 10 min incubation at 
70°C.  Fitzy told me that formamide is better at denaturing RNA that urea, and 
that formamide loading buffers can be used with urea gels.  I'm also trying a 
prolonged incubation at 70°C, although 10 min at 70°C still falls within the 
range recommended by NEB.

.. figure:: 20200128_denature_controls.svg

- The formamide loading dye significantly reduces the amount of annealed 
  linker.  Formamide might even completely eliminate annealing for linker-N, 
  although this conclusion is complicated by the fact that I used Gel Red 
  (discussed more below) and by the fact that linker-N is less bright in this 
  gel because I used the 520 nm laser (the 488 nm laser/filter matches FITC 
  better).

  I can quantify how much better formamide is than urea (for o93), because the 
  annealed bans in the urea lanes are saturated.  But it appears that the 
  improvement is on the order of 10x.

  I also like that formamide has less of a shadow.

- It would've been nice to take an image with lower intensity, so I could do 
  densiometry on the bound vs unbound lanes.  That might be the best way to 
  quantify ligation moving forward.

  The problem with only looking at the linker is that I'd either have to have 
  the mRNA in excess or have an exact 1:1 mRNA:linker ratio.  I'm already 
  aiming for a 1:1 ratio, but any deviations from that would look just like 
  different ligation efficiencies, and I don't want to depend too much on being 
  able to accurately pipet volumes as small as 0.5 µL.  Furthermore, I'd 
  ultimately rather have the linker is slight excess, since I can get rid of 
  unreacted linker.
  
- I shouldn't have used GelRed, because it makes it hard to see what's 
  happening with linker-N, and I should've used the 488 nm laser.  I can see 
  that formamide reduces the amount of annealed product, because there's some 
  separation between the annealed and un-annealed bands, but I can't say that 
  formamide completely removed the annealed product.

- I don't see a significant difference between incubating at 70°C for 3 min or 
  10 min.

- I think the imager was dirty, because there's a big green blot under most of 
  my gel.  I should remember to wipe it down before taking images in the 
  future.

Discussion
==========
- The formamide loading buffer is much better than the urea buffer for 
  preventing annealing between the mRNA and the linker during gel 
  electrophoresis.

- It may not be possible to completely eliminate such annealing, in which case 
  I'll have to account for it in analysis (e.g. by gel densiometry).
- 
