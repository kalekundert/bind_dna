***********************************
Express Zif268 via [Naimudden2016]_
***********************************
I'm modeling my effort to implement cDNA display on two papers: 
[Naimudden2016]_ and [Barendt2013]_.  [Naimudden2016]_ does cDNA display using 
linker-N (as I want to do), but expresses the proteins in rabbit reitc lysate 
and wheat germ extract.  [Barendt2013]_ does mRNA display, but uses PURExpress 
(as I want to do). 

Can a hybrid of these methods give efficient expression mRNA-displayed Zif268?

Considerations
==============

Controls
--------
- No mRNA.

- mRNA without linker.

mRNA purification
-----------------
- [Barendt2013]_ used an Amicon 100 kDa spin filter to remove unligated linker 
  (which contains puromycin) from the mRNA.

   - "a 100 kDa ultrafiltation device can effectively remove the splint and 
     unligated linkers (9.5 kDa) while retaining the mRNA-DNA-puromycin 
     molecules (170–210 kDa) with a yield of 85–90%."

   - My Zif268 (p49) mRNA is 239 kDa, so this approach should work well.  
     Including this step will also let me dilute the mRNA as much as I want in 
     the ligation reaction, which might be necessary to avoid intermolecular 
     ligation.

   - [Barendt2013]_ added urea to the reaction before using the spin filter, 
     then washed with urea and water.  I assume this is to separate oligos that 
     are annealed but not covalently ligated to the mRNA.  This is probably a 
     good idea.

- [Naimudden2016]_ doesn't specify any purification of the mRNA.  However, 
  [Naimudden2011]_ (a earlier version of the same method) specifies that an 
  RNeasy kit is used.

- I think the Amicon filter is a better idea.  I'm not sure about RNeasy, but 
  my Zymo kits will probably retain linker-N.  Plus the urea wash seems like a 
  good idea.

In vitro translation
--------------------
- [Barendt2013]_ used PURExpress:

   - 60 pmol mRNA, 15 pmol ribosomes per reaction.
     
     These reactions are smaller even than the 10 µL reaction I normally do!  
     In :expt:`20190626_purify_zif268_repa_via_ribosome_pull_down`, I work out 
     that a 25 µL PURExpress reaction contains 60 pmol of ribosomes.  Assuming 
     [Barendt2013]_ isn't diluting the reaction, 15 pmol of ribosomes then 
     corresponds to a reaction volume of 6.25 µL.  This is consistent with the 
     5.5 µL of master mix specified in the methods, and implies that the 60 
     pmol of mRNA are added in a volume of 0.75 µL (i.e. 80 µM stock).
     
   - Incubate at 37°C for 30 min, then 25°C for 10 min.

   - No salts added.

- [Naimudden2016]_:

   - Rabbit reticulocyte lysate (RRL):

      - Incubate 10 min at 30°C.

      - Add 65 mM MgCl₂, 750 mM KCl

      - Incubate 2h at 37°C

   - Wheat germ extract (WGE):

      - Incubate 10 min at 25°C

      - Add 65 mM MgCl₂, 750 mM KCl

      - Incubate 1h at 25°C

   - The amount of mRNA to use is not specified.  But [Naimudden2011]_ 
     specifies 3-5 pmol of mRNA in 25 µL RRL.

   - The high-salt incubations are claimed to be important, although I'm not 
     sure why it wouldn't be important for [Barendt2013]_.  Probably I should 
     try both.

- My instinct is to follow the PURExpress protocol more closely here.

   - I only have ~5 pmol RNA today, maybe more if I pool some reactions.  I 
     also don't want to scale down the PURExpress reaction below 5 µL, so 
     hopefully the mRNA:ribosome ratio isn't critically important...

   - I'll have to test if the high-salt incubation is important.  For now, I 
     should just do it.

Results
=======

First try --- 2019/12/20
------------------------
Total failure.  A lot of things could've gone wrong; this was a cowboy 
experiment.  I'll have to step through things more carefully.

.. update:: 2019/12/31

   I noticed that my linker-N doesn't have the reverse transcription primer 
   arm.  This shouldn't have affected this experiment (the puromycin is still 
   there), but I won't be able to progress beyond this step until I get the 
   right linker.
