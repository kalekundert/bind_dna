***********************************
Express Zif268 via [Naimudden2016]_
***********************************
I'm modeling my effort to implement cDNA display on two papers: 
[Naimudden2016]_ and [Barendt2013]_.  [Naimudden2016]_ does cDNA display using 
linker-N (as I want to do), but expresses the proteins in rabbit reticulocyte 
lysate and wheat germ extract.  [Barendt2013]_ does mRNA display, but uses 
PURExpress (as I want to do). 

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

   .. update:: 2020/02/18

      Formamide might work even better, based on my loading buffer experiments: 
      :expt:`20200128_denature_unannealed_linker_n`.

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

.. figure:: 20191220_cdna_display.svg

- I didn't have urea to use for the ligation cleanup.

- I lost the −mRNA, −linker-N reaction due to SDS precipitation.

- I had to run the remaining reactions through a spin column an extra time to 
  remove KCl from the buffer.  I might have lost things.

- I don't know what the FITC signal at the bottom of the ladder is.

- The signal in the two PURExpress lanes is probably in the green channel only 
  due to crosstalk, and the FITC channel being brighter overall.  The signal is 
  also present in the red channel, although it's hard to see with the green 
  channel superimposed.  Still, it's bizarre that the proteins in the ladder 
  are clearly brighter in the red channel, while the PURExpress proteins are 
  not.

.. update:: 2019/12/31

   I noticed that my linker-N doesn't have the reverse transcription primer 
   arm.  This shouldn't have affected this experiment (the puromycin is still 
   there), but I won't be able to progress beyond this step until I get the 
   right linker.  See :expt:`20190403_order_linker_n` for more about this.

Second try --- 2020/02/18
-------------------------
.. protocol:: 20200218_anneal_ligate_wash_barendt_purex_page.txt

.. figure:: 20200218_express_f11_o93_10_10.svg

   "Annealing steps": Steps 1--4 in the above protocol.  "Filtration steps": 
   Step 5 in the above protocol.  "Expression steps": Steps 6--7 in the above 
   protocol.  Ladder: SeeBlue™ Plus2 Pre-stained Protein Standard.

- The mRNA seems to be badly degraded in the PURExpress reaction.  The mRNA is 
  clearly very present after the filtration step.  That eluate was diluted 
  about 10x into the PURExpress reaction, but should still be easily visible.  
  There's clearly no band corresponding to the full-length mRNA in the 
  PURExpress reaction, though.
  
  Instead, there is a low-MW smear not present in the −pseudo-linker or −mRNA 
  PURExpress reactions.  This must be the mRNA, but I don't know why it's so 
  small.  It doesn't have puromycin, so it shouldn't be reacting with its 
  protein product (although it's kinda near that MW).  I did mark down that I 
  added RNase inhibitor to the reaction.  That said, it's clear the mRNA is 
  breaking down even before the PURExpress reaction, so maybe the 37°C 
  incubation just accelerated this.  

  I think the root problem is that the mRNA is getting degraded, and that's why 
  I'm not seeing any protein.  It's actually kinda nice that I have Cy5-labeled 
  mRNA, because it makes it easier to see what's going on.  

- I don't see a band for Zif268 in the +mRNA +expression lanes.  I can think of 
  two explanations:
  
  - *Zif268 isn't being expressed.*  I haven't done PURExpress directly from 
    mRNA before, and I'm using much less than the recommended amount of mRNA.  
    So it could be that I just don't have much protein.  See 
    :expt:`20200224_express_mrna_via_purexpress`.
  
  - *Zif268 is being obscured by another protein.*  There are bands in the 
    PURExpress reaction at about the MW I'd expect for Zif268, so maybe it's 
    there and just not very highly expressed. 

    It's worth noting that this gel is much lower resolution than I was hoping 
    for, even though I ran this gel in the same way as I have previously.  See 
    :expt:`20190625_purify_zif268_repa_via_reverse_his`, for example.  The 
    exact volumes I loaded onto the gel for that experiment are given in the 
    binder, and I confirmed that they are the same as what I used here.  Here 
    the lanes actually seem overloaded.  Maybe the volumes I listed previously 
    were only for the more diluted elution/wash fractions, and I used less in 
    the crude fractions?  I don't know.  I should try using less in any case, 
    because SYPRO Orange does seem more sensitive than Coomassie.

    .. note::

       Actually, :expt:`20190603_express_zif268_in_vitro` is where I work out 
       how to run these PURExpress reaction on SDS-PAGE gels, and there I used 
       different (and lower) volumes.  Specifically I used 2.5 µL of each 
       reaction per lane, as recommended by NEB.  That's probably what I need 
       to do here.
    
- The filtration steps do help remove low-MW species, although most of the 
  pseudo-linker is reacted anyways.  It's interesting that you can easily see 
  the BSA and the T4 RNA ligase in the unfiltered reaction, and that both 
  proteins are depleted by the filtration steps.  I'm not sure what the >100 
  kDa protein bands are, though.

- I seemed to have a significant amount of RNA degradation in this reaction.  
  Whether or not I have degradation has seemed random.  I didn't think I did a 
  bad job of being careful with this reaction, but maybe I'm just not being 
  careful about the right things.  It might be smart in the future to include a 
  raw mRNA control.

  It's also possible that the mRNA is just not fully denatured by the SDS gel, 
  and this is just what maritally folded mRNA looks like.  I wonder what would 
  happen if I ran my proteins in an TB/urea gel.  It also seems like running 
  SDS-PAGE with up to 6M urea is a thing [Schagger2006]_:

    For unknown reasons, urea reduces the electrophoretic mobility of proteins 
    in general, but the migration of small proteins in particular. Therefore, 
    the resolution of proteins in the low mass range is improved at the cost of 
    a lower resolution for larger proteins

- It's very unexpected that the ladder is present in the red channel (Cy5) and 
  not the green channel (SYPRO Orange).  In fact, the green channel actually 
  has shadows where the ladder bands should be.  Because some of the ladder 
  bands are not easily visible in either channel, I verified the MW assignments 
  by comparing to old gels with the same ladder and gel percentage.
  
  I think the signal in the red channel is most likely due to the near-IR 
  fluorescence of Coomassie [Butt2013]_.  The ladder is prestained, which 
  explains why those bands would have Coomassie.  Note that the 17 and 98 kDa 
  bands are prestained with purple and orange dyes, not Coomassie, which 
  explains why they are much less visible in this channel.  In the future, I 
  could avoid this (if I want to) by using an unstained protein standard, e.g.  
  Novex™ Sharp Unstained Protein Standard (Invitrogen LC5801).
  
  I can think of two reasons why the ladder might not appear in the green 
  channel, but I'm not sure which (if either) is right:
  
  - *Maybe there's no SDS for SYPRO Orange to bind.*  SYPRO Orange binds SDS; 
    it doesn't bind protein directly.  It may be that the proteins in the 
    ladder are denatured by some means other than SDS.  There is still SDS in 
    the running buffer (and the gel, I think), but maybe it's not enough 
    without SDS in the loading buffer.  

    This hypothesis is not specific to Coomassie, which supports the 
    observations that even the non-Coomassie bands (17 and 98 kDa) do not 
    appear in the green channel.

  - *Maybe Coomassie is FRET-ing with SYPRO Orange.*  The emission maximum for 
    SYPRO Orange is 586 nm, see :download:`sypro_orange_ex_em.csv`  (This comes 
    from the `Thermo SpectraViewer 
    <https://www.thermofisher.com/us/en/home/life-science/cell-analysis/labeling-chemistry/fluorescence-spectraviewer.html>`_; 
    choose "SYPRO Orange protein gel stain" and then click "Export".) According 
    to :download:`this pamphlet <dnr_coomassie_blue.pdf>`, the absorption 
    maximum for Coomassie (when bound to protein) is 595 nm, so this is 
    certainly a possibility. 
    
    It's too bad I can't use the 488 nm laser with the 710BP40 filter, because 
    that would make it very clear if FRET was happening.  But this hypothesis 
    is supported by the shadows that I see in the green channel.

- SYPRO Orange is very faintly visible in the 658 nm channel, as can be seen in 
  the +expression lanes.  I don't think this will be a problem in practice, 
  because the true Cy5 bands are very significantly brighter.  I don't plan to 
  be using trace amounts of Cy5, so I don't think I'll have problems seeing it 
  or distinguishing it from SYPRO Orange cross-talk.
