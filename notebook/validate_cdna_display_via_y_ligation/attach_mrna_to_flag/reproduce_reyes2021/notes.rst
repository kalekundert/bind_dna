**********************
Reproduce [Reyes2021]_
**********************

My goal for this experiment is to try to establish a positive control by 
following the [Reyes2021]_ protocol as closely as I can.  My linker is 
different, but hopefully that will be the only difference.

My linker (o236) --- 2021/06/14
===============================
.. protocol:: 20210614_reproduce_reyes2021.pdf 20210614_reproduce_reyes2021.txt

.. figure:: 20210614_reproduce_reyes2021.svg

Observations:

- I don't see any evidence of mRNA/peptide coupling.

- I also don't see any evidence of FLAG expression.

  In :expt:`117`, I was able to visualize FLAG expression with PUREfrex 1.0 
  using SDS PAGE.  But I used TBE/urea PAGE for this experiment because I 
  didn't want to salt to interfere with the gel.  I thought about using 
  desalting columns, but FLAG is small enough that it'd be retained by the 
  column.

  [Reyes2021]_ used SDS/urea PAGE, and just pelleted the SDS precipitate.  In 
  my expreience, though, this causes some bands to not appear.

[Reyes2021] linker
==================
Since I still haven't seen any evidence of mRNA/peptide coupling, I want to try 
reproducing the [Reyes2021] protocol exactly---namely with the same linker.

.. protocol:: 20210810_make_f127.pdf 20210810_make.txt 20210813_dilute.txt

.. protocol:: 20210813_reproduce_reyes2021_o243.pdf 20210813_reproduce_reyes2021_o243.txt

.. figure:: 20210813_reproduce_reyes2021_o243.svg

Observations:

- I can only visualize the mRNA in this experiment, because I'm using the exact 
  same linker as [Reyes2021]_, which is FITC-labeled.

- This doesn't look like Fig 2a from [Reyes2021]_.  That figure shows two 
  distinct (and sharp) bands, representing the reacted and unreacted mRNA.  The 
  sharpness may be attributable to their using SDS-PAGE (although this doesn't 
  work in my hands...).

- The high-salt bands are slightly shifted relative to the mRNA control, but I 
  think this is probably just because of the PUREfrex mixture.

- I've repeatedly seen that the +salt −PUREfrex  control is much dimmer than 
  the corresponding −salt control.  I don't think this is a fluke anymore, but 
  I don't have an explanation.

- I need to do a western blot.  It'll be a lot easier to see the FLAG moving.


