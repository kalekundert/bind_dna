****************************
Check T4 DNA ligase activity
****************************

I've been having trouble with my cloning reactions recently, so I wanted to 
check whether or not my ligase was compromised.  I found a simple assay 
described by Thermo based on the ligation of HindIII-digested lambda DNA:

2021/07/21
==========
:download:`thermo_ligase_control.pdf`

.. protocol:: 20210721_check_ligase.txt

.. figure:: 20210721_check_ligase_activity.svg

Observations:

- The ligase is clearly active.  In fact, my results look very similar to the 
  example results shown in the Thermo PDF.

2022/04/13
==========
I decided to repeat this reaction, because I'm still having troubles getting a 
lot of colonies, and now my buffer is past its expiration date.

.. protocol:: 20220413_check_ligase.pdf 20220413_check_ligase.txt 

.. figure:: 20220413_check_ligase_buffer.svg

  ligase: whether or not ligase was included in the reaction.  buffer: whether 
  or not buffer was included in the reaction, and if so, the expiration date of 
  said buffer.

- This time I used the same ligase for both reactions, and only varied the 
  buffer.

- The enzyme seems fully active in both cases, so I guess there's no problem 
  with the buffer.

- I'll probably get rid of the old buffer anyways, since it's expired and I 
  have more of the new buffer than I'll ever use.

- I'm starting to think that my competent cells just aren't very good.
