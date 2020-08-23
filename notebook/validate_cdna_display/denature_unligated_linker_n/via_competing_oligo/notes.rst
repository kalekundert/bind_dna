*******************
Via competing oligo
*******************

Chemical denaturants such as formamide and urea do not seem to be enough to 
fully denature duplex DNA/RNA (see especially :expt:`17`, but also :expt:`49`).  
A different approach is to add a competing oligo that would prevent linker-N 
and the mRNA from re-annealing after being heat-denatured.

- Add free oligo
- Denature by heating.

Considerations
==============

Oligo target
------------
The competing oligo could bind either linker-N or the mRNA.  Both are present 
in the same molar concentrations, so the same amount of competing oligo would 
be required either way. 

Advantages of having the competing oligo bind linker-N:

- DNA/DNA duplexes may be stronger than DNA/RNA duplexes.  (This sounds 
  reasonable, but I have no idea if it's actually true or not.)
  
Advantages of having the competing oligo bind the mRNA:

- The oligo can be longer than linker-N, which would make it thermodynamically 
  favored.

- The Y-tag mismatch might serve as a kind of toe-hold that would help the 
  competing oligo displace linker-N.

- Don't need to use more oligo if I decide to use more linker-N.

I designed the competing oligo to bind the mRNA, but I actually hadn't even 
considered that it could target linker-N at the time.  I think this choice will 
end up being fine, though.

Results
=======

2020/08/14:

.. figure:: 20200814_titrate_o194_saturated.svg

.. datatable:: 20200814_titrate_o194.xlsx

- Both the linker (Cy5) and mRNA (GelGreen) channels agree in terms of how far 
  the ligation reaction progressed, which is reassuring.

- Although the −ligase reactions look very different with and without the 
  competing oligo (indicating that the oligo is preventing unligated binding), 
  the +ligase reactions actually look about the same.

Discussion
==========
- The competing oligo very effectively eliminates non-covalent annealing.

