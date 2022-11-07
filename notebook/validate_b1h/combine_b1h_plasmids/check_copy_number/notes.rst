*****************
Check copy number
*****************

The plasmids I got from addgene have 2-3 mutations in the Rep gene (part of the 
SC101 origin): R257H, L274I, H289Y.  I say 2-3 because the "reference" SC101 
sequences I've found aren't all the same.  Mutations in this protein can lead 
to increases in copy number [Thompson2018]_.  None of these mutations are among 
(or even within 100 aa of) those known to affect copy number.

For now, I think it's reasonable to just use the sequence I have.  Presumably 
[Noyes2008]_ used plasmids with these same mutations, and ultimately that's the 
work I'm trying to replicate.  However, it may be worth checking that the copy 
number is still low.  Here's the basic protocol:

- Measure relationship between A600 and cell density
- Grow cells
- Measure OD
- ddPCR/qPCR


