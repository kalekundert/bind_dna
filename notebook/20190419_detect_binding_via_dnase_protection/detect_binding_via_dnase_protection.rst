***********************************
Detect binding via DNase protection
***********************************

Kettner had the idea that I might be able to detect binding by encoding the 
protein and the target sequence on the same DNA molecule.  Protein binding to 
its target sequence would protect the DNA from DNase treatment (possibly in 
proportion to the strength of binding).  With the right design of barcodes and 
PCR primers, I could then amplify only the sequences that survive DNase 
treatment and identify those protein/target pairs via sequencing.

Advantages:

- The reaction would be unimolecular.  This would allow me to really dilute the 
  reaction and not have to worry about interactions between different library 
  members.
  
- I could tune the sensitivity of the assay by adding free target DNA.

- This assay should be compatible with any DNA display technology.

   - CIS display has an advantage in that I won't have to worry about 
     exonuclease activity from the other direction.
     
   - But I can accomodate cDNA display as well by using the right nuclease with 
     the right phosphorylation (e.g. T7 exonuclease).

Disadvantages:

- Indirect measurement of binding.  
  
   - Factors like the orientation of the DNA binding protein or the length of 
     the tether could be important (although this could be true in my ligation 
     assay, as well).

   - Weak binders may be hard to detect, if they are just too transient to 
     resist the DNase.

- Can't control the effective concentration of protein.  And the effective 
  concentration is unnaturally high because the protein is being confined in 
  such close proximity.

   - Well, I could increase the length of the tether...
