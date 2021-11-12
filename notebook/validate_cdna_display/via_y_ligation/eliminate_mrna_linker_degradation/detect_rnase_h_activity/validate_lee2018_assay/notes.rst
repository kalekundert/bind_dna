*************************
Validate [Lee2018]_ assay
*************************

This assay is based on RCA; a similar concept to the qPCR idea I've been trying 
to implement myself.  In this case, the RNA is contained in a primer with a 
blocked 3' end, such that cleavage by RNase H frees the 3' end and allows RCA 
to occur.

Considerations
==============

3' modification
---------------
[Lee2018]_ use a "3'-amino" amino modification to prevent polymerase activity.  
They tried several other 3' modifications, and found that this was the best.  

The way I read the paper, it's strongly implied that this modification is 
literally the 3' hydroxyl group being replaced by a 3' amino group.  However, 
this is not the case.  [Lee2018]_ ordered their oligos from Bioneer, and the 
only 3' amine that Bioneer offers is "3' C6 Amine" modification.  This 
modification has the amine attached to the 3' oxygen via a 6-carbon linker.

IDT does not seem to sell this exact modification.  The closest is the `3' 
amino modification`__, which has an alcohol and an ether in the 6-carbon 
linker.

__ https://www.idtdna.com/Site/Catalog/Modifications/Product/3299

Notably, the 3' C6 amine modification is not meant to prevent polymerase 
extension.  It's meant to provide an attachment site for dyes and so forth.  
The modification that's meant to stop polymerase is inverted-dT.  For that 
reason, I'm tempted to replace the 3'-amino with a 3'-invdT.

[Dames2007]_ actually compares the blocking efficacy of several 3' 
modifications, including 3' C6 amine and inverted-dT.  Both perform well, but 
the C6 amine performs slightly better.  Notably, [Dames2007]_ ordered their 
primers from IDT, so they probably had the same exact chemistry as I'll have, 
although the paper does not actually describe the chemistry (and the IDT 
modification is not exactly a C6 amine).  

All things considered, though, I think the C6-amino modification is the best 
choice.  It should have the best performance and should be what [Lee2018]_ 
used.  The only downside is that I'm not really certain if what I'm ordering is 
the same as what those groups used.  But worst case, I'll just order the 
inverted-dT oligos later.

Ribosomes
---------
The translation reactions have a large amount of RNA in the form of ribosomes.  
Since the ultimate readout is SBYR Green II, these ribosomes might generate a 
lot of signal.  Of course, this signal should not change with time, but I might 
still have problems with the background.  I'll just have to try it and see.

