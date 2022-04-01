****************
Use plate reader
****************

Currently, I'm measuring bacterial fitness by plating cells on selective 
plates, the same as [Noyes2008]_ did.  However, it would be easier to use a 
plate reader.  In this experiment, I want to find a good plate reader protocol, 
and test its accuracy.

Note that I can use a plate reader for both the auxotrophic and fluorescent 
reporters, and that the assay would be slightly different in both cases.  Here 
I'm focusing specifically on the auxotrophic reporters, but I expect that many 
of the details will apply to both.

Considerations
==============
I briefly looked for plate reader protocols, and decided that [Kurokawa2017]_ 
seems like a reasonable place to start.

Evaporation
-----------
It's important to prevent evaporation, because OD600 measurements depend on 
path length.  Here are the methods I've found to deal with this:

- Lid with parafilm

  - Specifically, the parafilm is wrapped around the edge of the plate to 
    prevent vapor from escaping.

  - It seems that condensation on the plate could be a problem.  The parafilm 
    might help with that too, though.

  - Important to get lids with "condensation rings" over each well, otherwise 
    condensation from one well can get into another, leading to 
    cross-contamination.

- Mineral oil

  - This is what I've done previously, but I'm sure it impairs oxygen exchange.

I'm going to try lids to start with.

Plates
------
I need plates with the following properties:

- 96 wells
- ≈300 µL wells
- flat bottom: so path length is constant across whole well
- clear bottom: so absorbance measurements can be made
- lids with condensation rings: to minimize cross-contamination
- sterile

One thing I wasn't sure about is what kind of coating, if any, is best:

- Untreated:

  - These are considered "medium-binding".  Polystyrene is hydrophobic, so 
    hydrophobic molecules may adsorb to it.  

- Corning ultra-low attachment

  - Inert hydro-gel covalently bonded to the well surfaces to minimize binding.
  - About 33% more expensive that untreated plates.
  - I get the impression that low-binding surfaces are more important for 
    enzymatic assays.

- Tissue-culture (TC) treated plates

  - Surface made net-negative to support cell attachment [Auld2020]_
  - "The tissue culture treatment process involves exposing a polystyrene 
    microplate to a plasma gas in order to modify the hydrophobic plastic 
    surface to make it more hydrophilic. The resulting surface carries a net 
    negative charge due to the presence of oxygen-containing functional groups 
    such as hydroxyl and carboxyl. In general, this will lead to increased cell 
    attachment." 
    https://www.perkinelmer.com/lab-products-and-services/application-support-knowledgebase/microplates/plate-treatments.html
    
I think that untreated plates are the way to go.  I just don't really think 
surface binding will be a problem.  I also found this quote from Thermo: "The 
non-treated surface is available for suspension cell culture or general 
solution-based assays".

I decided to order Nunc Edge non-treated 96-well plates (Thermo 267427).  These 
have a moat around the edge, which should help with evaporation.  I can't tell 
if the lids have condensation rings or not, but I guess I'll find out.


Analysis
========
I'm trying to figure out the best way to analyze these data.

Fit to exponential function
---------------------------
My thought was to fit just the beginning of the growth curve (e.g. data points 
with OD600 < 0.2) to an exponential curve, given that bacterial growth should 
be an exponential process.  However, the data are just not exponential:

- The growth is completely flat before the culture density exceeds ≈0.1.  
  Presumably this is just because the plate reader is not sensitive enough to 
  detect changes in this regime.

- When I plot the log of the growth values, there isn't any extended linear 
  part of the curve.

Fit to a Gaussian process
-------------------------
Instead of calculating a doubling time, I probably have to calculate the 
maximum slope.  [Kurokawa2017]_ does this by looking at the point in groups of 
five, but I don't like that approach because it seems too ad-hoc and dependent 
on the interval at which measurements are made.  Instead, I want to try fitting 
the data to a Gaussian process, then calculating the maximal slope from that.
