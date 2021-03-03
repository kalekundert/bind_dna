*********************
Ligate with RNA oligo
*********************

- If the linker is RNA, then there would be no substrate for RNase H.

  - Until the RT step...

- `RtcB ligase <https://www.neb.com/products/m0458-rtcb-ligase>`_ can ligate 
  RNA with 2',3'-cyclic phosphate to RNA with 5'-OH.  2',3'-cyclic phosphates 
  are the ends produced by HDV-ribozymes, see :expt:`64`.  That would be an 
  interesting way to make use of the ribozyme:

  - Normal nucleotides have 5' phosphates, and thus are would not be substrates 
    for RtcB.

  - HDV produces 2',3'-cyclic phosphates (which are not substrates for most RNA 
    ligases) and 5'-hydroxyls.  These 5'-hydroxyls would be substrates for 
    RtcB, but would simply reform the ribozyme.

  - The linker is totally synthetic, so I could easily order it with a 5'-OH.

  - Thus, the combination of RtcB and HDV-ribozyme would mean that the only 
    ligation that could possibly go forward would be the one between the linker 
    and the correct-length mRNA.

- How would I get the RT primer?

  - I think I can use RNA as an RT primer; that's what happens in nature.

- Can use 2'-O-methyl or 2'-O-methoxy-ethyl modifications to potentially 
  improve stability while still eliminating RNase H activity.
