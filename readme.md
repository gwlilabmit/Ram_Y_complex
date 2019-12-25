## Scripts and figures accompanying _Diverse Primary and Secondary Structural Features Are Associated With Y Complex-Dependent mRNA Maturation in Bacillus subtilis_

Below you will find detailed instructions for recreating each figure/means of analysis presented in the paper. If you have questions that are beyond the scope of the paper or this repository please contact me. 

### Figure 1
![Figure 1](figures/fig1.png)

Should you want rend-seq visualizations of the 200nt windows I've been using, look no further than rend\_seq/. I've included the following wigs [all libraries prepared from cells in LB exponential]: 
* WT _B. subtilis_ 168 ["wt\*"]
* ∆RNase J1 168 ["rnj\*"]
* 168 treated with a 5' exonuclease to identify 5' monophosphates [to determine if a 5' end is an alternative isoform or the product of cleavage] ["5exo\*"]
* ∆PNPase 168 [to identify 3' peaks that may result from cleavage] ["pnpA\*"]
* ∆RNase Y 168 ["rny\*"]
* ∆_ylbF_ 168 ["ylbF\*"]

Add the appropriate 5'f/r and 3'f/r wigs into plot\_rend\_seq.py [currently can only plot one set of 5'/3/ files but it's a really easy fix to modify if needed] and run!
If you want to visualize a different set of windows, just use a different bedfile as input--21\_cleavage\_final.bed12 contains all 21 windows that I use for analysis. 

### Figure 3
![Figure 3](figures/fig3.png)

Use sequences/rna\_21\_cleavage\_for\_seqlogo.txt [This is a fasta-format file containing RNA sequences of all 21 windows used. I've shortened my windows from the original ~200nt for visualization purposes so each window contains 25nt after the cleavage site. Since not all windows contain 25nt before the cleavage site I've standardized the length of all windows to be that of the shortest window, meaning that there are 12nt before the cleavage site. Please see methods if the rationale behind this is confusing. All cleavage sites occur between one-indexed nucleotides 11 and 12] as an input to [weblogo](https://weblogo.berkeley.edu/logo.cgi), making sure to click the "Frequency Plot" option. You can view the full 112nt window used for frequency plot analysis using rna\_shortened\_21\_cleavage\_final.txt. 

### Figure 4
![Figure 4a](figures/fig4a.png)
![Figure 4b](figures/fig4b.png)




