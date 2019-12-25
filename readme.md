## Scripts and figures accompanying _Diverse Primary and Secondary Structural Features Are Associated With Y Complex-Dependent mRNA Maturation in Bacillus subtilis_

Below you will find detailed instructions for recreating each figure/means of analysis presented in the paper. 

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

Use sequences/rna\_shortened\_21\_cleavage\_final.txt [this is a fasta-format file containing RNA sequences of all 21 windows used--the length of each window is standardized to be the length of the shortest window out of all 21 and all cleavage sites occur between one-indexed nucleotides 11 and 12] as an input to [weblogo](https://weblogo.berkeley.edu/logo.cgi), making sure to click the "Frequency Plot" option.

### Figure 4
![Figure 4a](figures/fig4a.png)
![Figure 4b](figures/fig4b.png)




