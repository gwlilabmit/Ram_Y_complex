## Scripts and figures accompanying _Diverse Primary and Secondary Structural Features Are Associated With Y Complex-Dependent mRNA Maturation in Bacillus subtilis_

Below you will find detailed instructions for recreating each figure/means of analysis presented in the paper. If you have questions that are beyond the scope of the paper or this repository [e.g. technical questions like "What shell commands and scripts would I need to run to grab sequence windows for all 21 sites?"] please contact me. 

Please note that all figures used in the paper can also be found under figures/. 

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
You can verify these sequences in Mochiview if you'd like--go to mochiview/subtilis\_mochiview. There are wigs for [all using derivatives of _B. subtilis_ strain 168 in LB exponential] wt, ∆RNase Y, ∆_ylbF_, ∆_yaaT_ and ∆_ymcA_. You can use Mochiview's find function to search for the 200nt windows, i.e. 21\_cleavage\_final.txt. 

### Figure 3
![Figure 3](figures/fig3.png)

Use sequences/rna\_21\_cleavage\_for\_seqlogo.txt [This is a fasta-format file containing RNA sequences of all 21 windows used. I've shortened my windows from the original ~200nt for visualization purposes so each window contains 25nt after the cleavage site. Since not all windows contain 25nt before the cleavage site I've standardized the length of all windows to be that of the shortest window, meaning that there are 12nt before the cleavage site. Please see methods if the rationale behind this is confusing. All cleavage sites occur between one-indexed nucleotides 11 and 12] as an input to [weblogo](https://weblogo.berkeley.edu/logo.cgi), making sure to click the "Frequency Plot" option. You can view the full 112nt window used for frequency plot analysis using rna\_shortened\_21\_cleavage\_final.txt. 

### Figure 4
![Figure 4a](figures/fig4a.png)
![Figure 4b](figures/fig4b.png)

Run staph/values\_from\_rend\_seq\_staph.py. I'm only plotting the 5 seemingly Y complex-dependent sites. The sixth one corresponds to the sequence TACTTACTAAATTTTATTTAACCTAAAAATGAACCACCTGGATGTGTGGG and doesn't seem to be Y complex-dependent. The locations of all staph cleavage sites from [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4608709/) can be found in staph/staph\_rny\_sites.txt. The wig files under staph/ need slight modification to run in mochiview. If this is your goal, go to mochiview/staph\_mochiview/. There you can find wt and ∆rny data for our staph strain [wt files are under wt\* and ∆ylbF are under ylbF\*]. We map reads to the [NC\_007795 genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_007795.1). You can grab the NC\_007795 cds from the NCBI site [here](https://www.ncbi.nlm.nih.gov/genome/proteins/154?genome_assembly_id=299272). This requires some modification to be converted to the Mochiview format, so you can use staph\_mochiview/convert\_mochiview\_location.py to do the conversion for you. Or you can just use the already-converted CDS I have, namely NC\_007795\_mochi\_cds.txt. The one benefit of using the script is that it should work for any CDS that needs conversion from NCBI-\>Mochiview, not just NC\_007795. 

### Figure 5
![Figure 5](figures/fig5.png)

I use [Clustalw2](http://www.clustal.org/omega/#Download) in shell to run all my MSAs. Under msa\_for\_paper I have directories for each cleavage site producing a MSA with alignment-to-subtilis score \>30. All genomes are available under genome\_fasta\_files [this also has a file for the staph genome]. I first grab sequences from my fasta files [atpi\_correct\_evolutionary.txt], standardize the length of all sequences to be that of the shortest sequences [equal\_size\_sequences\_seqlogo.py] and convert all T's to U's [dna\_to\_rna\_sequences.py]. From there I run the MSA on my equal-length RNA sequences [rna\_atpi\_correct\_evolutionary.txt] using the slow/accurate option. If you want to further shorten sequences to visualize their sequence logos I've included shorten\_for\_seqlogo.py--this produces 50nt windows that you can use in weblogo. That's really more for fun, though. Those sequence logos don't really tell you much. 
Unfortunately there's no way to denote the cleavage site in the MSA itself so you have to search for it manually. 

![Figure 14](figures/fig14.png)

If you want to recalculate the alignment-to-subtilis scores, take the average of all pairwise alignment scores of the form "Sequences (1:n) Aligned. Score: " [as seen above]. Note that ClustalW2 will output these scores _prior_ to producing the MSA, so be on the lookout. Sequence 1 in all \*\_evolutionary.txt files is always the _subtilis_ sequence. This was easy enough to do manually that I didn't bother automating it. 




