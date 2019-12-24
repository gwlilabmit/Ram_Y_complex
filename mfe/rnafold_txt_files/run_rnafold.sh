#!/bin/zsh

for filename in *_downstream_hairpins.txt; do 
	if [[ $filename =~ ' re=(^[a-zA-Z0-9]*)' ]]; then
		re=$match[1]
		echo $re
	fi
	#rnafold --noPS --infile=$filename --outfile=$name_hairpins_rnafold.txt
done 
