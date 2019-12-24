# from https://stackoverflow.com/questions/17749058/combine-multiple-text-files-into-one-text-file-using-python
# I had to get windows for each species individually so I'm just combining them here

import glob

read_files = glob.glob("*.fasta")

with open("all_genomes.fasta", "wb") as outfile:
    for f in read_files:
        with open(f, "rb") as infile:
            outfile.write(infile.read())
