# from https://stackoverflow.com/questions/17749058/combine-multiple-text-files-into-one-text-file-using-python
# I had to get windows for each species individually so I'm just combining them here

import glob

read_files = glob.glob("*.txt1")

with open("sequences_for_seqlogo.txt", "wb") as outfile:
    for f in read_files:
        with open(f, "rb") as infile:
            outfile.write(infile.read()) #dump everything in infile to the outfile
