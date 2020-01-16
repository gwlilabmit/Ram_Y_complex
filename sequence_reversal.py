from string_reversal import * 
from base_complement import * 

print "This script will only slice windows of a user-specified size in the genome file provided--if you need to reverse complement or anything like that, run that before running this script."
print "For now it'll get you windows symmetric around the cleavage site. I'll build in more options to get different types of windows as needed later."

#middle = int(raw_input("At what index [you should be able to find this on Mochiview] is the cleavage site? \n"))

###I'm getting some weird index offsets between what I'm seeing in mochiview and the file itself...bear with me
genetype = raw_input("What cleavage site is this? \n")
genetype = genetype.lower()

#bear with me here
#finding the index in NC_000964 for each cleavage site
 
if 'cggr' in genetype:
	middle = 3482769
	offset = 49756 #why the weird offsets? I'm basing this off what I see in Mochiview...if I don't do this my window is not at the indices I'd expect
	middle += offset
elif 'atpi' in genetype: 
	middle = 3783696
	offset = 49756+8237+127 
	middle += offset

elif 'cwlo' in genetype: 
	middle = 3575886
	offset = 49756+8237+127-6925-98	
	middle += offset

elif 'daca' in genetype: 
	middle = 17476
	offset = 263
	middle += offset

elif 'ddl' in genetype: 
	middle = 508162
	offset = 7169+103
	middle += offset

elif 'glna' in genetype: 
	middle = 1878263
	offset = 2650+23853+341
	middle += offset

elif 'mena' in genetype: 
	middle = 3951710
	offset = 55669+796
	middle += offset

#defining window size 
windowsize = int(raw_input("How big should be the window be?\n"))

if windowsize%2==1:
	print "Your window size is odd so I'm going to assume you want two even-sized windows flanking the cleavage site."
	windowsize -=1
else: 
	print "Your window size is even so I'm going to assume you want two even-sized windows flanking the cleavage site."

flanking = windowsize/2

start = middle - flanking
if start<0: 
	print "Your window size is out of bounds so I'm going to start it at the beginning of the file."
if windowsize%2==0:
	end = middle+flanking+2 #adding the +1 because of how Python slices arrays


#open genome file and extract window

genomefile = "../B_subtilis_NC_000964.txt"  #edit path as needed
#genomefile = "hello.txt" #use for debugging purposes
genome = open(genomefile, 'r')
genomelist = []
for line in genome: 
	if '>' in line:
		pass
	for char in line:
		genomelist.append(char)
if end>len(genomelist):
	print "Your window size is out of bounds so I'm going to end it at the end of the file."
	end = len(genomelist)
genome.close()
window = genomelist[start:end]
endwindow = "".join(window)
finalwindow = endwindow.rstrip('\n')

#reverse complement if the transcript is on the minus strand
#either way, print out the window here
strandedness = raw_input("Is this on the + strand?\n")
if 'n' in strandedness or 'N' in strandedness: 
	reverse = reversal(finalwindow)
	comp = seqcomplement(reverse, 'DNA')
	print "Here's the window for the minus strand: \n"
	print comp
else: 
	print "Here's the window for the plus strand: \n"
	print finalwindow
