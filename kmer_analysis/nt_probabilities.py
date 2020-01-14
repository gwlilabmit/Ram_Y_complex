'''
Takes the 168 genome from subtilis.fa and calculates p(a), p(g), etc.
'''

infile = open("subtilis.fa", "r")

total = 0

a = 0
g = 0
c = 0
t = 0

for line in infile: 
	total += len(line[:-1])
	a += line.count("A")
	t += line.count("T")
	c += line.count("C")
	g += line.count("G")

infile.close()

proba = float(a)/float(total)
probt = float(t)/float(total)
probg = float(g)/float(total)
probc = float(c)/float(total)


print "p(a)= "+ str(proba)
print "p(t)= "+ str(probt)
print "p(g)= "+ str(probg)
print "p(c)= "+ str(probc)

