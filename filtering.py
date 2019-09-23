import sys

# Filter by length of V sequence and J gene assignment
nf = open("filteredIGKL.csv", 'w')
with open(sys.argv[1], 'r') as file:
    next(file)
    for line in file:
        x = line.split(",")
        try:
            jgene = str(x[23])
            length = int(x[13])
        except ValueError as err:
            pass
        if jgene == (str("-")) and length < 50:
            continue
        else:    
            nf.write(line)

nf.close()


        
