#!/usr/bin/python

from heapq import nlargest
import sys
import os
 
## Get Top 3 Genes ##
abundance = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        x = line.split("\t")
        try:
            counts = x[4]
            gene = x[7]
            noallele = gene.split("*", 1)[0]
            k = noallele
            v = counts
            abundance[k]=(v)
        except IndexError as err:
            pass
alltop3 = nlargest(3, abundance, key=abundance.get)
print("Top 3 genes: ", alltop3)

## Get Top Pseudogenes ##
allpseudo = []
with open("../pseudoV.list", 'r') as pseudoV:
    for line in pseudoV:
        x = line.strip("\n")
        allpseudo.append(x)
with open("../pseud.K", 'r') as pseudoK:
    for line in pseudoK:
        x = line.strip("\n")
        allpseudo.append(x)
with open("../pseud.L", 'r') as pseudoL:
    for line in pseudoL:
        x = line.strip("\n")
        allpseudo.append(x)
matches = {}
for k in allpseudo:
    if k in abundance:
        key = k
        v = abundance[k]
        if v == 0:
            continue
        matches[key]=(v)
top3pseudo = nlargest(3, matches, key=matches.get)
print("Top 3 pseudogenes:", top3pseudo)

## Get Top 3 Heavy Chain- functional and pseudo ##
Hdict = {}
Hdictpseudo = {}
H = 'IGHV'
heavy = [i for i in abundance if H in i]
for k in heavy:
    if k in abundance:
        key = k
        v = abundance[k]
        Hdict[key] = v
top3 = nlargest(3, Hdict, key=Hdict.get)
print("Top 3 Heavy:", top3)

heavy = [i for i in matches if H in i]
for k in heavy:
    if k in matches:
        key = k
        v = matches[k]
        if v ==0:
            continue
        Hdictpseudo[key] = v
top3 = nlargest(3, Hdictpseudo, key=Hdictpseudo.get)
print("Top 3 Heavy Pseudogenes:", top3)

## Get Top 3 Kappa- functional and pseudo ##
Kdict = {}
Kdictpseudo = {}
K = 'IGKV'
kappa = [i for i in abundance if K in i]
for k in kappa:
    if k in abundance:
        key = k
        v = abundance[k]
        Kdict[key] = v
top3 = nlargest(3, Kdict, key=Kdict.get)
print("Top 3 Kappa:", top3)

kappa = [i for i in matches if K in i]
for k in kappa:
    if k in matches:
        key = k
        v = matches[k]
        if v ==0:
            continue
        Kdictpseudo[key] = v
top3 = nlargest(3, Kdictpseudo, key=Kdictpseudo.get)
print("Top 3 Kappa Pseudogenes:", top3)

## Get Top 3 Lambda- functional and pseudo ##
Ldict = {}
Ldictpseudo = {}
L = 'IGLV'
Lambda = [i for i in abundance if L in i]
for k in Lambda:
    if k in abundance:
        key = k
        v = abundance[k]
        Ldict[key] = v
top3 = nlargest(3, Ldict, key=Ldict.get)
print("Top 3 Lambda:", top3)

Lambda = [i for i in matches if L in i]
for k in Lambda:
    if k in matches:
        key = k
        v = matches[k]
        Ldictpseudo[key] = v
        if v ==0:
            continue
top3 = nlargest(3, Ldictpseudo, key=Ldictpseudo.get)
print("Top 3 Lambda Pseudogenes:", top3)

## Check if top hit is pseudogene ##
print(set(alltop3) & set(top3pseudo))


            
