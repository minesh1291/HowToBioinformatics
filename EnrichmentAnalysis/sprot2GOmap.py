
from Bio import SwissProt
import sys

UniprotIDsFile="../0_Data/Selected.Protein.IDs.txt"
FH=open(UniprotIDsFile,'r')

lines=FH.readlines()

IDs= map(str.strip,lines)


from collections import defaultdict
#goComp = defaultdict(list)
#goFunc = defaultdict(list)
#goProc = defaultdict(list)
goAll = defaultdict(list)

c=0

for record in SwissProt.parse(open('../0_DB/uniprot_sprot.dat')):
	#print record.__dict__
	pid = record.accessions[0] #can be modified , only takng 1st identifier
	if pid in IDs :
		c=c+1
		print >>sys.stderr, c 
		for feature in record.cross_references:
			if feature[0] in ['GO']:
				goAll[pid].append(feature[1])
				
#				if feature[2][0] == "C":
#					goComp[pid].append(feature[1])
#				if feature[2][0] == "F":
#					goFunc[pid].append(feature[1])
#				if feature[2][0] == "P":
#					goProc[pid].append(feature[1])

for i in sorted(goAll.keys()):
	print i+"\t"+", ".join(goAll[i])
