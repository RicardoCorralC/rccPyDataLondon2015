import networkx as nx
import RCCdata as rcd
import os
from scipy.spatial.distance import euclidean,cosine
from numpy import zeros


TM_align_bin = '/home/rcc/Desktop/ClusterIngredients/RCCpackage/TMalignc/TMalign' #WARNING hard path

def resCompare(xx,yy):
	"""Given a three letter code residue name, it says what is first in sequence. 'ASN37 > GLN15'"""
	x = int(xx[3:])
	y = int(yy[3:])
	return x-y


def getFirstChain(pdbfile):
	"""Returns first one letter code chain identifier encountered in a PDB"""
	pdb = open(pdbfile,'r')
	for l in pdb:
		if l[:4] != 'ATOM': continue
		chain = l[21]
		break
	pdb.close()
	return chain


def findIsomophicRegions(P1,P2,sizeRegion=11): 
	"""Generator of isomorphic regions """
	graphs1 = dict()
	graphs2 = dict()
	for res1 in P1.G.nodes():
		graphs1[res1] = P1.G.subgraph(P1.regions[res1])

	for res2 in P2.G.nodes():
		graphs2[res2] = P2.G.subgraph(P2.regions[res2])
		
	for res1 in P1.G.nodes():
		for res2 in P2.G.nodes():
			if len(graphs1[res1].nodes()) < sizeRegion or len(graphs2[res2].nodes()) < sizeRegion : continue
			if nx.is_isomorphic(graphs1[res1],graphs2[res2]):
				yield dict(residues=(res1,res2),graphs=(graphs1[res1],graphs2[res2]))

def AlignFragmentPairs(Proteina,Proteina2,sizeRegion=10):
	"""Dynamic programming aligner of residues with same environment"""
	matches = findIsomophicRegions(Proteina,Proteina2,sizeRegion=sizeRegion)
	amatches, bmatches = list(), list()
	pairmatches = list()
	for mm in matches:
		mr = mm['residues']
		a,b = mr[0], mr[1]
		pairmatches.append(tuple([a,b]))
		amatches.append(a)
		bmatches.append(b)

	amatches.sort(cmp=resCompare)
	bmatches.sort(cmp=resCompare)
	anames, bnames = dict(list(enumerate(amatches,start=1))), dict(list(enumerate(bmatches,start=1)))
	apos, bpos = new_dict = dict (zip(anames.values(),anames.keys())) , dict (zip(bnames.values(),bnames.keys()))

	#Aligned Fragment Pairs Cost Matrix; 1 if ri and rj have isomorphic environments
	AFPM = zeros([len(amatches)+1,len(bmatches)+1])
	AFPwalk = zeros([len(amatches)+1,len(bmatches)+1])
	for p in pairmatches:
		AFPM[apos[p[0]],bpos[p[1]]]=1

	for i in xrange(1,len(amatches)+1):
		for j in xrange(1,len(bmatches)+1):
			d = AFPM[i-1,j-1]+AFPM[i,j]
			gapi = AFPM[i-1,j]
			gapj = AFPM[i,j-1]
			AFPM[i,j] = max(d,gapi,gapj)
			if d > max(gapi, gapj):
				AFPwalk[i,j] = 3
			elif gapi > gapj:
				AFPwalk[i,j] = 2
			else:
				AFPwalk[i,j] = 1

	maxalign = AFPM[i,j]
	print 'maximo align: ',maxalign

	while(i>0 and j>0):
		if AFPwalk[i,j] == 3:
			yield (anames[i],bnames[j])
			i-=1
			j-=1
		elif AFPwalk[i,j] == 2:
			i-=1
		else:
			j-=1


def setToSequece(residueSet,gaps=True,gap_char='_',withNum=False,threeLetter=False): #three letter code amino name and position to only a one letter code sequence
		positions = []
		aaPositionName = dict()
		for res in residueSet:
			position = int(res[3:])
			positions.append(position)
			if threeLetter:
				aaPositionName[position] = res[:3]
			else:
				aaPositionName[position] = str(rcd.AA_three_to_one[res[:3]])
		positions = sorted(positions)
		amino_sequence = str(aaPositionName[positions[0]])
		if withNum: amino_sequence += str(positions[0]) + ' '
		for i in xrange(1,len(positions)):
			if gaps and (positions[i] != (positions[i-1]+1)):
				amino_sequence += gap_char + ' '
			amino_sequence += aaPositionName[positions[i]] 
			if withNum: amino_sequence += str(positions[i]) + ' '
		return amino_sequence

def iter_directory_files(maindir):
	for dirname, dirnames, filenames in os.walk(maindir):
		for filename in filenames:
			yield os.path.join(dirname,filename)


def disassemble_NMR(PDBfile):
	"""Create a directory cointaining each conformation reported in PDBfile in a different file"""
	fin = open(PDBfile,'r')
	basename = PDBfile[:PDBfile.find('.')]
	if not os.path.exists(basename):
	    os.makedirs(basename)
	print basename
	confnum = 1
	tmpfn = basename + str(confnum) + '.pdb'
	first = True
	prevresnum = 0
	for linea in fin:
		fout = open(os.path.join(basename,tmpfn),'a')
		if linea.split()[0] == 'ATOM':
			resnum = int(linea[23:27])
			print prevresnum, resnum
			if not first:
				if prevresnum <= resnum:
					fout.write(linea)
				else:
					confnum += 1
					tmpfn = basename + str(confnum) + '.pdb'
					fout = open(os.path.join(basename,tmpfn),'a')
					fout.write(linea)
			else:
				fout.write(linea)
				first = False
			prevresnum = resnum



def pdb_residue_offset(pdbname,chain):
	fin = open(pdbname,'r')
	offset = 0
	for f in fin:
		if f[0:6] == "ATOM  "  and  f[21]==chain:
			offset = int(f[23:26])
			return offset

def TM_align(pdb1,pdb2):
	"""This call TMalign from Zhang lab University of Michigan, returns RMSD"""
	os.system(TM_align_bin + ' -A ' + pdb1 + ' -B ' + pdb2 + ' > TMalign.tmp')
	os.system('cat TMalign.tmp | grep RMSD > TMalignline.txt')
	fin = open('TMalignline.txt','r')
	RMSD = fin.readline().strip().split(',')[1]
	RMSD = RMSD[RMSD.find('=')+1:]
	return dict(RMSD=RMSD)
	
	
def TM_aligned_residues(pdb1,pdb2,offste1=0, offset2=0):
	"""This call TMalign from Zhang lab University of Michigan, returns tuples with aligned residues"""
	os.system('TMalign -A ' + pdb1 + ' -B ' + pdb2 + ' > TMalign.tmp')
	fin = open('TMalign.tmp','r')
	while True:
		f = fin.readline().strip()
		if f.startswith('(":" denotes'):
			sA = fin.readline().strip()
			matchs = fin.readline().replace(' ','@').strip()
			sB = fin.readline().strip()
			break
	print len(sA), len(sB), len(matchs)
	idxa = offste1
	idxb = offset2
	aligns1, aligns2 = [], []
	NOTaligns1, NOTaligns2 = [], []
	for i in xrange(len(matchs)):
		if matchs[i] == ':':
			ra = rcd.AA_one_to_three[sA[i]]+str(idxa)
			rb = rcd.AA_one_to_three[sB[i]]+str(idxb)
			aligns1.append( ra )
			aligns2.append( rb )
		if sA[i]!= '-':
			ra = rcd.AA_one_to_three[sA[i]]+str(idxa)
			NOTaligns1.append( ra )
			idxa += 1
		if sB[i]!= '-':
			rb = rcd.AA_one_to_three[sB[i]]+str(idxb)
			NOTaligns1.append( ra )
			idxb += 1
	return dict(alignedList1=aligns2,NOTalignedList1=NOTaligns2,alignedList2=aligns1,NOTalignedList2=NOTaligns1,seqA=sA,seqB=sB,matchs=matchs) #Sure, TM-aling shows the query sequence first

def distance(pdb1,pdb2,metric='euclidean'):
	"""Gives the distance betwen RCCs distribution vectors"""
	if metric=='cosine':
		return cosine(pdb1.RCCvector, pdb2.RCCvector)
	return euclidean(pdb1.RCCvector, pdb2.RCCvector)

def distance2(pdb1,pdb2,metric='euclidean'):
	"""Gives the distance betwen RCCs distribution vectors"""
	if metric=='cosine':
		return cosine(pdb1.RCCvector2, pdb2.RCCvector2)
	return euclidean(pdb1.RCCvector2, pdb2.RCCvector2)
