import os, sys

def getFirstChain(pdbfile):
	"""Returns first one letter code chain identifier encountered in a PDB"""
	pdb = open(pdbfile,'r')
	for l in pdb:
		if l[:4] != 'ATOM': continue
		chain = l[21]
		break
	pdb.close()
	return chain
	
def main():
	pdb = sys.argv[1]
	print getFirstChain(pdb)

if __name__ == '__main__':
	main()
