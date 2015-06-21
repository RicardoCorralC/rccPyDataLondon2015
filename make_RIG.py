from pymin_pdb import *
from pymin_graph import *
import networkx as nx
import sys

def print_RIG(PDB, CHAIN, DISTANCE_CUTOFF=5.0):
	residues = get_aa_residues(PDB, CHAIN)
	network = build_unweighted_psn(residues, DISTANCE_CUTOFF)

	for a, b in network.edges():
		print get_aa_name(a), get_aa_name(b)

def main():
	PDB_ = sys.argv[1] #'1hiv.pdb'
	CHAIN_ = sys.argv[2] #'A'
	DISTANCE_CUTOFF_ = float(sys.argv[3]) #5.0
	print_RIG(PDB_, CHAIN_, DISTANCE_CUTOFF_)

if __name__ == '__main__':
	main()
