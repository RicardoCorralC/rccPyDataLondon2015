
setSignatures = ((list( [0 ])),
				(list( [1 ])),
				(list( [1,1 ]),list( [2 ])),
				(list( [1,1,1 ]),list( [1,2 ]),list( [3 ])),
				(list( [1,1,1,1 ]),list( [1,1,2 ]),list( [2,2 ]),list( [1,3 ]),list( [4 ])),
				(list( [1,1,1,1,1 ]),list( [1,1,1,2 ]),list( [1,1,3 ]),list( [1,4 ]),list( [1,2,2 ]),list( [2,3 ]),list( [5 ])),
				(list( [1,1,1,1,1,1 ]),list( [1,1,1,1,2 ]),list( [1,1,1,3 ]),list( [1,1,2,2 ]),list( [1,1,4 ]),list( [1,2,3 ]),list( [1,5 ]),list( [2,2,2 ]),list( [2,4 ]),list( [3,3 ]),list( [6 ])),
	)

def iter_classes():
	for i in xrange(3,7):
		for c in setSignatures[i]:
			yield c

PyMol_header_string="""
hide all
set ray_shadows,0
set antialias = 1
set ray_trace_fog, 5
set fog, on
set fog, 2
set depth_cue, 1
set cartoon_fancy_helices, 1
bg_color white
color white, all
show cartoon, chain  *X*
"""

"""String lambda binding PyMol header"""
PyMol_header = lambda s : PyMol_header_string.replace('*X*',s)

AA_three_to_one = {
				'ALA':'A',
				'ARG':'R',
				'ASN':'N',
				'ASP':'D',
				'CYS':'C',
				'GLU':'E',
				'GLN':'Q',
				'GLY':'G',
				'HIS':'H',
				'ILE':'I',
				'LEU':'L',
				'LYS':'K',
				'MET':'M',
				'PHE':'F',
				'PRO':'P',
				'SER':'S',
				'THR':'T',
				'TRP':'W',
				'TYR':'Y',
				'VAL':'V',
			}

AA_one_to_three = {
				'A':'ALA',
				'R':'ARG',
				'N':'ASN',
				'D':'ASP',
				'C':'CYS',
				'E':'GLU',
				'Q':'GLN',
				'G':'GLY',
				'H':'HIS',
				'I':'ILE',
				'L':'LEU',
				'K':'LYS',
				'M':'MET',
				'F':'PHE',
				'P':'PRO',
				'S':'SER',
				'T':'THR',
				'W':'TRP',
				'Y':'TYR',
				'V':'VAL',
			}

