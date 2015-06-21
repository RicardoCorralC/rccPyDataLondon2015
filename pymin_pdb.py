from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Residue

def get_aa_residues(pdb, chain):
    '''
    pdb: Protein Data Bank file.
    chain: Chain of the PDB file.

    Get the amino acids from a protein.

    returns: List of Biopython PDB Residue objects representing the amino acids
    of the specified protein.
    '''
    parser = PDBParser()
    structure = parser.get_structure('prot', pdb)
    model=structure[0]
    chain=model[chain]

    # Get a list of all residues in the specified protein model.
    residue_list = list(chain.get_residues())
    to_remove_list = []

    for res in residue_list:
        # Store non-amino acid residues in PDB in another list.
        if res.get_id()[0] != ' ':
            to_remove_list.append(res)
    
    # Remove non-amino acid residues from original list.
    for res in to_remove_list:
        residue_list.remove(res)

    return residue_list

def general_center_of_mass(pdb, chain):
    '''
    pdb: Protein Data Bank file.
    chain: Chain of the PDB file.
    
    Get the General Center of Mass from a specified protein.

    returns: Tuple of x,y,z coordinates of the GCM.
    '''

    # Molecular weight of atoms in proteins.
    mw = {'C':12.017, 'N':14.0067, 'O':15.9994, 'S':32.065, 'Se':78.96}
    
    # Get a list of all heavy atoms in the protein.
    residues = get_aa_residues(pdb, chain)
    heavy_atoms = []

    for res in residues:
        for atom in res:
            # PDB files usually do not include the hydrogen atoms of a protein, so 
            # there's no need to filter atoms by type.
            heavy_atoms.append(atom)

    # Coordinates for GCM numerator and denominator.
    xn, yn, zn = [0.0] * 3
    d = 0 

    step = 0

    for atom in heavy_atoms:
        
        if atom.get_name()[:2] == 'SE':
            element = 'Se'
        else:
            element = atom.get_name()[0]
            
        mass = mw[element]
        coords = atom.get_coord()

        xn += coords[0] * mass
        yn += coords[1] * mass   
        zn += coords[2] * mass
        
        d += mass

    return (xn/d, yn/d, zn/d)

def distance_to_point(res, point):
    '''
    res: Biopython Residue object.
    point: x,y,z coordinates

    Calculates the distance from a residue to a point (e.g. the gcm of the protein)
    The coordinates of the residue are taken to be the coordinates of the alpha 
    carbon of the amino acid.

    returns: Distance to point (float). If the coordinates of the alpha carbon of the
    residue cannot be found, returns None.
    '''
    import math
    coords = res['CA'].get_coord()

    if len(coords) != 3:
        print 'Could not determine coordinates of alpha carbon for residue:'
        print str(res)
        return None
    
    distanceX = coords[0] - point[0]
    distanceY = coords[1] - point[1]
    distanceZ = coords[2] - point[2]

    return math.sqrt((distanceX*distanceX) + (distanceY*distanceY)
            + (distanceZ*distanceZ))

def list_dist_to_point(reslist, point):
    '''
    reslist: List of Biopython PDB Residues.

    Get the distances from every residue in the list (alpha carbon) to a point.

    returns: Dictionary with residues as keys and distance to GCM as values. If coords
    for alpha carbon of a residue cannot be found, None is returned as the value
    in the dictionary for that residue.

    '''

    dist = {}

    for res in reslist:
        dist[res] = distance_to_point(res, point)

    return dist

def get_aa_name(res):
    '''
    res: Biopython PDB Residue object representing an amino acid.

    returns: Name of residue in three letter code + residue number format (e.g. LYS23)
    '''
    return res.get_resname() + str(res.get_id()[1])
