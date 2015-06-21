import pymin_pdb as pm_pdb
import networkx as nx
from networkx import Graph
from Bio.PDB.PDBParser import PDBParser

def dk(graph):
    '''
    graph: A networkx Graph object to calculate the dk from.

    Calculates the dynamic connectivity (dk) of a network. This is a centrality
    measure that is related to betweenness centrality. The shortest path
    between every pair of nodes in the network is calculated with Dijkstra's 
    algorithm. The dk of a node is the number of times it lies within one of 
    these paths.

    This centrality measure was used in the program JAMMING, which predicts
    functional residues of a protein through network analysis. 
    
    Reference:
    Efficient identification of critical residues based only on protein structure
    by network analysis
    Cusack M, Thibert B, Bredesen DE and del Rio G. 2007. PLoS ONE 2(5):e421
    
    returns: Dictionary with nodes as key and dk for each node as value.
    '''
    assert type(graph) is nx.Graph
    
    # Get a list of the nodes of the graph.
    nodes = graph.nodes()

    # Only for debugging. Sort amino acids by key.
    nodes.sort(key=lambda x: x.get_id()[1])

    # Create list of same length as node list with zero as all entries.
    dk_list = [0] * len(nodes)

    # Calculate the dijkstra shortest path for all node pairs.
    sp = nx.all_pairs_dijkstra_path(graph)

    # For all nodes in the graph
    for i in range(0, len(nodes)):
        node1 = nodes[i]

        # For all nodes later in list than first node
        for j in range(i+1, len(nodes)):
            node2 = nodes[j]
            # Get shortest path between the first and second node
            path = sp[node1][node2]
            # Add one to the dk_list at an index for every time the node of that index
            # is inside the shortest path between a node pair.
            for n in path[1:-1]:
                dk_list[nodes.index(n)] += 1

    # Store dk values in a dictionary.
    dictionary = dict(zip(nodes, dk_list))
    return dictionary

def jamming_dk(graph):
    '''
    graph: A networkx Graph object to calculate the dk from.

    Calculates the dynamic connectivity (dk) of a network. This is a centrality
    measure that is related to betweenness centrality. The shortest path
    between every pair of nodes in the network is calculated with Dijkstra's 
    algorithm. The dk of a node is the number of times it lies within one of 
    these paths. This function attempts to replicate the results of JAMMING
    (including quirks).

    This centrality measure was used in the program JAMMING, which predicts
    functional residues of a protein through network analysis. 
    
    Reference:
    Efficient identification of critical residues based only on protein structure
    by network analysis
    Cusack M, Thibert B, Bredesen DE and del Rio G. 2007. PLoS ONE 2(5):e421
    
    Note: Currently, this implementation assumes the graph has Bio.PDB Residue
    objects as nodes, so this doesn't work for general graphs. 
    
    Besides, this implementation is terribly inefficient. It might be better to
    somehow force NetworkX to have the residues with sequential ids and calculate 
    the Dijkstra shortest paths.

    returns: Dictionary with nodes as key and dk for each node as value.
    '''
    # TODO
    # Still doesn't give same results as JAMMING.
    def get_jamming_shortest_paths(graph, nodelist):
        '''
        Chooses the shortest path JAMMING would choose.
        '''
        # Sort list by sequence number. This implementation only works if
        # the nodes are a list of Bio.PDB Residue objects representing amino acids.
        nodelist.sort(key=lambda x: x.get_id()[1])

        jamming_paths = []

        for i in range(len(nodelist)):
            for j in range(i+1, len(nodelist)):
                shortest_paths = list(nx.all_shortest_paths(graph, nodelist[i], nodelist[j]))
                if len(shortest_paths) == 1:
                    jamming_paths.append(shortest_paths[0])
                elif len(shortest_paths) > 1:
                    jamming_path = min(shortest_paths, key=lambda x: x[-2].get_id()[1])
                    jamming_paths.append(jamming_path)

        return jamming_paths

    assert type(graph) is nx.Graph
    
    # Get a list of the nodes of the graph.
    nodes = graph.nodes()

    # Only for debugging. Sort amino acids by key.
    nodes.sort(key=lambda x: x.get_id()[1])

    # Create list of same length as node list with zero as all entries.
    dk_list = [0] * len(nodes)

    # Get the shortest paths for the graph as JAMMING would. 
    sp = get_jamming_shortest_paths(graph, nodes)

    # Add one to dk_list every time a node is in a shortest path at that node's position.
    for path in sp:
        for n in path:
            dk_list[nodes.index(n)] += 1

    # Store dk values in a dictionary.
    dictionary = dict(zip(nodes, dk_list))
    return dictionary


def get_distance(coords1, coords2):
    '''
    coords1: A set of (x,y,z) coordinates.
    coords2: A set of (x,y,z) coordinates.

    From the (x,y,z) coordinates of two objects (e.g. pair of atoms), 
    calculates the distance between them.

    returns: Distance between pair of coordinates.
    '''
    import math

    distancex = coords1[0] - coords2[0]
    distancey = coords1[1] - coords2[1]
    distancez = coords1[2] - coords2[2]
    return math.sqrt((distancex*distancex) + (distancey*distancey) 
            + (distancez*distancez))

def build_unweighted_psn(residues, distance_cutoff, verbose=False):
    '''
    residues: Biopython.PDB Residue objects representing the amino acids of a 
    protein.

    distance_cutoff: Distance cutoff in Angstroms for edge (interaction) inclusion 
    between two nodes (amino acids).
    
    Builds a protein structure network. That is, a network which nodes are
    the amino acids of a protein and an edge between two nodes represents an 
    interaction between the corresponding amino acids. Amino acids that have at 
    least one atom pair between them within the specified cutoff distance are assumed
    to interact and an edge between them is included in the network.

    returns: NetworkX Graph object representing the protein structure network.
    '''

    def within_cutoff(res1, res2, distance_cutoff, verbose=False):
        '''
        res1: BioPython PDB module Residue object.
        res2: BioPython PDB module Residue object.

        Calculates whether two residues have at least an atom pair within the 
        specified cutoff distance.

        returns: True if residues are within cutoff distance, False if they are not.
        '''
        for atom1 in res1:
            atom1coords = atom1.get_coord()
            for atom2 in res2:
                atom2coords = atom2.get_coord()
                distance = get_distance(atom1coords, atom2coords)
                if (distance <= distance_cutoff):
                    return True
        return False
        
    ###########################################################################

    distance_cutoff = float(distance_cutoff)

    # Create the network and add amino acid residues from the pdb file as nodes.
    if verbose == True:
        print 'Initializing network...'            
    network = nx.Graph()

    if verbose == True:
        print 'Adding amino acids as nodes to network...'
    network.add_nodes_from(residues)    

    # Add an edge between every pair of amino acids that has at least one atom pair
    # within the distance cutoff.
    if verbose == True:
        print 'Adding edges to network based on distance criterion...'
    for i in range(0, len(residues)):
        residue1 = residues[i]
        for j in range(i+1, len(residues)):
            residue2 = residues[j]
            if within_cutoff(residue1, residue2, distance_cutoff):
                network.add_edge(residue1, residue2)
    
    if verbose == True:
        print 'Done building network.'

    return network

def build_weighted_psn_CA(residues, distance_dependence, cutoff, verbose):
    '''
    Builds a weighted PSN with distance dependent edge weights and distance
    between amino acids = distance between alpha carbons. This is an internal 
    function not intended for direct use. Please use build_weighted_psn() instead.
    '''

    # Create the network and add amino acid residues as nodes.
    if verbose == True:
        print 'Initializing network...'            
    network = nx.Graph()

    if verbose == True:
        print 'Adding amino acids as nodes to network...'
    network.add_nodes_from(residues)

    if verbose == True:
        print 'Adding edges to network with a distance dependence of',
        print str(distance_dependence) + '.'
        print 'Distance for edge weight calculation is distance between alpha carbons.'

    # For every amino acid in the list:
    for i in range(0, len(residues)):
        residue1 = residues[i]
        # Get coordinates of the alpha carbon of the first amino acid.
        coords1 = residue1['CA'].get_coord()

        # For every amino acid in the list downstream in sequence to the first one:
        for j in range(i+1, len(residues)):
            residue2 = residues[j]
            # Get coordinates of the alpha carbon of the second amino acid:
            coords2 = residue2['CA'].get_coord()
            # Get distance between alpha carbons.
            dist = get_distance(coords1, coords2)
            # If a cutoff was specified
            if cutoff != None:
                # Check if distance is below cutoff, add weight to network only if
                # it is.
                if dist < cutoff:
                    # Set edge weight to the distance between the amino acids and add
                    # edge to network.
                    #
                    # Note: Centrality calculations assume lower edge weight values
                    # mean the residues are more closely connected, so the actual
                    # edge weight depends inversely on the assumed interaction strength.
                    #
                    edge_weight = dist**distance_dependence
                    network.add_edge(residue1, residue2, weight=edge_weight)
            else:
                # Same as above, but don't check for distance cutoff.
                edge_weight = dist**distance_dependence
                network.add_edge(residue1, residue2, weight=edge_weight)
            
    if verbose == True:
        print 'Done building network.'

    return network

def build_weighted_psn_closest(residues, distance_dependence, cutoff, verbose):
    '''
    Builds a weighted PSN with distance dependent edge weights and distance
    between amino acids = distance between closest atom pair. Internal function not
    intended for direct use. Please use build_weighted_psn() instead.
    '''
    def min_dist(res1, res2):
        '''
        res1: BioPython.PDB Residue object.
        res2: BioPython.PDB Residue object.

        Calculates whether two residues have at least an atom pair within the 
        specified cutoff distance.

        returns: True if residues are within cutoff distance, False if they are not.
        '''
        min_dist = 1000.0
        for atom1 in res1:
            atom1coords = atom1.get_coord()
            for atom2 in res2:
                atom2coords = atom2.get_coord()
                distance = get_distance(atom1coords, atom2coords)
                if (distance <= min_dist):
                    min_dist = distance
        return min_dist

    # Create the network and add amino acid residues as nodes.
    if verbose == True:
        print 'Initializing network...'            
    network = nx.Graph()

    if verbose == True:
        print 'Adding amino acids as nodes to network...'
    network.add_nodes_from(residues)

    if verbose == True:
        print 'Adding edges to network with a distance dependence of',
        print str(distance_dependence) + '.'
        print 'Distance for edge weight calculation is distance between closest atoms.'

    # For every amino acid in the list:
    for i in range(0, len(residues)):
        residue1 = residues[i]

        # For every amino acid in the list downstream in sequence to the first one:
        for j in range(i+1, len(residues)):
            residue2 = residues[j]
            # Get distance between closest atom pair of the two amino acids.
            dist = min_dist(residue1, residue2)
            # If a cutoff was specified
            if cutoff != None:
                # Check if distance is below cutoff, add weight to network only if
                # it is.
                if dist < cutoff:
                    # Set edge weight to the distance between the amino acids and add
                    # edge to network.
                    #
                    # Note: Centrality calculations assume lower edge weight values
                    # mean the residues are more closely connected, so the actual
                    # edge weight depends inversely on the interaction strength.
                    #
                    edge_weight = dist**distance_dependence
                    network.add_edge(residue1, residue2, weight=edge_weight)
            else:
                # Same as above, but don't check for distance cutoff.
                edge_weight = dist**distance_dependence
                network.add_edge(residue1, residue2, weight=edge_weight)
            
    if verbose == True:
        print 'Done building network.'

    return network

def build_weighted_psn(residues, distance_dependence=2, Ca=True, cutoff=None, verbose=False):
    '''
    residues: List of Biopython.PDB Residue objects representing the amino acids of a 
    protein.

    distance_dependence: Distance dependence of the strength of an interaction between
    amino acids. (e.g. if distance_dependence=2, then strength of interaction is assumed to
    decay as the square of the distance).

    Ca: If True, consider the distance between a pair of amino acids the distance
    between their alpha carbons. If False, consider the distance between a pair of 
    amino acids the distance between their closest atom pair.

    cutoff: Optional parameter. Floating point value specifying the distance cutoff
    between two amino acids for edge inclusion in the network (regardless of edge 
    weight). Default is None, which means the network constructed is a complete 
    (weighted) network. e.g. if cutoff=12.0, no amino acids with a distance higher
    than 12 angstroms are connected by an edge.

    verbose: Optional parameter. Boolean value, default=True. Output intermediate 
    steps to terminal as text. Useful for debugging.

    Builds a weighted protein structure network from the coordinates in the specified 
    pdb and chain. The weight of the nodes is distance dependent. 

    The network constructed by this method is a complete network by default (every amino
    acid pair is connected by an edge). A non-complete network can be constructed with
    the optional parameter cutoff. The weight of the edges is assumed to be dependent
    on a power of the distance. The optional parameter distance_dependence (default=2) 
    is the distance dependence of the strength of the interaction (e.g. dist_dep=3 means 
    the strength of the interaction decays as the third power of the distance).

    returns: A NetworkX Graph object representing the protein structure network.
    '''
    assert type(residues) == list
    assert type(distance_dependence) == int or type(distance_dependence) == float 
    assert type(Ca) == bool
    assert type(cutoff) == float or cutoff == None

    if Ca == True:
        psn = build_weighted_psn_CA(residues, distance_dependence, cutoff, verbose)
    else:
        psn = build_weighted_psn_closest(residues, distance_dependence, cutoff, verbose)
    return psn
