import uuid, string, random, datetime, difflib
import math, re
from features import RNA, DNA
from bson.objectid import ObjectId

def get_points(x1, y1, x2, y2, distance):
    opposite_side = get_distance(x1, y1, x1, y2)
    adjacent_side = get_distance(x1, y2, x2, y2)
    if not adjacent_side:
        return []
    angle = get_angle(opposite_side, adjacent_side)
    new_x1 = None
    new_y1 = None
    new_x2 = None
    new_y2 = None
    
    if x1 >= x2:
        new_x2 = x2 + get_adjacent_side(angle, distance)
        new_x1 = x1 - get_adjacent_side(angle, distance)
    else:
        new_x2 = x2 - get_adjacent_side(angle, distance)
        new_x1 = x1 + get_adjacent_side(angle, distance)
    
    if y1 >= y2:
        new_y2 = y2 + get_opposite_side(angle, distance)
        new_y1 = y1 - get_opposite_side(angle, distance)
    else:
        new_y2 = y2 - get_opposite_side(angle, distance)
        new_y1 = y1 + get_opposite_side(angle, distance)
        
    return [[new_x1, new_y1], [new_x2, new_y2]]
    
def get_angle(opposite_side, adjacent_side):
    return math.atan(opposite_side/adjacent_side)
    
def get_distance(x1, y1, x2, y2):
    horizontal = x1 - x2
    vertical = y1 - y2
    return math.sqrt(horizontal*horizontal + vertical*vertical)

def get_adjacent_side(angle, hypothenuse):
    return math.cos(angle) * hypothenuse

def get_opposite_side(angle, hypothenuse):
    return math.sin(angle) * hypothenuse

"""
This methods finds clusters from a list of genomic annotations.
Args:
- genomic_annotations: a list of genomic annotations. Each genomic annotation should be described as a dict like {'genomicStart':int, 'genomicEnd':int}
- threshold: this methods will only return the clusters whose size is >= threshold (default: 1)
- fill_cluster_with_genomic_annotations: if True, all the genomic annotations making the cluster are stored in a list linked to the key 'genomic_annotations' (default: False).
Returns:
- a list of cluster. Each cluster is a dict like {'genomicStart': x:int, genomicEnd': x:int, 'annotations_count': x:int, 'genomic_annotations': x:list}
"""
def cluster_genomic_annotations(genomic_annotations, threshold = 1, fill_cluster_with_genomic_annotations = False):
    clusters = []
    sorted_genomic_annotations = sorted(genomic_annotations, key=lambda genomic_annotation: genomic_annotation['genomicStart'])
    first_genomic_annotation = sorted_genomic_annotations[0]
    cluster  = { 
        'genomicStart': first_genomic_annotation['genomicStart'],
        'genomicEnd': first_genomic_annotation['genomicEnd'],
        'annotations_count': 1,
        'genomic_annotations': [first_genomic_annotation]
    }

    for i in xrange(1, len(sorted_genomic_annotations)):
        genomic_annotation = sorted_genomic_annotations[i]
        if genomic_annotation['genomicStart'] <= cluster['genomicEnd']:
            cluster['annotations_count'] += 1
            if fill_cluster_with_genomic_annotations:
                cluster['genomic_annotations'].append(genomic_annotation)    
            if genomic_annotation['genomicEnd'] > cluster['genomicEnd']:
                cluster['genomicEnd'] = genomic_annotation['genomicEnd']                       
        else:
            if cluster['annotations_count'] >= threshold:
                clusters.append(cluster)
            cluster = { 
                'genomicStart': genomic_annotation['genomicStart'],
                'genomicEnd': genomic_annotation['genomicEnd'],
                'annotations_count':1,
                'genomic_annotations': [genomic_annotation]
            }
    # last cluster
    if cluster['annotations_count'] >= threshold:
        clusters.append(cluster)    

    return clusters

def dataframe_to_json(dataframe):
    import ujson
    d = [ 
        dict([
            (colname, row[i]) 
            for i, colname in enumerate(dataframe.columns)
        ])
        for row in dataframe.values
    ]
    return ujson.dumps(d)

def make_random_molecule(size, name="unnamed", type="RNA"):
    """
    Create a molecule with a random sequence of size given as argument
    """
    residues = []
    if type =='RNA':
        residues += ['A','U','G','C']
    elif type == 'DNA':
        residues += ['A','T','G','C']
    seq = list(random.choice(residues) for x in range(size))
    if type =='RNA':
        return RNA(sequence = ''.join(seq), name = name)
    elif type == 'DNA':
        return DNA(sequence = ''.join(seq), name = name)

def get_file_as_source(absolute_path):
    import socket
    return ':'.join(['file', socket.gethostname(), absolute_path])

def remove_spaces(s, substitution_char='_'):
    """
    Replace all the spaces characters with the char of substitution
    """
    import re
    return re.sub(r"\s+", substitution_char, s)

def get_time(uuid1):
    """
    Return the date time of an uuid1.
    uuid1: the uuid1 as a string
    """
    t = uuid.UUID(uuid1)
    return datetime.datetime.fromtimestamp((t.time - 0x01b21dd213814000L) * 100 / 1e9)

def chunks(l, n):
    """
     Return l as a list of tuples of size n
     """
    return [tuple(l[i:i + n]) for i in range(0, len(l), n)]

def chuncks_with_overlap(l,n,overlap):
    return [l[i:i+n] for i in range(0, len(l), n - overlap)]

def generate_random_name(n):
    """
    Generate a random name of n letters
    """
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(n))

def is_canonical(residue_1, residue_2, orientation, edge_1, edge_2):
    return (residue_1.upper() == 'A' and residue_2.upper() == 'U' or residue_1.upper() == 'U' and residue_2.upper() == 'A' or\
            residue_1.upper() == 'G' and residue_2.upper() == 'C' or residue_1.upper() == 'C' and residue_2.upper() == 'G' or\
            residue_1.upper() == 'G' and residue_2.upper() == 'U' or residue_1.upper() == 'U' and residue_2.upper() == 'G') and\
            orientation.lower() == 'c' and edge_1.upper() == '(' and edge_2.upper() == ')'

def get_atoms_distance(a1,a2):
    """
    Return the distance between two atoms.
    """
    x_dist = (a1[0] - a2[0])**2
    y_dist = (a1[1] - a2[1])**2
    z_dist = (a1[2] - a2[2])**2
    return math.sqrt(x_dist + y_dist + z_dist)

def get_levenshtein_distance(a,b):
    """
    Return the Levenshtein distance between two strings
    """
    return difflib.SequenceMatcher(None, a,b).ratio()

def get_distances(molecules, ref):
    """
    Compute the Levenshtein distance between ref and all the other molecules. Return a list of tuples (mol2, distance) as values sorted according to the distance
    - molecules: an array of Molecule objects
    - ref: the reference molecule
    """
    distances=[]
    for i in range(0,len(molecules)):
        mol1 = molecules[i]
        if mol1 != ref:
            continue
        for mol2 in molecules[0:len(molecules)]:
            if mol2 == mol1:
                continue
            distances.append((mol2,get_levenshtein_distance(mol1.sequence,mol2.sequence)))
        if mol1 == ref:
            break
    return sorted(distances, key= itemgetter(1), reverse=True)

def find_lcs(a,b):
    """
    Return the longest common substring between two strings
    """
    return difflib.SequenceMatcher(None, a,b).find_longest_match(0,len(a),0,len(b))

def get_molecule_by_name(molecules, molecule_name):
    hits = [i for i in molecules if i.name == molecule_name]
    if len(hits) == 0:
        return None
    return hits.pop()

def renumber_pdb_atoms(pdb_file):
    """
    Renumber all atoms in a pdb file, starting from 1.

    Args:
    - pdb_file: the path of the pdb file to renumber

    Returns:
    The renumbered pdb file content as a String
    """

    new_lines = []
    counter = 1
    with open(pdb_file) as h:
        for line in h:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                new_lines.append("%s%5s%s" % (line[0:6],counter,line[11:]))
                counter += 1
            else:
                new_lines.append(line)
    return ''.join(new_lines)



