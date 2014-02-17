import uuid, string, random, datetime, difflib
import math, re
from features import RNA, DNA
from bson.objectid import ObjectId

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

def migrate(db_name, db_host = "localhost", db_port = 27017):
    """
    Using a MongoDB instance, this function migrates all the ncRNAs annotations from the table annotations to the table ncRNAs.
    If a ncRNA in the annotations table is not well "characterized", it is associated to the class '' and name ''.
    TODO: add the ncRNA to the alignment 
    """
    from pymongo import MongoClient
    client = MongoClient(db_host, db_port)
    db = client[db_name]

    before_ncRNAs = db['ncRNAs'].count()
    before_annotations = db['annotations'].count()

    print "Processing the migration of:"
    print "%i rRNAs from Genbank files"%db['annotations'].find({'class':'rRNA'}).count()
    print "%i tmRNAs from Genbank files"%db['annotations'].find({'class':'tmRNA'}).count()
    print "%i miscRNAs from Genbank files"%db['annotations'].find({'class':'misc_RNA'}).count()
    print "%i ncRNAs from Genbank files"%db['annotations'].find({'class':'ncRNA'}).count()

    print "At now, %i entries are stored in the ncRNAs table."%before_ncRNAs

    #the rRNAs have not been extracted from RFAM, so we add them directly to the ncRNAs table (and remove them from the annotations table)

    for rRNA in db['annotations'].find({'class':'rRNA'}):
        ncRNA_description = {
            '_id': str(ObjectId()),            
        }

        for key in rRNA.keys():
             ncRNA_description[key] = rRNA[key]

        ncRNA_description['class'] = 'Gene, rRNA'

        name = ""
        if rRNA.has_key('note'):
            name += rRNA['note']+" "

        if rRNA.has_key('product'):
            name += rRNA['product']

        if re.findall("5\.8", name):
            ncRNA_description['name'] = '5.8S'
            #What about the alignment??????????
            db['annotations'].remove(rRNA)
            db['ncRNAs'].insert(ncRNA_description)
        elif re.findall("5", name):
            ncRNA_description['name'] = '5S'
            #What about the alignment??????????
            db['annotations'].remove(rRNA)
            db['ncRNAs'].insert(ncRNA_description)
        elif re.findall("12", name) or re.findall("16", name) or re.findall("17", name) or re.findall("18", name) or re.findall("small", name, re.IGNORECASE) or re.findall("SSU", name, re.IGNORECASE):
            ncRNA_description['name'] = 'SSU_rRNA_eukarya'
            #What about the alignment??????????
            db['annotations'].remove(rRNA)
            db['ncRNAs'].insert(ncRNA_description)
        elif re.findall("23", name) or re.findall("25", name) or re.findall("26", name) or re.findall("28", name) or re.findall("large", name, re.IGNORECASE) or re.findall("LSU", name, re.IGNORECASE):
            ncRNA_description['name'] = 'LSU_rRNA_eukarya'
            #What about the alignment?????????? 
            db['annotations'].remove(rRNA)
            db['ncRNAs'].insert(ncRNA_description)
        else: 
            ncRNA_description['class'] = 'Misc RNA'
            ncRNA_description['name'] = 'Misc_RNA_from_NCBI'
            db['annotations'].remove(rRNA)
            db['ncRNAs'].insert(ncRNA_description)

    #the tmRNAs, miscRNAs and ncRNAs could overlap those extracted from RFAM. If so, they're simply removed from the annotations table. If not, they're new ncRNAs to add to the ncRNAs table
    
    for tmRNA in db['annotations'].find({'class':'tmRNA'}):
        ncRNA_description = {
            '_id': str(ObjectId()),            
        }

        for key in tmRNA.keys():
             ncRNA_description[key] = tmRNA[key]

        overlap = False

        for ncRNA in db['ncRNAs'].find({'genome':tmRNA['genome']}):
            if ncRNA['genome'] == tmRNA['genome'] and ncRNA['genomicStrand'] == tmRNA['genomicStrand'] and (ncRNA['genomicPositions'][0] >= tmRNA['genomicPositions'][0] and ncRNA['genomicPositions'][0] <= tmRNA['genomicPositions'][1] or ncRNA['genomicPositions'][1] >= tmRNA['genomicPositions'][0] and ncRNA['genomicPositions'][1] <= tmRNA['genomicPositions'][1] or tmRNA['genomicPositions'][0] >= ncRNA['genomicPositions'][0] and tmRNA['genomicPositions'][0] <= ncRNA['genomicPositions'][1] or tmRNA['genomicPositions'][1] >= ncRNA['genomicPositions'][0] and tmRNA['genomicPositions'][1] <= ncRNA['genomicPositions'][1]):
                overlap = True
                break
        
        if not overlap: #its a new ncRNA, we add it to the ncRNAs table and remove it from the annotations table
            ncRNA_description['class'] = 'Gene'
            ncRNA_description['name'] = 'tmRNA'
            
            #What about the alignment??????????
            db['annotations'].remove(tmRNA)
            db['ncRNAs'].insert(ncRNA_description)
        else:
            ncRNA_description['class'] = 'Overlapping RNA'
            ncRNA_description['name'] = 'Overlapping_RNA_from_NCBI'
            db['annotations'].remove(tmRNA)
            db['ncRNAs'].insert(ncRNA_description)                          
 
    for misc_RNA in db['annotations'].find({'class':'misc_RNA'}):
        ncRNA_description = {
            '_id': str(ObjectId()),            
        }

        for key in misc_RNA.keys():
             ncRNA_description[key] = misc_RNA[key]

        overlap = False

        for ncRNA in db['ncRNAs'].find({'genome':misc_RNA['genome']}):
            if ncRNA['genome'] == misc_RNA['genome'] and ncRNA['genomicStrand'] == misc_RNA['genomicStrand'] and (ncRNA['genomicPositions'][0] >= misc_RNA['genomicPositions'][0] and ncRNA['genomicPositions'][0] <= misc_RNA['genomicPositions'][1] or ncRNA['genomicPositions'][1] >= misc_RNA['genomicPositions'][0] and ncRNA['genomicPositions'][1] <= misc_RNA['genomicPositions'][1] or misc_RNA['genomicPositions'][0] >= ncRNA['genomicPositions'][0] and misc_RNA['genomicPositions'][0] <= ncRNA['genomicPositions'][1] or misc_RNA['genomicPositions'][1] >= ncRNA['genomicPositions'][0] and misc_RNA['genomicPositions'][1] <= ncRNA['genomicPositions'][1]):
                overlap = True
                break
        
        if not overlap: #its a new ncRNA, we add it to the ncRNAs table and remove it from the annotations table
            ncRNA_description['class'] = 'Misc RNA'
            ncRNA_description['name'] = 'Misc_RNA_from_NCBI'
            #What about the alignment??????????
            db['annotations'].remove(misc_RNA)
            db['ncRNAs'].insert(ncRNA_description)

        else:
            ncRNA_description['class'] = 'Overlapping RNA'
            ncRNA_description['name'] = 'Overlapping_RNA_from_NCBI'
            db['annotations'].remove(misc_RNA)
            db['ncRNAs'].insert(ncRNA_description)
    
    from pyrna.db import Rfam

    rfam = Rfam()
    rfamFamiliesDetails = rfam.get_families_details()
    haca_box_names = []
    cd_box_names = []
    snrna_splicing_names= []
    for index, row in rfamFamiliesDetails.iterrows():
        if row['family'] == "Gene, snRNA, snoRNA, HACA-box":
            haca_box_names.append(row['id'].upper())
        elif row['family'] == "Gene, snRNA, snoRNA, CD-box":
            cd_box_names.append(row['id'].upper())
        elif row['family'] == "Gene, snRNA, splicing":
            snrna_splicing_names.append(row['id'].upper())

    for annotation in db['annotations'].find({'class':'ncRNA'}):

        ncRNA_description = {
            '_id': str(ObjectId()),            
        }

        for key in annotation.keys():
             ncRNA_description[key] = annotation[key]

        overlap = False
        
        for ncRNA in db['ncRNAs'].find({'genome':annotation['genome']}):
            if ncRNA['genome'] == annotation['genome'] and ncRNA['genomicStrand'] == annotation['genomicStrand'] and (ncRNA['genomicPositions'][0] >= annotation['genomicPositions'][0] and ncRNA['genomicPositions'][0] <= annotation['genomicPositions'][1] or ncRNA['genomicPositions'][1] >= annotation['genomicPositions'][0] and ncRNA['genomicPositions'][1] <= annotation['genomicPositions'][1] or     annotation['genomicPositions'][0] >= ncRNA['genomicPositions'][0] and annotation['genomicPositions'][0] <= ncRNA['genomicPositions'][1] or annotation['genomicPositions'][1] >= ncRNA['genomicPositions'][0] and annotation['genomicPositions'][1] <= ncRNA['genomicPositions'][1]):
                overlap = True
                break
        
        if not overlap: #its a new ncRNA, we add it to the ncRNAs table and remove it from the annotations table
            ncRNA_class = annotation['ncRNA_class']
            if ncRNA_class == "antisense_RNA":
                ncRNA_description['class'] = 'Gene, antisense'
                #What about the alignment??????????
                #db['annotations'].remove(annotation)
                #db['ncRNAs'].insert(ncRNA_description)
            elif ncRNA_class == "autocatalytically_spliced_intron":
                #What about the alignment??????????
                #db['annotations'].remove(annotation)
                #db['ncRNAs'].insert(ncRNA_description)
                pass
            elif ncRNA_class == "telomerase_RNA":
                ncRNA_description['class'] = 'Gene'
                ncRNA_description['name'] = 'Sacc_telomerase'
                #What about the alignment??????????
                db['annotations'].remove(annotation)
                db['ncRNAs'].insert(ncRNA_description)
            elif ncRNA_class == "hammerhead_ribozyme":
                ncRNA_description['class'] = 'Gene, ribozyme'
                #What about the alignment??????????
                #db['annotations'].remove(annotation)
                #db['ncRNAs'].insert(ncRNA_description)
            elif ncRNA_class == "RNase_P_RNA":
                ncRNA_description['class'] = 'Gene, ribozyme'
                ncRNA_description['name'] = 'RNaseP_nuc' #could be RNaseP_bact_a or RNaseP_bact_b, who knows? RFAM found some RNaseP_bact_a in fungal genomes
                #What about the alignment??????????
                db['annotations'].remove(annotation)
                db['ncRNAs'].insert(ncRNA_description)
            elif ncRNA_class == "RNase_MRP_RNA":
                ncRNA_description['class'] = 'Gene, ribozyme'
                ncRNA_description['name'] = 'Rnase_MRP'
                #What about the alignment??????????
                db['annotations'].remove(annotation)
                db['ncRNAs'].insert(ncRNA_description)
            elif ncRNA_class == "guide_RNA":
                #What about the alignment??????????
                #db['annotations'].remove(annotation)
                #db['ncRNAs'].insert(ncRNA_description)
                pass
            elif ncRNA_class == "rasiRNA":
                ncRNA_description['class'] = 'rasiRNA'
                ncRNA_description['name'] = 'rasiRNA_from_NCBI'
                #What about the alignment??????????
                #db['annotations'].remove(annotation)
                #db['ncRNAs'].insert(ncRNA_description)
            elif ncRNA_class == "scRNA":
                ncRNA_description['class'] = 'scRNA'
                ncRNA_description['class'] = 'scRNA_from_NCBI'
                #What about the alignment??????????
                #db['annotations'].remove(annotation)
                #db['ncRNAs'].insert(ncRNA_description)
            elif ncRNA_class == "siRNA":
                ncRNA_description['class'] = 'siRNA'
                ncRNA_description['class'] = 'siRNA_from_NCBI'
                #What about the alignment??????????
                #db['annotations'].remove(annotation)
                #db['ncRNAs'].insert(ncRNA_description)
            elif ncRNA_class == "miRNA":
                ncRNA_description['class'] = 'Gene, miRNA'
                #What about the alignment??????????
                db['annotations'].remove(annotation)
                db['ncRNAs'].insert(ncRNA_description)
            elif ncRNA_class == "snoRNA":
                description = ""
                if annotation.has_key('note'):
                    description += annotation['note']+" "

                if annotation.has_key('product'):
                    description += annotation['product']

                name1, name2, name3 = None, None, None

                for token in description.split(' '):
                    if re.findall("snr[0-9]+", token, re.IGNORECASE): #we search for SNRXX word
                        name1 = token.upper()
                        name2 = "SNO%s"%name1 #SNRXX => SNOSNRXX
                        name3 = "SNO%s"%name1.split('SN')[-1] #SNRXX => SNORXX
                        break

                if re.findall("H/ACA", description):
                    ncRNA_description['class'] = 'Gene, snRNA, snoRNA, HACA-box'
                    found = False
                    for index, row in rfamFamiliesDetails.iterrows():
                        if row['family'] == "Gene, snRNA, snoRNA, HACA-box" and (row['id'].upper() == name1 or row['id'].upper() == name2 or row['id'].upper() == name3):
                            found = True
                            ncRNA_description['name'] = row['id']
                            #What about the alignment??????????
                            db['annotations'].remove(annotation)
                            db['ncRNAs'].insert(ncRNA_description)
                            break
                    if not found:
                        ncRNA_description['class'] = 'Misc RNA'
                        ncRNA_description['name'] = 'Misc_RNA_from_NCBI'
                        db['annotations'].remove(annotation)
                        db['ncRNAs'].insert(ncRNA_description)
                    
                elif re.findall("C/D", description):
                    ncRNA_description['class'] = 'Gene, snRNA, snoRNA, CD-box'
                    found = False
                    for index, row in rfamFamiliesDetails.iterrows():
                        if row['family'] == "Gene, snRNA, snoRNA, CD-box" and (row['id'].upper() == name1 or row['id'].upper() == name2 or row['id'].upper() == name3):
                            found = True
                            ncRNA_description['name'] = row['id']
                            #What about the alignment??????????
                            db['annotations'].remove(annotation)
                            db['ncRNAs'].insert(ncRNA_description)
                            break
                    if not found:
                        ncRNA_description['class'] = 'Misc RNA'
                        ncRNA_description['name'] = 'Misc_RNA_from_NCBI'
                        db['annotations'].remove(annotation)
                        db['ncRNAs'].insert(ncRNA_description) 
                else:
                    found = False
                    for index, row in rfamFamiliesDetails.iterrows():
                        if row['family'] == "Gene, snRNA, snoRNA, HACA-box" and (row['id'].upper() == name1 or row['id'].upper() == name2 or row['id'].upper() == name3):
                            found = True
                            ncRNA_description['class'] = 'Gene, snRNA, snoRNA, HACA-box'
                            ncRNA_description['name'] = row['id']
                            #What about the alignment??????????
                            db['annotations'].remove(annotation)
                            db['ncRNAs'].insert(ncRNA_description)
                            break 

                    for index, row in rfamFamiliesDetails.iterrows():
                        if row['family'] == "Gene, snRNA, snoRNA, CD-box" and (row['id'].upper() == name1 or row['id'].upper() == name2 or row['id'].upper() == name3):
                            found = True
                            ncRNA_description ['_id'] =  str(ObjectId()) #without any explicit C/D or H/ACA word in the description, the RNA_description could be added twice. Then new _id.            
                            ncRNA_description['class'] = 'Gene, snRNA, snoRNA, CD-box'
                            ncRNA_description['name'] = row['id']
                            #What about the alignment??????????
                            db['annotations'].remove(annotation)
                            db['ncRNAs'].insert(ncRNA_description)
                            break

                    if not found:
                        ncRNA_description['class'] = 'Misc RNA'
                        ncRNA_description['name'] = 'Misc_RNA_from_NCBI'
                        db['annotations'].remove(annotation)
                        db['ncRNAs'].insert(ncRNA_description)

            elif ncRNA_class == "snRNA":
                description = ""
                if annotation.has_key('note'):
                    description += annotation['note']+" "

                if annotation.has_key('product'):
                    description += annotation['product']

                matchObj = re.search("(U[1-9]).+(spliceosomal)?.+(sn)?RNA",description, re.IGNORECASE)

                name = None
                if matchObj:
                    name =  matchObj.group(1)
                
                found = False
                for index, row in rfamFamiliesDetails.iterrows():
                    if row['family'] == "Gene, snRNA, splicing" and row['id'] == name:
                        found = True
                        ncRNA_description['class'] = 'Gene, snRNA, splicing'
                        ncRNA_description['name'] = row['id']
                        #What about the alignment??????????
                        db['annotations'].remove(annotation)
                        db['ncRNAs'].insert(ncRNA_description)
                        break

                if not found:
                    ncRNA_description['class'] = 'Misc RNA'
                    ncRNA_description['name'] = 'Misc_RNA_from_NCBI'
                    db['annotations'].remove(annotation)
                    db['ncRNAs'].insert(ncRNA_description)

            elif ncRNA_class == "SRP_RNA":
                ncRNA_description['class'] = 'Gene'
                ncRNA_description['name'] = 'Fungi_SRP'
                #What about the alignment??????????
                db['annotations'].remove(annotation)
                db['ncRNAs'].insert(ncRNA_description)
            elif ncRNA_class == "vault_RNA":
                ncRNA_description['class'] = 'Gene'
                ncRNA_description['name'] = 'Vault'
                #What about the alignment??????????
                db['annotations'].remove(annotation)
                db['ncRNAs'].insert(ncRNA_description)
            elif ncRNA_class == "Y_RNA":
                ncRNA_description['class'] = 'Gene'
                ncRNA_description['name'] = 'Y_RNA'
                #What about the alignment??????????
                db['annotations'].remove(annotation)
                db['ncRNAs'].insert(ncRNA_description)
            elif ncRNA_class == "other": #other RNA are not migrated, they have to be checked manually
                ncRNA_description['class'] = 'Misc RNA'
                ncRNA_description['name'] = 'Misc_RNA_from_NCBI'
                db['annotations'].remove(annotation)
                db['ncRNAs'].insert(ncRNA_description)
            else: # RNA with unknown ncRNA_class are not migrated, they have to be checked manually
                ncRNA_description['class'] = 'Misc RNA'
                ncRNA_description['name'] = 'Misc_RNA_from_NCBI'
                db['annotations'].remove(annotation)
                db['ncRNAs'].insert(ncRNA_description)

        else:
            ncRNA_description['class'] = 'Overlapping RNA'
            ncRNA_description['name'] = 'Overlapping_RNA_from_NCBI'
            db['annotations'].remove(annotation)
            db['ncRNAs'].insert(ncRNA_description)

    print "%i ncRNAs added."%(db['ncRNAs'].count()-before_ncRNAs)
    print "%i annotations removed "%(before_annotations-db['annotations'].count())

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
    o = open(pdb_file,'r')
    for line in o.readlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            new_lines.append("%s%5s%s" % (line[0:6],counter,line[11:]))
            counter += 1
        else:
            new_lines.append(line)

    return ''.join(new_lines)



