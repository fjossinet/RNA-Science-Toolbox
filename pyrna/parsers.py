import re
from pandas import DataFrame
from pyrna.features import RNA, DNA, TertiaryStructure, SecondaryStructure
from pyrna import utils
from pysam import Samfile

def consensus2d_to_base_pairs(aligned_rna, consensus_2d):
    """
    Parameters:
    ---------
    - aligned_rna: an RNA object (see pyrna.features) whose sequence can contain gap symbols ('-')
    - consensus_2d: the consensus secondary structure and described as a list of base_pairs in a pandas Dataframe

    Returns:
    ------
    the secondary structure as a list of base-pairs in a pandas DataFrame
    """

    aligned_sequence = aligned_rna.sequence
    rna = RNA(name = aligned_rna.name, sequence = aligned_sequence.replace('-',''))
    new_base_pairs = []

    for index, bp in consensus_2d.iterrows():
        if aligned_sequence[bp['pos1']-1] != '-' and aligned_sequence[bp['pos2']-1] != '-':
            new_pos1 = bp['pos1']-aligned_sequence[0:bp['pos1']].count('-')
            new_pos2 = bp['pos2']-aligned_sequence[0:bp['pos2']].count('-')
            new_base_pairs.append({
                'pos1': new_pos1,
                'pos2': new_pos2,
                'orientation': bp['orientation'],
                'edge1': bp['edge1'],
                'edge2': bp['edge2']
                })

    return DataFrame(new_base_pairs)

def secondary_structure_to_base_pairs(secondary_structure, keep_tertiaries = False):
    """
    Parameters:
    ---------
    - secondary_structure: a SecondaryStructure object (see pyrna.features)
    - keep_tertiaries (default: False) : if True, the tertiary interactions will be exported.

    Returns:
    ------
    the secondary structure as a list of base-pairs in a pandas DataFrame
    """
    base_pairs = []
    for helix in secondary_structure.helices:
        length = helix['length']
        helix_start = helix['location'][0][0]
        helix_end = helix['location'][-1][-1]

        for i in range(0,length):
            found_non_canonical = False
            for interaction in helix['interactions']:
                if interaction['location'][0][0] == helix_start+i and interaction['location'][-1][-1] == helix_end-i:
                    base_pairs.append({
                        'pos1': interaction['location'][0][0],
                        'pos2': interaction['location'][-1][-1],
                        'orientation': interaction['orientation'],
                        'edge1': interaction['edge1'],
                        'edge2': interaction['edge2']
                    })
                    found_non_canonical = True
                    break

            if not found_non_canonical:
                base_pairs.append({
                    'pos1': helix_start+i,
                    'pos2': helix_end-i,
                    'orientation': 'c',
                    'edge1': '(',
                    'edge2': ')'
                })    


    if keep_tertiaries:
        for interaction in secondary_structure.tertiary_interactions:
            base_pairs.append({
                'pos1': interaction['location'][0][0],
                'pos2': interaction['location'][-1][-1],
                'orientation': interaction['orientation'],
                'edge1': interaction['edge1'],
                'edge2': interaction['edge2']
            })         

    return DataFrame(base_pairs)

def base_pairs_to_secondary_structure(rna, base_pairs):
    """
    Parameters:
    ---------
    - rna: an RNA object (see pyrna.features)
    - base_pairs: the base pairs listed in a pandas Dataframe

    Returns:
    ------
    a SecondaryStructure object (see pyrna.features)
    """

    ss = SecondaryStructure(rna)

    if not len(base_pairs):
        ss.add_single_strand("SS1", 1, len(rna))
        return ss

    new_helix = False
    helix_start = -1
    helix_end = -1
    helix_length = -1
    base_pairs = base_pairs.sort_index(by='pos1') #the base pairs are sorted according to the first position
    base_pairs = base_pairs.as_matrix()
    helix_count = 1

    next_pos1 = -1
    next_edge1 = None
    next_pos2 = -1
    next_edge2 = None
    next_orientation = None

    non_canonical_secondary_interactions = []

    for i in range(0, len(base_pairs)-1):
        bp = base_pairs[i]
        pos1 = bp[3]
        edge1 = bp[1]
        pos2 = bp[4]
        edge2 = bp[2]
        orientation = bp[0]

        next_bp = base_pairs[i+1]
        next_pos1 = next_bp[3]
        next_edge1 = next_bp[1]
        next_pos2 = next_bp[4]
        next_edge2 = next_bp[2]
        next_orientation = next_bp[0]

        if pos1+1 == next_pos1 and pos2-1 == next_pos2:
            if new_helix:
                helix_length += 1
                if not utils.is_canonical(rna[next_pos1-1], rna[next_pos2-1], next_orientation, next_edge1, next_edge2):
                    non_canonical_secondary_interactions.append((next_orientation, next_edge1, next_edge2, next_pos1, next_pos2))
            else:
                new_helix = True
                helix_length = 2
                helix_start = pos1
                helix_end = pos2
                if not utils.is_canonical(rna[pos1-1], rna[pos2-1], orientation, edge1, edge2):
                    non_canonical_secondary_interactions.append((orientation, edge1, edge2, pos1, pos2))
                if not utils.is_canonical(rna[next_pos1-1], rna[next_pos2-1], next_orientation, next_edge1, next_edge2):
                    non_canonical_secondary_interactions.append((next_orientation, next_edge1, next_edge2, next_pos1, next_pos2))
                        
        else:
            if new_helix:
                ss.add_helix("H"+str(helix_count), helix_start, helix_end, helix_length)
                helix_count += 1
            else:
                ss.add_tertiary_interaction(orientation, edge1, edge2, pos1, pos2)
            new_helix = False

    #the last helix
    if new_helix:
        ss.add_helix("H"+str(helix_count), helix_start, helix_end, helix_length)
        helix_count+=1
    else:
        ss.add_tertiary_interaction(next_orientation, next_edge1, next_edge2, next_pos1, next_pos2)

    #now we add the non-canonical interactions to their helices
    for non_canonical_secondary_interaction in non_canonical_secondary_interactions:
        ss.add_base_pair( non_canonical_secondary_interaction[0], non_canonical_secondary_interaction[1], non_canonical_secondary_interaction[2], non_canonical_secondary_interaction[3], non_canonical_secondary_interaction[4] )

    #we construct the single-strands
    ss_count = 1
    ss_start = -1
    ss_length = 0
    for i in range(1, len(rna)+1):
        paired_residue = ss.get_paired_residue(i)
        if paired_residue != -1 and ss_length > 0:
            ss.add_single_strand("SS"+str(ss_count), ss_start, ss_length)
            ss_length = 0
            ss_count += 1
        elif paired_residue == -1:
            if ss_length == 0:
                ss_start = i
            ss_length += 1

    #the last single-strand
    if ss_length > 0:
        ss.add_single_strand("SS"+str(ss_count), ss_start, ss_length)

    return ss

def to_pdb(tertiary_structure, export_numbering_system = False):
    """
    Convert a TertiaryStructure object into PDB data

    Parameters:
    ---------
    - tertiary_structure: a TertiaryStructure object (see pyrna.features)
    - export_numbering_system (default: False): export the numbering system. If False, the residues are numbered from 1 to the length of the molecular chain

    Returns:
    ------
    the PDB data as a String
    """
    lines= []
    i = 1
    keys = []

    for k in tertiary_structure.residues.keys():
        keys.append(k)

    keys.sort() #the absolute position are sorted

    for key in keys:
        atoms = tertiary_structure.residues[key]['atoms']
        for atom in atoms:
            if export_numbering_system:
                lines.append("%-6s%5u  %-4s%3s %s%4s    %8.3f%8.3f%8.3f"%("ATOM", i, atom['name'], tertiary_structure.rna.sequence[key-1], tertiary_structure.rna.name[0], tertiary_structure.get_residue_label(key), atom['coords'][0], atom['coords'][1], atom['coords'][2]))
            else:
                lines.append("%-6s%5u  %-4s%3s %s%4u    %8.3f%8.3f%8.3f"%("ATOM", i, atom['name'], tertiary_structure.rna.sequence[key-1], tertiary_structure.rna.name[0], key, atom['coords'][0], atom['coords'][1], atom['coords'][2]))
            i += 1

    lines.append("END")

    return '\n'.join(lines)

def to_ct(base_pairs, rna):
    """
    Convert a list of base pairs into CT data

    Parameters:
    ---------
    - base_pairs: the base pairs listed in a pandas Dataframe
    - rna : an RNA object (see pyrna.features)

    Returns:
    ------
    the CT data as a String
    """
    lines=["ENERGY"]
    base_pairs = base_pairs.sort_index(by='pos1') #the base pairs are sorted according to the first position
    for molecular_pos in range (1, len(rna)+1):
        if len(base_pairs[base_pairs.pos1 == molecular_pos]):
            for index, bps in base_pairs[base_pairs.pos1 == molecular_pos].iterrows():
                lines.append("%i\t%s\t%i\t%i\t%i\t%i"%(molecular_pos, rna.sequence[molecular_pos-1], molecular_pos-1, molecular_pos+1, bps['pos2'], molecular_pos))
        elif len(base_pairs[base_pairs.pos2 == molecular_pos]):
            for index, bps in base_pairs[base_pairs.pos2 == molecular_pos].iterrows():
                lines.append("%i\t%s\t%i\t%i\t%i\t%i"%(molecular_pos, rna.sequence[molecular_pos-1], molecular_pos-1, molecular_pos+1, bps['pos1'], molecular_pos))
        else:
            lines.append("%i\t%s\t%i\t%i\t0\t%i"%(molecular_pos, rna.sequence[molecular_pos-1], molecular_pos-1, molecular_pos+1, molecular_pos))
        
    return '\n'.join(lines)    

def to_fasta(molecules, single_line=False):
    """
    Convert a list of Molecule objects into FASTA data

    Parameters:
    ---------
    - molecules: a list of Molecule objects (see pyrna.features)
    - single_line (default: False): if True, each molecular sequence will we exported into a single line

    Returns:
    ------
    the FASTA data as a String
    """
    outputs = []
    for molecule in molecules:
        outputs.append(molecule.to_fasta(single_line))
    return '\n'.join(outputs)

def to_vienna(base_pairs, molecules, single_line=False):
    """
    Convert lists of base pairs and Molecule objects into Vienna data.

    Parameters:
    ---------
    - base_pairs: a list of pandas Dataframes, each one listing base pairs. This list should contain either a single Dataframe or as many Dataframes as Molecule objects provided in the second parameter.
    - molecules: a list of Molecule objects (gapped or ungapped) (see pyrna.features)
    - single_line (default: False): if True, each molecular sequence will we exported into a single line

    Returns:
    ------
    the Vienna data as a String
    """
    if len(base_pairs) != 1 and len(molecules) != len(base_pairs):
        raise Exception("You need to provide either a single Dataframe or as many Dataframes as Molecule objects")
    
    if single_line:
        lines =[]
        if len(base_pairs) == 1:
            for molecule in molecules:
                lines.append(">"+molecule.name)
                lines.append(molecule.sequence)
            lines.append(to_bn(base_pairs[0], len(molecules[0])))
        else:
            for i in range(0, len(molecules)):
                lines.append(">"+molecules[i].name)
                lines.append(molecules[i].sequence)
                lines.append(to_bn(base_pairs[i], len(molecules[i])))       
        return '\n'.join(lines)
    else:
        if len(base_pairs) == 1:
            sequence_lines = []
            bn_lines = []
            bn = to_bn(base_pairs[0], len(molecules[0]))
            for molecule in molecules:
                c = 0
                sequence_lines.append(">"+molecule.name)
                while c < len(molecule):
                    d = min(len(molecule), c + 79)
                    sequence_lines.append(molecule[c:d])
                    c += 79
            c = 0
            while c < len(molecules[0]):
                d = min(len(molecules[0]), c + 79)
                bn_lines.append(bn[c:d])
                c += 79
            return '\n'.join(sequence_lines)+'\n'+'\n'.join(bn_lines)
        else:
            lines = []
            for i in range(0, len(molecules)):
                bn = to_bn(base_pairs[i], len(molecules[i]))

                c = 0
                lines.append(">"+molecules[i].name)
                while c < len(molecules[i]):
                    d = min(len(molecules[i]), c + 79)
                    lines.append(molecules[i][c:d])
                    c += 79

                c = 0
                while c < len(bn):
                    d = min(len(bn), c + 79)
                    lines.append(bn[c:d])
                    c += 79
            return '\n'.join(lines)

def to_stockholm(base_pairs, molecules, rfam_accession_number = None, family_id = None):
    """
    Convert a list of base pairs and a list of Molecule objects into Stockholm data. 

    Parameters:
    ---------
    - base_pairs: a pandas Dataframe listing the base pairs.
    - molecules: a list of Molecule objects (gapped or ungapped) (see pyrna.features)
    - rfam_accession_number (default: None): the RFAM ID to export (corresponding to the line starting with #=GF AC)
    - family_id (default: None): the family ID to export (corresponding to the line starting with #=GF ID)

    Returns:
    ------
    the Stockholm data as a String
    """
    bn = to_bn(base_pairs, len(molecules[0]))
    lines = []
    lines.append("# STOCKHOLM 1.0")
    if rfam_accession_number:
        lines.append("#=GF AC %s"%rfam_accession_number)
    if family_id:
        lines.append("#=GF ID  %s"%family_id)
    c = 0
    alignment_length = len(bn)

    while c < alignment_length:
        d = min(alignment_length, c + 80)
        for molecule in molecules:
            lines.append("%s\t%s"%(molecule.name, molecule.sequence[c:d]))
        lines.append("#=GC SS_cons\t%s"%bn[c:d])
        lines.append("")
        c += 80

    lines.append("//")
    return '\n'.join(lines)

def to_clustalw(base_pairs, molecules, curate = False):
    """
    Convert a list of base pairs and a list of Molecule objects into Clustalw data. 
    Parameters:
    ---------
    - base_pairs: a pandas Dataframe listing the base-pairs. This can be a consensus 2D.
    - molecules: an list of Molecule objects (gapped or ungapped) (see pyrna.features)
    - curate (default: False): remove the columns filled with gaps
    
    Returns:
    ------
    the clustalw data as a String. The name of the molecules will be non-redundant and will not contain any spaces characters.
    """

    sequence_lines = []
    bn = to_bn(base_pairs, len(molecules[0]))
    c = 0
    gaps_positions = None
    while c < len(molecules[0]):
        names = []
        for molecule in molecules:
            if not gaps_positions:
                gaps_positions = set(molecule.get_gaps_positions())
            else: #we search the same gap positions
                gaps_positions = gaps_positions.intersection(molecule.get_gaps_positions())
            d = min(len(molecule), c + 60)
            name = molecule.name.replace(' ', '_') #molecule name without any space
            if names.count(name): #and non-redundant
                sequence_lines.append('%s.%i'%(name,names.count(name))+"\t"+molecule[c:d]+'\n')
            else:
                sequence_lines.append(name+"\t"+molecule[c:d]+'\n')#if already non-redundant, not .0 as suffix
            names.append(name)
        sequence_lines.append('\n')
        c += 60

    if curate:
        sequence_lines = []
        chars = list(molecules[0].sequence)
        for pos in gaps_positions:
            chars[pos]=""
        curated_sequence = ''.join(chars)

        c = 0
        while c < len(curated_sequence):
            names = []
            for molecule in molecules:
                chars = list(molecule.sequence)
                for pos in gaps_positions:
                    chars[pos]=""
                curated_sequence = ''.join(chars)
                d = min(len(curated_sequence), c + 60)
                name = molecule.name.replace(' ', '_') #molecule name without any space
                if names.count(name): #and non-redundant
                    sequence_lines.append('%s.%i'%(name,names.count(name))+"\t"+curated_sequence[c:d]+'\n')
                else:
                    sequence_lines.append(name+"\t"+curated_sequence[c:d]+'\n') #if already non-redundant, not .0 as suffix
                names.append(name)
            sequence_lines.append('\n')
            c += 60

        chars = list(bn)

        i = 0
        left_positions = []
        base_pairs = []
        for char in chars:
            if char == '(':
                left_positions.append(i)
            elif char == ')':
                base_pairs.append([left_positions.pop(), i])
            i+=1

        for base_pair in base_pairs:
            if base_pair[0] in gaps_positions or base_pair[1] in gaps_positions: #this means that one of the paired column will be removed, then the partner column becomes a single-strand (otherwise the curated bracket notation will be unbalanced)
                chars[base_pair[0]] = '.'
                chars[base_pair[1]] = '.'

        for pos in gaps_positions:
            chars[pos]=""

        bn = ''.join(chars)

    return ''.join(sequence_lines)+"2D\t"+bn

def to_bn(base_pairs, length):
    """
    Convert a list of base pairs into a bracket notation

    Parameters:
    ---------
    - base_pairs: a pandas Dataframe listing the base pairs.
    - length: the length of the molecule linked to the secondary structure. This is mandatory, since the base_pairs parameter describes only paired positions in the sequence.
    
    Returns:
    ------
    the bracket notation as a String
    """
    bn = []
    if not len(base_pairs): #its a pure single_strand fold
        for pos in range(1, length+1):
            bn.append('.')
    else:
        for pos in range(1, length+1):
            if len(base_pairs.pos1[base_pairs.pos1 == pos]) != 0:
                bn.append(base_pairs[base_pairs.pos1 == pos].edge1.values[0])
            elif len(base_pairs.pos2[base_pairs.pos2 == pos]) != 0:
                bn.append(base_pairs[base_pairs.pos2 == pos].edge2.values[0])
            else:
                bn.append('.')

    return ''.join(bn)

def parse_genbank(genbank_data):
    """
    Parse Genbank data.

   Parameters:
   ---------
    - genbank_data: the Genbank data as a String

    Returns:
    ------
    a DNA object (see pyrna.features) and a pandas Dataframe listing the genomic features. The columns are:
    - type
    - genomicStrand ('+' or '-')
    - genomicPositions
    - a sequence if it is a ncRNA
    - a column for each qualifier attached to the feature (/qualifier=)
    """

    pieces_of_seq=[]
    start_of_sequence = False
    accession = None
    feature_type = None
    qualifer_type = None
    qualifier_content = None
    qualifiers = []
    genomic_strand = '+'
    genomic_positions = None
    features = []
    organism = None
    inOrganism = False
    lineage = ""
    location = None
    lines = genbank_data.strip().split('\n')
    if not lines[-1].strip() == '//':
        raise Exception("Uncomplete file")
    for line in lines:
        tokens = re.split('\s+', line)
        if line.startswith('ACCESSION'):
            accession = re.split('\s+', line)[1]
        elif line.strip().startswith('ORGANISM'):
            organism = line.strip().split('ORGANISM')[1].strip()
            inOrganism = True
        elif line.strip().startswith('REFERENCE'):
            inOrganism = False
        elif line.startswith('ORIGIN'): #the genomic sequence
            start_of_sequence = True
            #we store the last feature
            #the last
            if feature_type and not feature_type == "source" and not feature_type == "intron":
                if location.startswith("complement(join("):#a joined location on the Crick strand
                    genomic_strand = '-'
                    ends = location.split('complement(join(')[1][:-2].replace("join(", "").replace(')','').replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                    if feature_type == 'CDS':
                        for i in range(1, len(ends)-2, 2):
                            intron = {
                                'type': 'intron',
                                'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                'genomicStrand': genomic_strand,
                            }
                            features.append(intron)
                elif location.startswith("join(complement("):#a joined location on the Crick strand
                    genomic_strand = '-'
                    ends = location.split('join(complement(')[1][:-2].replace("complement(", "").replace(')', '').replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                    if feature_type == 'CDS':
                        for i in range(1, len(ends)-2, 2):
                            intron = {
                                'type': 'intron',
                                'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                'genomicStrand': genomic_strand,
                            }
                            features.append(intron)
                elif location.startswith("complement(order("):
                    genomic_strand = '-'
                    ends = location.split('complement(order(')[1][:-2].replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                elif location.startswith("order("):
                    ends = location.split('order(')[1][:-1].replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                elif location.startswith("complement("): #a location on the Crick strand
                    genomic_strand = '-'
                    ends = location.split('complement(')[1][:-1].split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                elif location.startswith("join("): #a joined location
                    ends = location.split('join(')[1][:-1].replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                    if feature_type == 'CDS':
                        for i in range(1, len(ends)-2, 2):
                            intron = {
                                'type': 'intron',
                                'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                'genomicStrand': genomic_strand,
                            }
                            features.append(intron)
                else: #a regular location
                    ends = location.split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]

                feature = {
                    'type': feature_type,
                    'genomicPositions': genomic_positions,
                    'genomicStrand': genomic_strand,
                }
                if qualifer_type and qualifier_content:
                    if qualifer_type == 'translation':
                        qualifier_content = qualifier_content.replace(" ","")
                    qualifiers.append({
                        "type": qualifer_type,
                        "content": qualifier_content
                         })
                for qualifier in qualifiers:
                    feature[qualifier['type']] = qualifier['content']
                features.append(feature)
        elif len(tokens) == 3 and re.findall('\.\.>?[0-9]+', tokens[2]):
            #new feature
            #we store the previous one (if any)
            if feature_type and not feature_type == "source" and not feature_type == "intron":
                if location.startswith("complement(join("):#a joined location on the Crick strand
                    genomic_strand = '-'
                    ends = location.split('complement(join(')[1][:-2].replace("join(", "").replace(')','').replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                    if feature_type == 'CDS':
                        for i in range(1, len(ends)-2, 2):
                            intron = {
                                'type': 'intron',
                                'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                'genomicStrand': genomic_strand,
                            }
                            features.append(intron)
                elif location.startswith("join(complement("):#a joined location on the Crick strand
                    genomic_strand = '-'
                    ends = location.split('join(complement(')[1][:-2].replace("complement(", "").replace(')', '').replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                    if feature_type == 'CDS':
                        for i in range(1, len(ends)-2, 2):
                            intron = {
                                'type': 'intron',
                                'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                'genomicStrand': genomic_strand,
                            }
                            features.append(intron)
                elif location.startswith("complement(order("):
                    genomic_strand = '-'
                    ends = location.split('complement(order(')[1][:-2].replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                elif location.startswith("order("):
                    ends = location.split('order(')[1][:-1].replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                elif location.startswith("complement("): #a location on the Crick strand
                    genomic_strand = '-'
                    ends = location.split('complement(')[1][:-1].split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                elif location.startswith("join("): #a joined location
                    ends = location.split('join(')[1][:-1].replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                    if feature_type == 'CDS':
                        for i in range(1, len(ends)-2, 2):
                            intron = {
                                'type': 'intron',
                                'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                'genomicStrand': genomic_strand,
                            }
                            features.append(intron)
                else: #a regular location
                    ends = location.split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]

                feature = {
                    'type': feature_type,
                    'genomicPositions': genomic_positions,
                    'genomicStrand': genomic_strand,
                }

                if qualifer_type and qualifier_content:
                    if qualifer_type == 'translation':
                        qualifier_content = qualifier_content.replace(" ","")
                    qualifiers.append({
                        "type": qualifer_type,
                        "content": qualifier_content
                         })
                for qualifier in qualifiers:
                    feature[qualifier['type']] = qualifier['content']
                features.append(feature)
            feature_type = None
            genomic_strand = '+'
            genomic_positions = None
            qualifer_type = None
            qualifier_content = None
            qualifiers = []
            feature_type = tokens[1].strip()
            location = tokens[2].strip()
        elif not qualifer_type and not qualifier_content and len(tokens) == 2 and re.findall('\.\.',tokens[1]): #still the content of the current location
            location += tokens[1].strip()
        elif re.findall('^\s+/.+=', line): # a new qualifier /bla_bla=
            if qualifer_type and qualifier_content:
                if qualifer_type == 'translation':
                    qualifier_content = qualifier_content.replace(" ","")
                qualifiers.append({
                    "type": qualifer_type,
                    "content": qualifier_content
                     })
            qualifer_type = line.strip()[1:].split('=')[0].strip()[0:]
            qualifier_content = line.strip()[1:].split('=')[1].strip().replace('"','')
        elif re.findall('^\s+/.+', line): # a qualifier like /manual => ignore
            pass
        elif not start_of_sequence and qualifer_type and qualifier_content : #still the content of the current qualifier
            qualifier_content += " "+line.strip().replace('"','')
        elif line.startswith('//'): #end of the genomic sequence
            start_of_sequence = False
        elif start_of_sequence:
            pieces_of_seq.append(''.join(re.split('\s+',line.strip())[1:]).upper())
        elif inOrganism:
            lineage += " "+line.strip()

    dna = DNA( name = accession, sequence=''.join(pieces_of_seq))

    dna.lineage = lineage.strip()

    if organism:
        dna.organism = organism

    for feature in features:
        if feature.has_key('ncRNA_class'):
            if feature['genomicStrand'] == '+':
                feature['sequence'] = dna.sequence[feature['genomicPositions'][0]-1:feature['genomicPositions'][-1]]
            else:
                feature['sequence'] = DNA(name = dna.name, sequence = dna.sequence[feature['genomicPositions'][0]-1:feature['genomicPositions'][-1]]).get_complement()[::-1]

    return dna, DataFrame(features)

def parse_embl(embl_data):
    """
    Parse EMBL data.

    Parameters:
   ---------
     - embl_data:  EMBL data as a String

    Returns:
    -------
    a DNA object (see pyrna.features) and a pandas Dataframe listing the genomic features. The columns are:
    - feature type
    - genomicStrand ('+' or '-')
    - genomicPositions
    - a sequence if it is a ncRNA
    - a column for each qualifier attached to the feature (/qualifier=)
    """

    pieces_of_seq=[]
    start_of_sequence = False
    accession = None
    feature_type = None
    qualifer_type = None
    qualifier_content = None
    qualifiers = []
    genomic_strand = '+'
    genomic_positions = None
    features = []
    organism = None
    location = None
    lines = embl_data.strip().split('\n')
    if not lines[-1].strip() == '//':
        raise Exception("Uncomplete file")
    for line in lines:
        tokens = re.split('\s+', line)
        if line.startswith('AC'):
            accession = re.split('\s+', line)[-1]
        elif line.strip().startswith('OS'):
            organism = line.strip().split('OS')[1].strip()
        elif line.startswith('SQ'): #the genomic sequence
            start_of_sequence = True
            #we store the last feature
            if feature_type and not feature_type == "source" and not feature_type == "intron":
                if location.startswith("complement(join("):#a joined location on the Crick strand
                    genomic_strand = '-'
                    ends = location.split('complement(join(')[1][:-2].replace("join(", "").replace(')','').replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                    if feature_type == 'CDS':
                        for i in range(1, len(ends)-2, 2):
                            intron = {
                                'type': 'intron',
                                'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                'genomicStrand': genomic_strand,
                            }
                            features.append(intron)
                elif location.startswith("join(complement("):#a joined location on the Crick strand
                    genomic_strand = '-'
                    ends = location.split('join(complement(')[1][:-2].replace("complement(", "").replace(')', '').replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                    if feature_type == 'CDS':
                        for i in range(1, len(ends)-2, 2):
                            intron = {
                                'type': 'intron',
                                'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                'genomicStrand': genomic_strand,
                            }
                            features.append(intron)
                elif location.startswith("complement(order("):
                    genomic_strand = '-'
                    ends = location.split('complement(order(')[1][:-2].replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                elif location.startswith("order("):
                    ends = location.split('order(')[1][:-1].replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                elif location.startswith("complement("): #a location on the Crick strand
                    genomic_strand = '-'
                    ends = location.split('complement(')[1][:-1].split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                elif location.startswith("join("): #a joined location
                    ends = location.split('join(')[1][:-1].replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                    if feature_type == 'CDS':
                        for i in range(1, len(ends)-2, 2):
                            intron = {
                                'type': 'intron',
                                'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                'genomicStrand': genomic_strand,
                            }
                            features.append(intron)
                else: #a regular location
                    ends = location.split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                feature = {
                    'type': feature_type,
                    'genomicPositions': genomic_positions,
                    'genomicStrand': genomic_strand,
                }
                if qualifer_type and qualifier_content:
                    if qualifer_type == 'translation':
                        qualifier_content = qualifier_content.replace(" ","")
                    qualifiers.append({
                        "type": qualifer_type,
                        "content": qualifier_content
                         })
                for qualifier in qualifiers:
                    feature[qualifier['type']] = qualifier['content']
                features.append(feature)
        elif re.findall('^FT\s+/.+=', line): # a new qualifier /bla_bla=
            if qualifer_type and qualifier_content:
                if qualifer_type == 'translation':
                    qualifier_content = qualifier_content.replace(" ","")
                qualifiers.append({
                    "type": qualifer_type,
                    "content": qualifier_content
                     })
            qualifer_type = line.strip()[2:].split('=')[0].strip()[1:]
            qualifier_content = line.strip()[2:].split('=')[1].strip().replace('"','')
        elif re.findall('^FT\s+/.+', line): # a qualifier like /manual => ignore
            pass
        elif len(tokens) == 3 and re.findall('\.\.>?[0-9]+',tokens[2]):
            #new feature
            #we store the previous one (if any)
            if feature_type and not feature_type == "source" and not feature_type == "intron" :
                if location.startswith("complement(join("):#a joined location on the Crick strand
                    genomic_strand = '-'
                    ends = location.split('complement(join(')[1][:-2].replace("join(", "").replace(')','').replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                    if feature_type == 'CDS':
                        for i in range(1, len(ends)-2, 2):
                            intron = {
                                'type': 'intron',
                                'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                'genomicStrand': genomic_strand,
                            }
                            features.append(intron)
                elif location.startswith("join(complement("):#a joined location on the Crick strand
                    genomic_strand = '-'
                    ends = location.split('join(complement(')[1][:-2].replace("complement(", "").replace(')', '').replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                    if feature_type == 'CDS':
                        for i in range(1, len(ends)-2, 2):
                            intron = {
                                'type': 'intron',
                                'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                'genomicStrand': genomic_strand,
                            }
                            features.append(intron)
                elif location.startswith("complement(order("):
                    genomic_strand = '-'
                    ends = location.split('complement(order(')[1][:-2].replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                elif location.startswith("order("):
                    ends = location.split('order(')[1][:-1].replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                elif location.startswith("complement("): #a location on the Crick strand
                    genomic_strand = '-'
                    ends = location.split('complement(')[1][:-1].split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                elif location.startswith("join("): #a joined location
                    ends = location.split('join(')[1][:-1].replace(',','..').split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]
                    if feature_type == 'CDS':
                        for i in range(1, len(ends)-2, 2):
                            intron = {
                                'type': 'intron',
                                'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                'genomicStrand': genomic_strand,
                            }
                            features.append(intron)
                else: #a regular location
                    ends = location.split('..')
                    ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                    genomic_positions = [min(ends), max(ends)]

                feature = {
                    'type': feature_type,
                    'genomicPositions': genomic_positions,
                    'genomicStrand': genomic_strand,
                }
                if qualifer_type and qualifier_content:
                    if qualifer_type == 'translation':
                        qualifier_content = qualifier_content.replace(" ","")
                    qualifiers.append({
                        "type": qualifer_type,
                        "content": qualifier_content
                         })
                for qualifier in qualifiers:
                    feature[qualifier['type']] = qualifier['content']
                features.append(feature)
            feature_type = None
            genomic_strand = '+'
            genomic_positions = None
            qualifer_type = None
            qualifier_content = None
            qualifiers = []
            feature_type = tokens[1].strip()
            location = tokens[2].strip()
        elif not qualifer_type and not qualifier_content and len(tokens) == 2 and re.findall('\.\.',tokens[1]): #still the content of the current location
            location += tokens[1].strip()
        elif not start_of_sequence and qualifer_type and qualifier_content : # still the content of the current qualifier
            qualifier_content += " "+line.strip().replace('"','')
        elif line.startswith('//'): #end of the genomic sequence
            start_of_sequence = False
        elif start_of_sequence:
            pieces_of_seq.append(''.join(re.split('\s+',line.strip())[:-1]).upper())

    dna = DNA( name = accession, sequence=''.join(pieces_of_seq))

    if organism:
        dna.organism = organism

    for feature in features:
        if feature.has_key('ncRNA_class'):
            if feature['genomicStrand'] == '+':
                feature['sequence'] = dna.sequence[feature['genomicPositions'][0]-1:feature['genomicPositions'][-1]]
            else:
                feature['sequence'] = DNA(name = dna.name, sequence = dna.sequence[feature['genomicPositions'][0]-1:feature['genomicPositions'][-1]]).get_complement()[::-1]

    return dna, DataFrame(features)

def parse_rnaml(rnaml_data, canonical_only = False):
    """
    Parse RNAML data. At now, this method handles only single molecular secondary structures.

    Parameters:
    ---------
     - rnaml_data: the RNAML data as a String
     - canonical_only (default: False): if True, the helices will be made exclusively with canonical base-pairs: AU c( ), GC c( ) or GU c( ).
    
    Returns:
    ------
    a list of SecondaryStructure objects (see pyrna.features)
    """
    secondary_structures = []
    import xml.etree.ElementTree as ET
    rnaml_tree = ET.fromstring(rnaml_data)

    for molecule in rnaml_tree.findall('molecule'):

        rna = RNA(name = molecule.get('id'), sequence = re.sub('\s+','', molecule.find('sequence').find('seq-data').text))
        
        secondary_structure = SecondaryStructure(rna)

        if not canonical_only:       

            for helix in molecule.find('structure').find('model').find('str-annotation').findall('helix'):
                secondary_structure.add_helix(helix.get('id'), int(helix.find('base-id-5p').find('base-id').find('position').text), int(helix.find('base-id-3p').find('base-id').find('position').text), int(helix.find('length').text));

            for single_strand in molecule.find('structure').find('model').find('str-annotation').findall('single-strand'):
                end5 = int(single_strand.find('segment').find('base-id-5p').find('base-id').find('position').text)
                end3 = int(single_strand.find('segment').find('base-id-3p').find('base-id').find('position').text)
                secondary_structure.add_single_strand(single_strand.find('segment').find('seg-name').text, end5, end3-end5+1);

            for base_pair in molecule.find('structure').find('model').find('str-annotation').findall('base-pair'):
                edge1 = '('
                edge2 = ')'
                if base_pair.find('edge-5p').text == 'H':
                    edge1 = '['
                elif base_pair.find('edge-5p').text == 'S':
                    edge1 = '{'
                elif base_pair.find('edge-5p').text == 's':
                    edge1 = '{'
                elif base_pair.find('edge-5p').text == '!':
                    edge1 = '!'

                if base_pair.find('edge-3p').text == 'H':
                    edge2 = ']'
                elif base_pair.find('edge-3p').text == 'S':
                    edge2 = '}'
                elif base_pair.find('edge-3p').text == 's':
                    edge2 = '}'
                elif base_pair.find('edge-3p').text == '!':
                    edge2 = '!'

                secondary_structure.add_base_pair(base_pair.find('bond-orientation').text.upper(), edge1, edge2, int(base_pair.find('base-id-5p').find('base-id').find('position').text), int(base_pair.find('base-id-3p').find('base-id').find('position').text));
        
        else:
            canonical_bps = []
            non_canonical_bps = []
            for base_pair in molecule.find('structure').find('model').find('str-annotation').findall('base-pair'):
                edge1 = '('
                edge2 = ')'
                if base_pair.find('edge-5p').text == 'H':
                    edge1 = '['
                elif base_pair.find('edge-5p').text == 'S':
                    edge1 = '{'
                elif base_pair.find('edge-5p').text == 's':
                    edge1 = '{'
                elif base_pair.find('edge-5p').text == '!':
                    edge1 = '!'

                if base_pair.find('edge-3p').text == 'H':
                    edge2 = ']'
                elif base_pair.find('edge-3p').text == 'S':
                    edge2 = '}'
                elif base_pair.find('edge-3p').text == 's':
                    edge2 = '}'
                elif base_pair.find('edge-3p').text == '!':
                    edge2 = '!'

                orientation = base_pair.find('bond-orientation').text.upper()

                pos1 = int(base_pair.find('base-id-5p').find('base-id').find('position').text)
                residue1 = secondary_structure.rna.sequence[pos1-1]
                pos2 = int(base_pair.find('base-id-3p').find('base-id').find('position').text)
                residue2 = secondary_structure.rna.sequence[pos2-1]

                canonical_bps.append([orientation, edge1, edge2, pos1, pos2]) if utils.is_canonical(residue1, residue2, orientation, edge1, edge2) else non_canonical_bps.append([orientation, edge1, edge2, pos1, pos2])
             
            secondary_structure = base_pairs_to_secondary_structure(secondary_structure.rna, DataFrame(canonical_bps, columns=['orientation', 'edge1', 'edge2', 'pos1', 'pos2']))

            for bp in non_canonical_bps: #the non-canonical interactions are tertiary ones
                secondary_structure.add_tertiary_interaction(bp[0], bp[1], bp[2], bp[3], bp[4])

        secondary_structures.append(secondary_structure)

    return secondary_structures


def parse_fasta(fasta_data, type='RNA'):
    """
    Parse FASTA data

    Parameters:
    ---------
    - fasta_data: the Fasta data as a String
    - type (default: 'RNA'): can be equal to 'DNA' or 'RNA'

    Returns:
    ------
    a list of RNA or DNA objects (according to the value of the parameter type) (see pyrna.features)
    """
    molecules = []
    pieces = []
    molecule_name = None
    for line in fasta_data.split('\n'):
        if re.match('>',line):
            if molecule_name and len(pieces) > 0:
                if type == 'RNA':
                    m = RNA(sequence = ''.join(pieces), name = molecule_name.strip())
                elif type == 'DNA':
                    m = DNA(sequence = ''.join(pieces), name = molecule_name.strip())
                if m != None:
                    molecules.append(m)
            molecule_name = line[1:]
            pieces = []
        else:
            pieces.append(line.strip().upper())
    #last molecule
    if molecule_name and len(pieces) > 0:
        if type == 'RNA':
            m = RNA(sequence = ''.join(pieces), name = molecule_name.strip())
        elif type == 'DNA':
            m = DNA(sequence = ''.join(pieces), name = molecule_name.strip())
        if m != None:
            molecules.append(m)
    return molecules

def parse_vienna(vienna_data):
    """
    Parse Vienna data

    Parameters:
    ---------
     - vienna_data: the Vienna data as a String

    Returns:
    ------
    list of RNA objects, list of pandas Dataframes (each Dataframe listing base pairs). If more than 1 RNA objects are returned along with a single Dataframe, this Dataframe describes a consensus secondary structure.
    """
    name = None
    secondary_structures = []
    rnas = []
    current_bn = []
    current_sequence = []
    for line in vienna_data.split('\n'):
        if re.match('^[\.()]+$', line):
            current_bn.append(line)
        elif re.match('^>', line):
            if len(current_bn) and len(current_sequence):
                rnas.append(RNA(name = name, sequence = ''.join(current_sequence)))
                secondary_structures.append(parse_bn(''.join(current_bn)))
            name = line[1:]
            current_bn = []
            current_sequence = []
        elif len(line.strip()):
            current_sequence.append(line.strip())

    #last one
    if len(current_bn) and len(current_sequence):
        rnas.append(RNA(name = name, sequence = ''.join(current_sequence)))
        secondary_structures.append(parse_bn(''.join(current_bn)))
    return rnas, secondary_structures

def parse_bn(bn):
    """
    Parse a bracket notation.
    
    Parameters:
    ---------
     - bn: the bracket notation as a String

    Returns:
    ------
    a pandas Dataframe listing the base pairs. Returns an empty Dataframe if no base-pairs are found.
    """
    i = 0
    lastPairedPos = []
    lastPairedSymbol = []
    basePairs = []
    for pos in list(bn):
        i+=1
        if pos == '(' or pos == '{' or pos == '[':
            lastPairedPos.append(i)
            lastPairedSymbol.append(pos)
        elif pos == ')' or pos == '}' or pos == ']':
            basePairs.append(['c', lastPairedSymbol.pop(), pos, lastPairedPos.pop(), i])

    if len(basePairs):
        return DataFrame(basePairs, columns=['orientation', 'edge1', 'edge2', 'pos1', 'pos2'])
    else:
        return DataFrame()

def parse_clustalw(clustalw_data):
    """
    Parse Clustalw data

    Parameters:
    ---------
     - clustalw_data: the Clustalw data as a String

    Returns:
    ------
    a tuple containing: 
    - a list of gapped or ungapped RNA objects
    - a pandas Dataframe listing the paired positions of consensus secondary structure)
    """

    bn = None
    alignedSequences = {}
    lines = clustalw_data.strip().split('\n')
    for line in lines:
        line = line.strip()
        if line.startswith('2D'):
            bn = line.split('\t')[1]
        elif len(line) and not line.startswith('#'):
            tokens = line.split('\t')
            alignedSequences[tokens[0]] = alignedSequences.get(tokens[0], "") + tokens[1] 

    rnas = []
    
    for key in alignedSequences.keys():
        rna = RNA(name=key, sequence=alignedSequences[key])
        rnas.append(rna)

    return rnas, parse_bn(bn)   

def parse_stockholm(stockholm_data):
    """
    Parse Stokholm data

    Parameters:
    ---------
     - stockholm_data: the Stockholm data as a String

    Returns:
    ------
    a tuple containing: 
    - a list of gapped or ungapped RNA objects
    - a dict of organism names (keys)  and accession numbers/start-end (values)
    -a pandas Dataframe listing the paired positions of consensus secondary structure)
    """
    alignedSequences = {}
    organisms={}
    aligned2D = ""
    rfam_id = None
    lines = stockholm_data.strip().split('\n')
    for line in lines:
        tokens = re.split('\s+', line)
        if len(line) != 0 and not re.match('^#', line) and len(tokens) == 2:
            if alignedSequences.has_key(tokens[0]):
                alignedSequences[tokens[0]] = alignedSequences[tokens[0]]+tokens[1]
            else:
                alignedSequences[tokens[0]] = tokens[1]
        elif len(line) != 0 and re.match('^#=GC SS_cons', line):
            aligned2D += re.sub('>', ')', re.sub('<', '(', tokens[2]))
        elif len(line) != 0 and re.match('^#=GF AC', line):
            rfam_id = tokens[2].strip()
        elif len(line) != 0 and re.match('^#=GS', line):
            organisms[tokens[1]] = tokens[-1]

    rnas = []

    for key in alignedSequences.keys():
        rna = RNA(name=key, sequence=alignedSequences[key])
        rna.source = 'db:rfam:'+rfam_id
        if not key.split('/') == 2:
            rna.organism = key
        rnas.append(rna)

    return (rnas, organisms, parse_bn(aligned2D))

def parse_pdb(pdb_data):
    """
    Parse PDB data.

    Parameters:
    ---------
     - pdb_data: the PDB data as a String

    Returns:
    ------
    a list of TertiaryStructure objects (see pyrna.features). if the PDB data describes a tertiary structure made with several molecular chains, this method will return one TertiaryStructure object per chain.
    """
    rnas = []
    tertiary_structures = []
    current_chain = None
    current_residue = None
    current_residue_pos = None
    absolute_position = -1
    current_rna = None
    current_3D = None
    title = "N.A."

    for line in pdb_data.split('\n'):
        header = line[0:6].strip()
        atom_name = line[12:16].strip()
        residue_name = line[17:20].strip()
        chain_name = line[21:22].strip()
        residue_pos = line[22:27].strip()

        if (header == "ATOM" or header == "HETATM") and not residue_name in ["FMN","PRF","HOH","MG","OHX","MN","ZN", "SO4", "CA", "UNK"] and not atom_name in ["MG","K", "NA", "SR", "CL", "CD", "ACA"] and len(chain_name):
            if chain_name != current_chain: #new chain
                current_residue = residue_name
                current_residue_pos = residue_pos
                current_chain = chain_name
                absolute_position = 1
                current_rna = RNA(sequence="", name=current_chain)
                current_rna.add_residue(current_residue)
                current_3D = TertiaryStructure(current_rna)
                current_3D.title = re.sub(' +', ' ', title)
                current_3D.numbering_system[str(absolute_position)] = current_residue_pos

            elif current_residue_pos != residue_pos: # new residue
                current_residue = residue_name
                current_residue_pos = residue_pos
                current_rna.add_residue(current_residue)
                absolute_position += 1
                current_3D.numbering_system[str(absolute_position)] = current_residue_pos

            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            current_3D.add_atom(atom_name, absolute_position, [x,y,z])

            if (atom_name == "O4'" or atom_name == "O4*") and not current_rna in rnas:
                rnas.append(current_rna)
                tertiary_structures.append(current_3D)

        elif header == 'TITLE':
            title += line[10:]

    return tertiary_structures

def parse_sam(sam_file):
    """
    This method parses a SAM file. It delegates the low-level parsing to the pysam library (https://code.google.com/p/pysam/). At now, this method is not able to handle oriented and paired-ends reads.

    Parameters:
    ---------
     - sam_file: the absolute path of the SAM file as a String

    Returns:
    ------
    a tuple containing:
    - a list of aligned reads (each read is described as a dict like: {'tid':int, 'genomicStart':int, 'genomicEnd':int} )
    - the total number of reads described in the SAM file
    - a dictionary providing the correspondance between the name of the genomic sequences and the tids available in the SAM file
    """
    reads = []
    tid_dic = {}
    total_read_nb = 0
    aligned_read_nb = 0
    sam_file_content = Samfile(sam_file, 'r') #samfile: <type 'csamtools.Samfile'> ; 'r' for sam files and 'rb' for bam files

    sequence_nb = len(sam_file_content.header['SQ'])
    sequence_names = []

    for i in range(0, sequence_nb):
        reads.append([])
        sequence_names.append(sam_file_content.header['SQ'][i]['SN'])
    
    for alignedread in sam_file_content.fetch(): #alignedread: <type 'csamtools.AlignedRead'>
        total_read_nb += 1
        if alignedread.rname != -1:
            aligned_read_nb +=1
            for name in sequence_names:
                if name == sam_file_content.getrname(alignedread.rname):
                    if tid_dic.has_key(alignedread.rname) == 0:
                        tid_dic[alignedread.rname] = name

            c = 0 #counter for the read sequence which is incremented with each M/I/S
            r = alignedread.pos #counter for the reference sequence which is incremented with each M/D/S ; POS=Integer=position(-1) of the first base of the reference sequence that gives a match (M)
            tag = 'no'
            read = {'tid':alignedread.rname}
            #read = [alignedread.rname]
            aligned_seq = []
            for tu in alignedread.cigar:
                if tu[0] == 0: #0=M=alignment Match (match or mismatch)
                    tag = 'yes'
                    read_pos_list = range(c,c+tu[1],1)
                    for p in read_pos_list: #p+1=position on the read sequence
                        r += 1 #position on the reference sequence
                        aligned_seq.append(r)
                    c += tu[1]
                elif tu[0] == 1: #1=I=Insertion
                    c += tu[1]
                elif tu[0] == 2: #2=D=Deletion
                    r += tu[1]
                elif tu[0] == 4: #4=S=Soft clipping (match or mismatch or no alignment)
                    if tag == 'no':
                        r -= tu[1]
                    c += tu[1]
                    r += tu[1]
            read['genomicStart'] = aligned_seq[0] 
            #read.append(aligned_seq[0])
            read['genomicEnd'] = aligned_seq[-1]
            #read.append(aligned_seq[-1])
            reads[alignedread.rname-1].append(read)
    sam_file_content.close()
    return reads, total_read_nb, tid_dic
