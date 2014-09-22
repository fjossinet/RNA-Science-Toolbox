#!/usr/bin/env python

"""
Script to import and manage charnDB data in Mongodb
===================================================

Step 1: each RFAM family is recovered (seed alignment) and imported in the database
Step 1: fungal genomic sequences and annotations are downloaded from the NCBI FTP 
Step 3: fungal genomic sequences and annotations are imported in the database

Basic usage:
-----------

charnDB_mongo.py: creates the database (default behavior) and stores the data in a Mongodb instance localhost:27017 (options -p and -h are available to change the database parameters)

Other options available with this script:
-----------------------------------------
charnDB_mongo.py -add genbank_or_embl_files: add a new species to the database
charnDB_mongo.py -remove species_name: removes everything related to the species given as argument
charnDB_mongo.py -drop: drop the database
charnDB_mongo.py -dump: dump data into files (JSON files for the charnDB website and GFF3 files to be used in genome browsers).
"""

import sys, re, os, json, commands
from pyrna.features import DNA
from pyrna.db import NCBI
from pyrna.db import Rfam
from pyrna import parsers
from pyrna.parsers import parse_genbank, to_clustalw, parse_stockholm, to_fasta, parse_clustalw
from pyrna.utils import get_file_as_source
from pymongo import MongoClient
from bson.objectid import ObjectId
from pandas import isnull, DataFrame
import numpy as np

restrict_to_fungal_species = False

def create(base_path, db_host = "localhost", db_port = 27017):

    client = MongoClient(db_host, db_port)
    databases_names  = client.database_names()

    #STEP1: fetch and store the RFAM data (seed alignments and covariance models)

    rfam = Rfam()
    rfam.generate_CMs()
    rfam.generate_seed_alignments()
    familiesDetails = rfam.get_families_details()
    rfam_ids_to_process = []
    for index, family in familiesDetails.iterrows():
        rfam_ids_to_process.append(int(family['accession'].split('RF')[1]))
    rfam_ids_to_process.sort()

    for id in rfam_ids_to_process:
        id = "RF%05u"%id
        if len(familiesDetails[familiesDetails['accession'] == id]) and not client['comparative_genomics']['alignments'].find_one({'source': "db:rfam:%s" % id}): #some RFAM families are not available
            stockholm_content = rfam.get_entry(rfam_id = id, format='stockholm')
            (aligned_molecules, organisms, consensus2D) = parse_stockholm(stockholm_content)
            in_house_alignment = {
                '_id': str(ObjectId()),
                'source': "db:rfam:%s" % id,
                'name': familiesDetails[familiesDetails['accession'] == id].iloc[0,4],
                'class': familiesDetails[familiesDetails['accession'] == id].iloc[0,2],
                'alignment': to_clustalw(consensus2D, aligned_molecules),
                'cm': rfam.get_CM(id)
            }

            client['comparative_genomics']['alignments'].insert(in_house_alignment)
    
    #STEP2: fetch the genomic data from NCBI (through Eutils)

    if not os.path.exists(base_path):
        os.makedirs(base_path)

    ncbi = NCBI()

    ncbi_taxonomic_accession_numbers = []

    species2accession_numbers = get_fungal_genomic_entries_from_ncbi(ncbi, base_path) #we recover the acc numbers for the fungal entries recovered from the NCBI FTP

    ncbi_path = base_path+"/genomes/NCBI/"

    if not os.path.exists(ncbi_path):
        os.makedirs(ncbi_path)

    for species in species2accession_numbers:
        remaining_genomic_entries_to_dl = []
        species_path = ncbi_path+species+"/"
        if not os.path.exists(species_path):
            os.makedirs(species_path) 
        for accession_number in species2accession_numbers[species]:
            genomic_entry_path = species_path+accession_number+".gb"       
            if not os.path.exists(genomic_entry_path):
                remaining_genomic_entries_to_dl.append(accession_number)
        print "Remaining genome entries to dl for %s:"%species, len(remaining_genomic_entries_to_dl)

        downloaded = 0
        from pyrna import utils
        chunks = utils.chunks(remaining_genomic_entries_to_dl, 200) #with efetch, if more than about 200 UIDs are to be provided, the request should be made using the HTTP POST method
        for chunk in chunks:
            print len(chunk),  ' genomic entries to sent to efetch:', ','.join(chunk)
            genbank_content = ncbi.efetch(ids = [','.join(chunk)], db='nucleotide', rettype='gbwithparts')
            import re
            for file_content in re.split('\n//\n\n', genbank_content):
                match = re.findall("\nVERSION\s+(.+?)\s+",file_content)
                if len(match):
                    genomic_entry_path = species_path+match[0]+".gb"
                    with open(genomic_entry_path, 'w') as h:
                        file_content.strip()
                        file_content += "\n//"
                        h.write(file_content)

            downloaded += len(chunk)
            print "Downloaded %i/%i genomic entries for species %s"%(downloaded, len(remaining_genomic_entries_to_dl), species) 

    #STEP3: parse and store the genomic sequences and annotations for each genomic file in a species directory    
    for species in os.listdir(ncbi_path):
        i = 0
        genomic_files = os.listdir(ncbi_path+"/"+species)

        print "Parsing of %i genomic entries for species %s"%(len(genomic_files),species)

        for genomic_file in genomic_files:
            accession_number = genomic_file.split('.gb')[0]
            genomic_entry_path = ncbi_path+"/"+species+"/"+genomic_file
            i += 1
            print "Processing genomic entry %s (%s/%s)"%(genomic_file, i, len(genomic_files))

            genomic_sequence = None

            with open(genomic_entry_path) as h:
                genbank_content = h.read()

            try:
                genomic_sequence, features = parse_genbank(genbank_content)
            except Exception, e:
                print e
                os.remove(genomic_entry_path)
            
            db = client[genomic_sequence.organism.replace(' ', '_').replace('.', '_')]

            if db['genomes'].find_one({'source':"db:ncbi:%s"%accession_number}):
                print "%s already stored!!"%genomic_sequence.name
            else:
                #first we store the genome
                genome_description = {
                    "name": genomic_sequence.name,
                    "_id": genomic_sequence._id,
                    "sequence": genomic_sequence.sequence,
                    "source":'db:ncbi:%s'%accession_number,
                    "organism": genomic_sequence.organism,
                    'lineage': genomic_sequence.lineage
                }            
                db['genomes'].insert(genome_description)

                #then its annotations
                for (index, row) in features.iterrows():
                    annotation = {
                        '_id': str(ObjectId()),
                        'class': row['type'],
                        'source':'db:ncbi:%s'%accession_number,
                        'organism': genomic_sequence.organism,
                        'genome': "%s@genomes"%genomic_sequence._id,
                        'genomeName': genomic_sequence.name,
                        'genomicStrand': row['genomicStrand'],
                        'genomicPositions': row['genomicPositions'],
                        'lineage': genomic_sequence.lineage
                    }

                    for key in row.keys():
                        if not key in ['type', 'genomicPositions', 'genomicStrand'] and not isnull(row[key]):
                            annotation[key] = row[key]

                    #we extract the sequence
                    if row['genomicStrand'] == '+':
                        annotation['sequence'] = genomic_sequence.sequence[row['genomicPositions'][0]-1:row['genomicPositions'][1]]
                    else:
                        annotation['sequence'] = DNA(name = genomic_sequence.name, sequence = genomic_sequence.sequence[row['genomicPositions'][0]-1:row['genomicPositions'][1]]).get_complement()[::-1]

                    db['annotations'].insert(annotation)
    
    client.disconnect()

def to_newick_format(db_host = "localhost", db_port = 27017):
    """
    This function produces a phylogenetic tree for all the species available in mongodb. This tree is described using the Newick format (http://en.wikipedia.org/wiki/Newick_format). 
    """
    client = MongoClient(db_host, db_port)
    databases_names  = client.database_names()
    taxonomy = {}
    taxons = None

    print "Total species: ", len(databases_names)
    
    for species in databases_names:
        db = client[species] 

        genome  = db['genomes'].find_one()

        if genome: #some databases can be about something else than a species
            if genome.get('lineage', None):
                taxons = genome['lineage'][:-1].split('; ')

                for i in range(len(taxons)):
                    taxon = taxons[i]
                    if taxonomy.has_key(taxon) and i < len(taxons)-1 and not taxons[i+1] in taxonomy[taxon]:
                        taxonomy[taxon].append(taxons[i+1])
                    elif not taxonomy.has_key(taxon) and i < len(taxons)-1:
                        taxonomy[taxon] = [taxons[i+1]]

                if taxonomy.has_key(taxons[-1]) and not genome['organism'] in taxonomy[taxons[-1]]:
                    taxonomy[taxons[-1]].append(genome['organism'])
                elif not taxonomy.has_key(taxons[-1]):
                    taxonomy[taxons[-1]] = [genome['organism']]
            else:
                print "No lineage for", species
    import re
    newick = walk(taxonomy, taxons[0], s="", current_depth = 1, max_depth = 1)
    print '('+re.sub(r'\)\(', r'),(', newick)+')'

    client.close()

def walk(taxonomy, taxon, s, current_depth, max_depth): 
    if taxonomy.has_key(taxon):
        children = taxonomy[taxon]
        if not children:
            s += taxon
        else:
            for child in children:
                s += '('+taxon
                current_depth += 1
                s = walk(taxonomy, child, s, current_depth, max_depth)
                s += ')'
                current_depth -= 1
        return s  
    else:
        max_depth = current_depth
        s += '('+taxon+')'
        return s    

def dump(db_host = "localhost", db_port = 27017):
    """
    This function dumps data for the ncRNA website
    """
    client = MongoClient(db_host, db_port)
    databases_names  = client.database_names()
    species = []

    for database_name in databases_names:
        if not database_name in ["PDB", "RNA3DHub", "Candida_hispaniensis_reannotated", "logs", "comparative_genomics", "local", "admin"]:
            db = client[database_name]
            species.append(database_name)
    
    #### all the ncRNA and annotations as GFF3 files (to be used with a Genome browser for example)
    directory = "genomes"
    if not os.path.isdir(directory):
        print "Processing genomes and annotations..."
        os.makedirs(directory)
        for database_name in databases_names:
            if not database_name in ["PDB", "RNA3DHub", "Candida_hispaniensis_reannotated", "logs", "comparative_genomics", "local", "admin"]:
                db = client[database_name]
                molecules = []
                for genome in db['genomes'].find():
                    molecules.append(DNA(name = genome['name'], sequence = genome['sequence']))    
                with open(os.path.join(directory, '%s.fasta'%database_name), 'w') as h:
                    h.write(to_fasta(molecules))
                id = 0
                with open(os.path.join(directory, '%s.gff3'%database_name), 'w') as h:
                    for annotation in db['annotations'].find():
                        id += 1 
                        h.write("%s\t%s\t%s\t%i\t%i\t%s\t%s\t.\tID=%s\n"%(annotation['genomeName'], annotation.get('source', '?'), annotation['class'], annotation['genomicPositions'][0], annotation['genomicPositions'][-1], str(annotation.get('score', '.')),  annotation['genomicStrand'],  annotation.get('product', annotation.get('locus_tag', 'unknown'))))
                    for ncrna in db['ncRNAs'].find():
                        id += 1 
                        h.write("%s\t%s\t%s\t%i\t%i\t%s\t%s\t.\tID=%s\n"%(ncrna['genomeName'], ncrna.get('source', '?'), ncrna['class'], ncrna['genomicPositions'][0], ncrna['genomicPositions'][-1], str(ncrna.get('score', '.')),  ncrna['genomicStrand'], ncrna['name']))

    ##### the data for the webpage species.html
    if not os.path.exists('species.csv'):
        print "Processing species.csv..."
        data = {}
        for s in species:
            db = client[s]
            lineage = db['genomes'].find_one()['lineage'].split(';')

            if lineage[2] == " Microsporidia":
                continue
            
            levels1 = data.get('Level1', [])
            levels1.append(lineage[3].strip())
            data['Level1'] = levels1

            levels2 = data.get('Level2', [])
            levels2.append(lineage[4].strip())
            data['Level2'] = levels2

            levels3 = data.get('Level3', [])
            levels3.append(lineage[5].strip())
            data['Level3'] = levels3

            levels4 = data.get('Level4', [])
            levels4.append(s)
            data['Level4'] = levels4

            cdbox_counts = 0
            haca_counts = 0
            srna_counts = 0
            snrna_counts = 0
            riboswitch_counts = 0
            ribozyme_counts = 0
            gene_counts = 0
            cisreg_counts = 0
            all_ncRNAs_counts = 0
            
            for alignment in client["comparative_genomics"]['alignments'].find():  
                aligned_rnas, consensus_2d = parse_clustalw(alignment['alignment'])
                for aligned_rna in aligned_rnas:
                    tokens = aligned_rna.name.split('@')
                    if len(tokens) == 2 and tokens[1] == s:
                        if re.match('^.*CD-box.*$', alignment['class']):
                            cdbox_counts += 1
                        elif re.match('^.*HACA-box.*$', alignment['class']):
                            haca_counts += 1
                        elif re.match('^.*sRNA.*$', alignment['class']):
                            srna_counts += 1
                        elif re.match('^.*snRNA.*$', alignment['class']):
                            snrna_counts += 1
                        elif re.match('^.*riboswitch.*$', alignment['class']):
                            riboswitch_counts += 1
                        elif re.match('^.*ribozyme.*$', alignment['class']):
                            ribozyme_counts += 1
                        elif re.match('^.*Gene$', alignment['class']):
                            gene_counts += 1
                        elif re.match('^.*Cis-reg$', alignment['class']):
                            cisreg_counts += 1
                        all_ncRNAs_counts += 1
            
            all_ncrnas = data.get('all_ncrnas', [])
            all_ncrnas.append(all_ncRNAs_counts)
            data['all_ncrnas'] = all_ncrnas

            riboswitches = data.get('riboswitches', [])
            riboswitches.append(riboswitch_counts)
            data['riboswitches'] = riboswitches

            cdsnos = data.get('cdsnos', [])
            cdsnos.append(cdbox_counts)
            data['cdsnos'] = cdsnos

            hacasnos = data.get('hacasnos', [])
            hacasnos.append(haca_counts)
            data['hacasnos'] = hacasnos

            ribozyme = data.get('ribozymes', [])
            ribozyme.append(ribozyme_counts)
            data['ribozymes'] = ribozyme

            gene = data.get('genes', [])
            gene.append(gene_counts)
            data['genes'] = gene

            cisreg = data.get('cisregs', [])
            cisreg.append(cisreg_counts)
            data['cisregs'] = cisreg

            srna = data.get('srnas', [])
            srna.append(srna_counts)
            data['srnas'] = srna

            snrna = data.get('snrnas', [])
            snrna.append(snrna_counts)
            data['snrnas'] = snrna
        
        df = DataFrame.from_dict(data)
        df.to_csv('species.csv', index_label = "Species")

    ##### data for the webpage families.html
    directory = "families_counts"

    if not os.path.isdir(directory):
        print "Processing families counts..."
        os.makedirs(directory)

        db = client["comparative_genomics"]
        cdbox_all_species_counts = {}
        haca_all_species_counts = {}
        srna_all_species_counts = {}
        snrna_all_species_counts = {}
        riboswitch_all_species_counts = {}
        ribozyme_all_species_counts = {}
        gene_all_species_counts = {}
        cisreg_all_species_counts = {}
        all_ncRNAs_all_species_counts = {}
        all_predictions_sizes = {}
        i = 0
        for alignment in db['alignments'].find(timeout = False):
            i += 1
            aligned_rnas, consensus_2d = parse_clustalw(alignment['alignment'])
            cdbox_species_counts = []
            haca_species_counts = []
            srna_species_counts = []
            snrna_species_counts = []
            riboswitch_species_counts = []
            ribozyme_species_counts = []
            gene_species_counts = []
            cisreg_species_counts = []
            all_ncRNAs_species_counts = []
            sizes = []
            for s in species:
                count = 0
                for aligned_rna in aligned_rnas:
                    tokens = aligned_rna.name.split('@')
                    if len(tokens) == 2 and tokens[1] == s:
                        sizes.append(len(aligned_rna.sequence.replace('-','')))
                        count += 1
                if re.match('^.*CD-box.*$', alignment['class']):
                    cdbox_species_counts.append(count)
                elif re.match('^.*HACA-box.*$', alignment['class']):
                    haca_species_counts.append(count)
                elif re.match('^.*sRNA.*$', alignment['class']):
                    srna_species_counts.append(count)
                elif re.match('^.*snRNA.*$', alignment['class']):
                    snrna_species_counts.append(count)
                elif re.match('^.*riboswitch.*$', alignment['class']):
                    riboswitch_species_counts.append(count)
                elif re.match('^.*ribozyme.*$', alignment['class']):
                    ribozyme_species_counts.append(count)
                elif re.match('^.*Gene$', alignment['class']):
                    gene_species_counts.append(count)
                elif re.match('^.*Cis-reg$', alignment['class']):
                    cisreg_species_counts.append(count)
                all_ncRNAs_species_counts.append(count)
            if re.match('^.*CD-box.*$', alignment['class']):
                cdbox_species_counts += [np.std(sizes)]
                cdbox_all_species_counts[alignment['name']] = cdbox_species_counts
            elif re.match('^.*HACA-box.*$', alignment['class']):
                haca_species_counts += [np.std(sizes)]
                haca_all_species_counts[alignment['name']] = haca_species_counts
            elif re.match('^.*sRNA.*$', alignment['class']):
                srna_species_counts += [np.std(sizes)]
                srna_all_species_counts[alignment['name']] = srna_species_counts
            elif re.match('^.*snRNA.*$', alignment['class']):
                snrna_species_counts += [np.std(sizes)]
                snrna_all_species_counts[alignment['name']] = snrna_species_counts
            elif re.match('^.*riboswitch.*$', alignment['class']):
                riboswitch_species_counts += [np.std(sizes)]
                riboswitch_all_species_counts[alignment['name']] = riboswitch_species_counts
            elif re.match('^.*ribozyme.*$', alignment['class']):
                ribozyme_species_counts += [np.std(sizes)]
                ribozyme_all_species_counts[alignment['name']] = ribozyme_species_counts
            elif re.match('^.*Gene$', alignment['class']):
                gene_species_counts += [np.std(sizes)]
                gene_all_species_counts[alignment['name']] = gene_species_counts
            elif re.match('^.*Cis-reg$', alignment['class']):
                cisreg_species_counts += [np.std(sizes)]
                cisreg_all_species_counts[alignment['name']] = cisreg_species_counts
            all_ncRNAs_species_counts += [np.std(sizes)]
            all_ncRNAs_all_species_counts[alignment['name']] = all_ncRNAs_species_counts
        
        column_names = species + ["standard_deviation"]

        df = DataFrame.from_dict(cdbox_all_species_counts, orient='index')
        df.columns = column_names
        df = df.loc[df.sum(axis=1) != 0] #we remove all the rows (the rfam families) with no hit in the species
        df = df.sort(axis = 0)
        df.to_csv(os.path.join(directory,'cdbox_counts.csv'))

        df = DataFrame.from_dict(haca_all_species_counts, orient='index')
        df.columns = column_names
        df = df.loc[df.sum(axis=1) != 0] #we remove all the rows (the rfam families) with no hit in the species
        df = df.sort(axis = 0)
        df.to_csv(os.path.join(directory,'haca_counts.csv'))

        df = DataFrame.from_dict(srna_all_species_counts, orient='index')
        df.columns = column_names
        df = df.loc[df.sum(axis=1) != 0] #we remove all the rows (the rfam families) with no hit in the species
        df = df.sort(axis = 0)
        df.to_csv(os.path.join(directory,'srna_counts.csv'))

        df = DataFrame.from_dict(snrna_all_species_counts, orient='index')
        df.columns = column_names
        df = df.loc[df.sum(axis=1) != 0] #we remove all the rows (the rfam families) with no hit in the species
        df = df.sort(axis = 0)
        df.to_csv(os.path.join(directory,'snrna_counts.csv'))

        df = DataFrame.from_dict(riboswitch_all_species_counts, orient='index')
        df.columns = column_names
        df = df.loc[df.sum(axis=1) != 0] #we remove all the rows (the rfam families) with no hit in the species
        df = df.sort(axis = 0)
        df.to_csv(os.path.join(directory,'riboswitch_counts.csv'))

        df = DataFrame.from_dict(ribozyme_all_species_counts, orient='index')
        df.columns = column_names
        df = df.loc[df.sum(axis=1) != 0] #we remove all the rows (the rfam families) with no hit in the species
        df = df.sort(axis = 0)
        df.to_csv(os.path.join(directory,'ribozyme_counts.csv'))

        df = DataFrame.from_dict(gene_all_species_counts, orient='index')
        df.columns = column_names
        df = df.loc[df.sum(axis=1) != 0] #we remove all the rows (the rfam families) with no hit in the species
        df = df.sort(axis = 0)
        df.to_csv(os.path.join(directory,'gene_counts.csv'))

        df = DataFrame.from_dict(cisreg_all_species_counts, orient='index')
        df.columns = column_names
        df = df.loc[df.sum(axis=1) != 0] #we remove all the rows (the rfam families) with no hit in the species
        df = df.sort(axis = 0)
        df.to_csv(os.path.join(directory,'cisreg_counts.csv'))

        df = DataFrame.from_dict(all_ncRNAs_all_species_counts, orient='index')
        df.columns = column_names
        df = df.loc[df.sum(axis=1) != 0] #we remove all the rows (the rfam families) with no hit in the species
        df = df.sort(axis = 0)
        df.to_csv(os.path.join(directory,'all_ncrnas_counts.csv'))
    
    ##### data for the webpage synthenies.html
    directory = "synthenies"

    if not os.path.isdir(directory):
        print "Processing synthenies..."
        os.makedirs(directory)

        db = client["comparative_genomics"]

        for alignment in db['alignments'].find(timeout = False):
            aligned_rnas, consensus_2d = parse_clustalw(alignment['alignment'])
            orthologs = []
            for s in species:
                db = client[s]
                for aligned_rna in aligned_rnas:
                    tokens = aligned_rna.name.split('@')
                    if len(tokens) == 2 and tokens[1] == s:
                        ncRNA = db['ncRNAs'].find_one({'_id':tokens[0].split('[New]')[-1]})
                        close_annotations = []
                        upstream = ncRNA['genomicPositions'][0] - 2000
                        downstream = ncRNA['genomicPositions'][1] + 2000

                        for annotation in db['annotations'].find({'genome':ncRNA['genome']}):
                            if annotation['genomicPositions'][0] >= upstream and annotation['genomicPositions'][1] <= downstream :
                                annotation.pop("lineage", None)
                                annotation.pop("sequence", None)
                                annotation.pop("translation", None)
                                annotation.pop("genome", None)
                                annotation.pop("_id", None)

                                close_annotations.append(annotation)

                        print "%i close annotations in %s for %s"%(len(close_annotations), ncRNA['organism'], alignment['name'])

                        ncRNA.pop("sequence", None)
                        ncRNA.pop("genome", None)
                        ncRNA.pop("alignment", None)
                        ncRNA.pop("_id", None)

                        orthologs.append({
                            'ortholog': ncRNA,
                            'close_annotations': close_annotations,
                            'genomic_range': [upstream, downstream]
                        })
            orthologs = sorted(orthologs, key=lambda ortholog: ortholog['ortholog']['organism'])
            with open(os.path.join(directory,'%s.json'%alignment['name']), 'w') as h:
                h.write(json.dumps(orthologs, indent=4))

    ##### data for the webpage consensus_structure.html
    directory = "alignments"

    if not os.path.isdir(directory):
        print "Processing alignments..."
        os.makedirs(directory)

        db = client["comparative_genomics"]

        from pyrna.db import NCBI
        
        ncbi = NCBI()
        for alignment in db['alignments'].find(timeout = False):
            aligned_rnas, consensus_2d = parse_clustalw(alignment['alignment'])
            bona_fide_sequences = []
            fungal_sequences = []
            for aligned_rna in aligned_rnas:
                if aligned_rna.name.find('@') != -1: #its a fungal sequence. We need to fits its name to the bone fide sequences (imported from RFAM) => species_name/accession_number/pos1-pos2
                    tokens = aligned_rna.name.split('@')
                    _id = tokens[0].split('[New]')[-1]
                    species = tokens[-1]
                    ncRNA = client[species]['ncRNAs'].find_one({'_id':_id})
                    name = None
                    if ncRNA['genomicStrand'] == '+':
                        name =  "%s/%s/%i-%i"%(tokens[-1], ncRNA['genomeName'], ncRNA['genomicPositions'][0], ncRNA['genomicPositions'][-1])
                    else:
                        name =  "%s/%s/%i-%i"%(tokens[-1], ncRNA['genomeName'], ncRNA['genomicPositions'][-1], ncRNA['genomicPositions'][0])
                    aligned_rna.name = name
                    fungal_sequences.append(aligned_rna) 
                else:
                    bona_fide_sequences.append(aligned_rna)   
            #now we dump the alignment
            if len(bona_fide_sequences):
                with open(os.path.join(directory,'%s_%s.aln'%(alignment['name'], 'bona_fide')), 'w') as h:
                    h.write(to_clustalw(consensus_2d, bona_fide_sequences, curate = True))
            if len(fungal_sequences):
                with open(os.path.join(directory,'%s_%s.aln'%(alignment['name'], 'fungi')), 'w') as h:
                    h.write(to_clustalw(consensus_2d, fungal_sequences, curate = True))

def drop(db_host = "localhost", db_port = 27017):
    client = MongoClient(db_host, db_port)

    databases_names  = client.database_names()

    for database_name in databases_names:
        if not database_name in ["PDB", "RNA3DHub", "Candida_hispaniensis_reannotated", "logs"]:
            client.drop_database(database_name)

    client.close()

def add_species(species_name, files, db_host = "localhost", db_port = 27017):
    """
    This method adds a new species to the database. it accept EMBL or Genbank files. This species will be annotated with the next iteration of annotation.
    """
    client = MongoClient(db_host, db_port)
    db = client[species_name]

    genomic_sequences = []
    for file in files:
        print "Processing %s..."%file

        with open(file) as f:
            content = f.read()

        if file.endswith(".embl") or file.endswith(".gb"):
            genomic_sequence, features =  parsers.parse_embl(content) if file.endswith(".embl") else parsers.parse_genbank(content)
            #is this genomic sequence already stored in the database? (can be useful if we have imported the FASTA file and then, afterwards, we want to import the annotations from an EMBL file)
            genome_stored = None
            for genome in db['genomes'].find():
                if genome['sequence'] == genomic_sequence:
                    genome_stored = genome
                    print "Genomic sequence %s already stored. I will only import the annotations"%genome_stored['name']
                    break
            if not genome_stored:
                genomic_sequence.source = get_file_as_source(file) 
                genomic_sequences.append(genomic_sequence)
                
                for (index, row) in features.iterrows():
                    annotation = {
                        '_id': str(ObjectId()),
                        'class': row['type'],
                        'source': genomic_sequence.source,
                        'organism': genomic_sequence.organism,
                        'genomicPositions': row['genomicPositions'],
                        'genomicStrand': row['genomicStrand'],
                        'genome': "%s@genomes"%genomic_sequence._id,
                        'genomeName': genomic_sequence.name
                    }

                    for key in row.keys():
                        if not key in ['type', 'genomicPositions', 'genomicStrand'] and not isnull(row[key]):
                            annotation[key] = row[key]

                    db['annotations'].insert(annotation)
            else:
                for (index, row) in features.iterrows():
                    annotation = {
                        '_id': str(ObjectId()),
                        'class': row['type'],
                        'source': genome_stored['source'],
                        'organism': genome_stored['organism'],
                        'genomicPositions': row['genomicPositions'],
                        'genomicStrand': row['genomicStrand'],
                        'genome': "%s@genomes"%genome_stored['_id'],
                        'genomeName': genome_stored['name']
                    }

                    for key in row.keys():
                        if not key in ['type', 'genomicPositions', 'genomicStrand'] and not isnull(row[key]):
                            annotation[key] = row[key]

                    db['annotations'].insert(annotation) 
    
    if not len(genomic_sequences):
        print "No genomic sequences found"

    for genomic_sequence in genomic_sequences:
        
        description = {
            "name": genomic_sequence.name,
            "_id": genomic_sequence._id,
            "sequence": genomic_sequence.sequence,
            "source": genomic_sequence.source
        }

        if genomic_sequence.organism:
            description["organism"] = genomic_sequence.organism
        
        db['genomes'].insert(description)

    client.disconnect()

    print "%i genomic sequences imported..."%len(genomic_sequences)

def remove_species(species, db_host = "localhost", db_port = 27017):
    """
    This method removes the database for the species given as argument and clean the comparative_genomics database.
    """
    client = MongoClient(db_host, db_port)
    db = client[species]

    #TO DO!!!!!
    
    #print "%i genomic sequences, %i annotations and %i ncRNAs have been removed. %i alignments have been curated!!"%(genomes_removed, annotations_removed, ncRNAs_removed, alignments_curated)

    client.disconnect()

def curate(db_host = "localhost", db_port = 27017):
    """
    This function fixes several features of the database:
    - the name of the sequences in the alignments. RFAM labels them as accession_number/pos1-pos2. This function renames them as organism_name/accession_number/pos1-pos2
    """
    annoying_families = []
    client = MongoClient(db_host, db_port)
    db = client['comparative_genomics']
    ncbi = NCBI()
    for alignment in db['alignments'].find(timeout = False):
        if alignment['source'] in annoying_families:
            continue
        print "Processing %s"%alignment['source']
        aligned_rnas, consensus_2d = parse_clustalw(alignment['alignment'])
        ids = []
        for aligned_rna in aligned_rnas:
            tokens = aligned_rna.name.split('/')
            if len(tokens) == 2 and len(tokens[1].split('-')) == 2: #we have an RFAM sequence name without any species name
                ids.append(tokens[0])
        ids = list(set(ids))
        if len(ids):
            print "%i sequence names to curate..."%len(ids)
            organisms = []
            if len(ids) > 50:
                print "chunks generated"
                from pyrna.utils import chunks
                for _ids in chunks(ids, 50):
                    print "%i sequence names to curate..."%len(_ids)
                    organisms += ncbi.get_organism(db='nucleotide', ids = _ids) 
                    print "Got the NCBI answer..."          
            else:
                organisms = ncbi.get_organism(db='nucleotide', ids = ids) 
                print "Got the NCBI answer..."
            if len(ids) != len(organisms):
                print 'Warning!! For alignment %s, %i species names to recover, got %i!!'%(alignment['source'], len(ids), len(organisms))
            else:
                for aligned_rna in aligned_rnas:
                    tokens = aligned_rna.name.split('/')
                    if len(tokens) == 2 and len(tokens[1].split('-')) == 2: #we have an RFAM sequence name without any species name
                        aligned_rna.name =  organisms[ids.index(tokens[0])]+"/"+aligned_rna.name
                alignment['alignment'] = to_clustalw(consensus_2d, aligned_rnas)
                db['alignments'].save(alignment)
        else:
            print "No sequence name to curate!!"

def check_if_fungal(ncbi, accession_number):
    return re.search("\s+Eukaryota; Fungi;", ncbi.efetch(ids = [accession_number], db='nucleotide', rettype='gb', header=500)) #the header param is to avoid to DL the full entry just to check the lineage

def get_fungal_genomic_entries_from_ncbi(ncbi, base_path):
    if not os.path.exists(base_path+'/ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/Eukaryotes/fungi/'):
        commands.getoutput('cd %s ; wget --retr-symlinks -qr "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/Eukaryotes/fungi/"'%base_path)

    species2accession_numbers = {}
    assemblies_count = 0
    assemblies = base_path+'/ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/Eukaryotes/fungi/'
    found_high_quality_assembly = False
    for species in os.listdir(assemblies):
        species_dir = os.path.join(assemblies, species)
        accession_numbers = set()
        for dir in os.listdir(species_dir):
            if dir == "representatives":
                found_high_quality_assembly = False
                dir = os.path.join(species_dir, dir)
                for f in os.listdir(dir):
                    if f.endswith(".assembly.txt"):
                        scaffold_count = 0
                        contig_count = 0
                        molecule_count = 0
                        organism = None
                        with open(os.path.join(dir,"%s.stats.txt"%f.split(".assembly.txt")[0])) as h:
                            lines = h.readlines()

                        for line in lines:
                            match = re.search("^all.+scaffold-count\s+([0-9]+).+$",line)
                            if match:
                                scaffold_count = int(match.group(1))
                            match = re.search("^all.+contig-count\s+([0-9]+).+$",line)
                            if match:
                                contig_count = int(match.group(1))
                            match = re.search("^all.+molecule-count\s+([0-9]+).+$",line)
                            if match:
                                molecule_count = int(match.group(1))
                            match = re.search("^#.+Organism name:\s+(.+)$",line)
                            if match:
                                organism = match.group(1)
                        if molecule_count:
                            assemblies_count += 1
                            if not found_high_quality_assembly:
                                found_high_quality_assembly = True
                            
                            with open(os.path.join(dir,f)) as h:
                                lines = h.readlines()

                            found_a_file = False
                            for line in lines:
                                if not line.startswith('#') and re.findall("Chromosome", line):
                                    tokens = line.split('\t')
                                    refseq = tokens[6]
                                    genbank = tokens[4]
                                    if refseq != 'na':
                                        accession_numbers.add(refseq)
                                    elif genbank != 'na':
                                        accession_numbers.add(genbank)
                        #uncomment the next block if you want to use assemblies made only with scaffolds (less quality)
                        # elif scaffold_count:
                        #     assemblies_count += 1
                        #     if not found_high_quality_assembly:
                        #         found_high_quality_assembly = True
                        
                        #     with open(os.path.join(dir,f)) as h:
                        #       lines = h.readlines()

                        #     found_a_file = False
                        #     for line in lines:
                        #         if not line.startswith('#') and re.findall("scaffold", line):
                        #              tokens = line.split('\t')
                        #              refseq = tokens[6]
                        #              genbank = tokens[4]
                        #              if refseq != 'na':
                        #                  accession_numbers.add(refseq)
                        #              elif genbank != 'na':
                        #                  accession_numbers.add(genbank)
        species2accession_numbers[species] = list(accession_numbers)
    print len(species2accession_numbers), "species recovered from the NCBI FTP"
    print assemblies_count, "assemblies with good quality recovered from the NCBI FTP"
    return species2accession_numbers

if __name__ == '__main__':
    db_host = "localhost"
    db_port = 27017
    species_to_remove = None
    species_to_add = None
    files = None
    do_dump = False
    do_drop = False
    do_phylogeny = False
    do_curation = False
    base_path = None

    if '-h' in sys.argv:
        db_host = sys.argv[sys.argv.index("-h")+1]
    if '-p' in sys.argv:
        db_port = int(sys.argv[sys.argv.index("-p")+1])
    if '-remove' in sys.argv:
        species_to_remove = sys.argv[sys.argv.index("-remove")+1]
    if '-add' in sys.argv:
        species_to_add = sys.argv[sys.argv.index("-add")+1]
        files = sys.argv[sys.argv.index("-add")+2:]
    if '-dump' in sys.argv:
        do_dump = True 
        if len(sys.argv)-1 == sys.argv.index("-dump")+1:
            source = sys.argv[sys.argv.index("-dump")+1]
    if '-dir' in sys.argv:
        base_path =  sys.argv[sys.argv.index("-dir")+1]
    do_curation = '-curate' in sys.argv
    do_drop = '-drop' in sys.argv
    do_phylogeny = '-phylo' in sys.argv 

    if do_phylogeny:
        to_newick_format(db_host = db_host, db_port = db_port)    
    elif do_dump:
        dump(db_host = db_host, db_port = db_port)
    elif species_to_add:
        add_species(species_name = species_to_add, files = files, db_host = db_host, db_port = db_port)
    elif species_to_remove:
        remove_species(species = species_to_remove, db_host = db_host, db_port = db_port)
    elif do_drop:
        drop(db_host = db_host, db_port = db_port)
    elif do_curation:
        curate(db_host = db_host, db_port = db_port)
    else:
        if not base_path:
            print "Usage: charnDB_mongo.py [-dir output_dir] [-h db_host (default: localhost)] [-p db_port (default: 27017)] [-remove database_name] [-add database_name genbank_file_1 embl_file_2 ...] [-dump] [-drop]"
            sys.exit(-1)
        create(base_path, db_host = db_host, db_port = db_port)

