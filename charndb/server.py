#!/usr/bin/env python

import uuid, ujson, sys, datetime, shutil, re, os, random, string, shlex, webbrowser, numpy as np, pickle

import tornado.httpserver
import tornado.ioloop
import tornado.options
import tornado.web
import tornado.websocket
from tornado.escape import json_encode, to_unicode

from pymongo import MongoClient
from bson.objectid import ObjectId
from bson.json_util import dumps
from pandas import DataFrame, isnull
from werkzeug import secure_filename
from subprocess import Popen

from pyrna.computations import Rnafold, Contrafold, Rnaplot, Rnaview, Mlocarna, Rnasubopt, RnaAlifold, Samtools
from pyrna.db import PDB, Rfam
from pyrna import parsers
from pyrna.parsers import parse_vienna, parse_fasta, parse_stockholm, base_pairs_to_secondary_structure, parse_pdb, to_clustalw, parse_clustalw, consensus2d_to_base_pairs, to_bn
from pyrna.features import RNA

mongodb = None
logs_db = None
enabled_algorithms = []
rfam = Rfam(cache_dir = "%s/RFAM"%os.getenv("HOME"))
families_details = rfam.get_families_details()
rfam.generate_seed_alignments()
counts = {}
system_databases = ['logs', 'PDB', 'RNA3DHub', 'local', 'test', 'comparative_genomics']

def get_consensus_2d(structural_alignment, junction_diameter = 20):
        
    ss_fungi = None
    fungi_ss_json = {}
    fungal_consensus_2d = None
    fungal_aligned_rnas = None
    
    fungal_aligned_rnas, fungal_consensus_2d = parse_clustalw(structural_alignment)
    ss_fungi = base_pairs_to_secondary_structure(fungal_aligned_rnas[0], fungal_consensus_2d)
    
    if ss_fungi:
        ss_fungi.find_junctions()
        
        for helix in ss_fungi.helices:
            sizes = []
            descriptions = []
            for aligned_rna in fungal_aligned_rnas:
                strand1 = aligned_rna.sequence[helix['location'][0][0]-1:helix['location'][0][-1]].replace('-','')
                sizes.append(len(strand1))
                strand2 = aligned_rna.sequence[helix['location'][-1][0]-1:helix['location'][-1][-1]].replace('-','')
                sizes.append(len(strand2))
                descriptions.append([aligned_rna.name, strand1, strand2])
            if len(sizes) > 2:
                helix['quantitative_value'] = np.std(sizes)
            helix['descriptions'] = descriptions

        for single_strand in ss_fungi.single_strands:
            sizes = []
            descriptions = []
            for aligned_rna in fungal_aligned_rnas:
                strand = aligned_rna.sequence[single_strand['location'][0]-1:single_strand['location'][-1]].replace('-','')
                sizes.append(len(strand))  
                descriptions.append([aligned_rna.name, strand])
            if len(sizes) > 1:
                single_strand['quantitative_value'] = np.std(sizes)
            single_strand['descriptions'] = descriptions

        for junction in ss_fungi.junctions:
            junction['location'].sort() #we need to be sure that the locations are sorted (so the display will be from the 5' to the 3' ends)
            sizes = []
            descriptions = []
            for aligned_rna in fungal_aligned_rnas:
                size = 0
                strands = []
                for single_strand in junction['location']:
                    strand = aligned_rna.sequence[single_strand[0]:single_strand[-1]-1].replace('-','')
                    size += len(strand)
                    strands.append(strand)
                sizes.append(size)
                strands.insert(0, aligned_rna.name)
                descriptions.append(strands)
            if len(sizes) > 1:
                junction['quantitative_value'] = np.std(sizes)
            junction['descriptions'] = descriptions

        junction_diameter = 20
        ss_fungi.compute_plot(step = 40, residue_occupancy = 10, junction_diameter = junction_diameter)

        #the secondary structures are dumped as JSON string

        #the fungal ss
        
        helices_descr = []
        for helix in ss_fungi.helices:
            descr = {
                'name': helix['name'],
                'location': helix['location'],
                'coords': helix['coords'],
                'descriptions': helix['descriptions']
            }
            if helix.has_key('quantitative_value'):
                descr['quantitative_value'] = helix['quantitative_value']
            helices_descr.append(descr) 
            
        fungi_ss_json['helices'] = helices_descr

        single_strands_descr = []
        for single_strand in ss_fungi.single_strands:
            descr =  {
                    'name': single_strand['name'],
                    'location': single_strand['location'],
                    'descriptions': single_strand['descriptions']
                }

            if single_strand.has_key('coords'):
                descr['coords'] = single_strand['coords']
            if single_strand.has_key('quantitative_value'):
                descr['quantitative_value'] = single_strand['quantitative_value']

            single_strands_descr.append(descr) 

        fungi_ss_json['single_strands'] = single_strands_descr

        junctions_descr = []
        for junction in ss_fungi.junctions:
            descr = {
                'location': junction['location'],
                'coords': junction['coords'],
                'descriptions': junction['descriptions']
            }

            if junction.has_key('quantitative_value'):
                descr['quantitative_value'] = junction['quantitative_value']
            
            junctions_descr.append(descr)

        fungi_ss_json['junctions'] = junctions_descr

        from pyrna.utils import get_points
        diagonals_descr = []
        for junction in ss_fungi.junctions:
            if len(junction['location']) >= 3:
                junction_location = sorted(junction['location'])
                for i in range(len(junction_location)-1):
                    for h in ss_fungi.helices: #next helices in junction
                        if h['location'][0][0] == junction_location[i][-1]:
                            if h['coords'][0][1] != junction['coords'][0][1]: #to avoid to redraw a vertical line
                                new_points = get_points(h['coords'][0][0], h['coords'][0][1], junction['coords'][0][0], junction['coords'][0][1], distance = (junction_diameter+10)/2)
                                if len(new_points) == 2:
                                    diagonal = {}
                                    if h.has_key('quantitative_value'):
                                        diagonal['quantitative_value'] = h['quantitative_value']
                                    diagonal['coords'] = [
                                            [h['coords'][0][0], h['coords'][0][1]], 
                                            [new_points[1][0], new_points[1][1]]
                                            ]
                                    diagonal['descriptions'] = h['descriptions']
                                    diagonals_descr.append(diagonal)
        fungi_ss_json['diagonals'] = diagonals_descr

        directly_linked_helices_descr = []
        previous_helix = ss_fungi.helices[0]
        currentPos = previous_helix['location'][-1][-1]
        while currentPos <= len(ss_fungi.rna):
            currentPos +=1
            for helix in ss_fungi.helices:
                if currentPos == helix['location'][0][0]:
                    if previous_helix['location'][-1][-1] +1 == helix['location'][0][0]:
                        helix = {
                            'coords': [
                                [previous_helix['coords'][0][0], previous_helix['coords'][0][1]],
                                [helix['coords'][0][0], helix['coords'][0][1]]
                            ]
                        }
                        directly_linked_helices_descr.append(helix)
                    currentPos = helix['location'][-1][-1]
                    previous_helix = helix
                    break

        fungi_ss_json['directly-linked-helices'] = directly_linked_helices_descr

        return fungi_ss_json, to_bn(fungal_consensus_2d, len(fungal_aligned_rnas[0]))

def get_ortholog_2d(structural_alignment, ortholog_name):
    ss_json = None
    orthologs = []
    aligned_sequences = []
    
    fungal_aligned_rnas, fungal_consensus_2d = parse_clustalw(structural_alignment)

    fungal_aligned_rnas = sorted(fungal_aligned_rnas, key = lambda fungal_aligned_rna: fungal_aligned_rna.name)

    aligned_rna = fungal_aligned_rnas[0]
    for fungal_aligned_rna in fungal_aligned_rnas:
        orthologs.append(fungal_aligned_rna.name)
        aligned_sequences.append(fungal_aligned_rna.sequence)
        if ortholog_name and ortholog_name == fungal_aligned_rna.name:
            aligned_rna = fungal_aligned_rna    
    
    rna = RNA(name = aligned_rna.name, sequence = aligned_rna.sequence.replace('-',''))

    base_pairs = consensus2d_to_base_pairs(aligned_rna, fungal_consensus_2d)
    ss = base_pairs_to_secondary_structure(rna, base_pairs)
    ss.find_junctions()

    rnaplot = Rnaplot()
    coords = rnaplot.plot(base_pairs, rna)

    ss_json = {
        'rna': {
            'name': rna.name,
            'sequence': rna.sequence
        }
    }

    helices_descriptions = []
    for helix in ss.helices:
        helices_descriptions.append(helix)
    ss_json['helices'] = helices_descriptions

    sstrand_descriptions = []
    for sstrand in ss.single_strands:
        sstrand_descriptions.append(sstrand)
    ss_json['single-strands'] = sstrand_descriptions

    tertiary_descriptions = []
    for tertiary in ss.tertiary_interactions:
        sstrand_descriptions.append(tertiary)
    ss_json['single-strands'] = sstrand_descriptions

    ss_json['coords'] = []
    for index, row in coords.iterrows():
        ss_json['coords'].append([int(row['x']), int(row['y'])]) 
    return ss_json, orthologs, aligned_sequences  

class IndexHandler(tornado.web.RequestHandler):
    def get(self):
        self.render('index.html')

class GenomeBrowserHandler(tornado.web.RequestHandler):
    def get(self):
        species = self.get_argument('species', None)
        genome_id = self.get_argument('genome_id', None)
        annotation_id = self.get_argument('annotation_id', '')
        start = self.get_argument('start', None)
        end = self.get_argument('end', None)

        species_list = [db for db in mongodb.database_names() if db not in system_databases]
        species_list.sort()

        if species is None:
            species = species_list[0]
        else:
            species = self.decode_argument(species.replace(" ", "_"))

        genome_list = mongodb[species]['genomes'].find({}, {"_id":1, "name":1})

        if genome_id is None or species is None:
            genome_id = genome_list[0]["_id"]

        if start is None and end is None:
            start = 1
            end = len(mongodb[species]['genomes'].find_one({"_id": genome_id}, {"sequence":1})['sequence'])
            center = (end+start)/2
            start = center -5000
            end = center + 5000

        self.render('browser.html',
            species = species,
            species_list = species_list,
            genome_id = genome_id,
            genome_list = genome_list,
            annotation_id = annotation_id,
            start = start,
            end = end)

class GetDataHandler(tornado.web.RequestHandler):
    def get(self):
        self.render('data.html')

class HelpHandler(tornado.web.RequestHandler):
    def get(self):
        self.render('help.html')

class FamiliesHandler(tornado.web.RequestHandler):
    def get(self):
        species_names = []
        #species_names = self.get_argument('species')       
        #print species_names
        self.render('families.html', species_names = species_names)

    def post(self):
        species_names = []
        #species_names = self.get_argument('species')       
        #print species_names
        self.render('families.html', species_names = species_names)

class SpeciesHandler(tornado.web.RequestHandler):
    def get(self):
        self.render('species.html')

class SyntheniesHandler(tornado.web.RequestHandler):
    def get(self):
        ncRNA_name = self.get_argument('ncRNA_name', None)
        self.render('synthenies.html', ncRNA_name = ncRNA_name, ncRNA_description = families_details[families_details['id'] == ncRNA_name].iloc[0,1])

class S2SHandler(tornado.web.RequestHandler):
    def get(self):
        junction_diameter = 20
        ncRNA_name = self.get_argument('ncRNA_name', None)
        aligned_sequences = None
        fungal_bn = None
        ss_json = None
        fungi_ss_json = None
        orthologs = None
        if ncRNA_name and os.path.exists('static/data/alignments/%s_fungi.aln'%ncRNA_name):
            with open('static/data/alignments/%s_fungi.aln'%ncRNA_name) as h:
                structural_alignment = h.read()
                ss_json, orthologs, aligned_sequences = get_ortholog_2d(structural_alignment, None)
                fungi_ss_json, fungal_bn = get_consensus_2d(structural_alignment, junction_diameter)
        self.render('s2s.html', ncRNA_name = ncRNA_name, orthologs = orthologs, aligned_sequences = aligned_sequences, bn = fungal_bn, ncRNA_description = families_details[families_details['id'] == ncRNA_name].iloc[0,1], secondary_structure = json_encode(ss_json), ss_fungi = json_encode(fungi_ss_json), junction_diameter = junction_diameter)

class BooquetHandler(tornado.web.RequestHandler):
    def get(self):
        self.render('booquet.html', rfam_families_details = families_details.to_json())

class BooquetSampleHandler(tornado.web.RequestHandler):
    def get(self):
        with open('static/data/alignments/Afu_182_fungi.aln') as h:
            structural_alignment = h.read()
            self.set_header("Content-Type", "text/plain")
            self.set_header('Content-Disposition', 'attachment; filename="booquet_sample"')
            self.write(structural_alignment)

class ConsensusStructureHandler(tornado.web.RequestHandler):
    def get(self):
        ncRNA_name = self.get_argument('ncRNA_name', None)
        junction_diameter = 20
        fungi_ss_json, fungal_bn = get_consensus_2d(ncRNA_name, junction_diameter)
        self.render('consensus_structure.html', ncRNA_name = ncRNA_name, ncRNA_description = families_details[families_details['id'] == ncRNA_name].iloc[0,1], ss_fungi = json_encode(fungi_ss_json), junction_diameter = junction_diameter)

class LetItSnoHandler(tornado.web.RequestHandler):

    ALLOWED_EXTENSIONS = set(['txt', 'fasta', 'fna', 'embl', 'gb'])

    def letitsno():
        log = {
            'path': request.path,
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }
        logs_db['websites'].insert(log)
        return render_template('letitsno.html')

    def allowed_file(filename):
        return '.' in filename and \
            filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

    #@app.route('/letitsno/upload', methods=['POST'])
    def upload():
        log = {
            'path': request.path,
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }
        logs_db['websites'].insert(log)
        try:
            file = request.files['file']
            organism = request.form['organism']
            if file and allowed_file(file.filename):
                secret_key = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(10))
                
                import_genomic_sequences(file.read(), organism, file.filename.rsplit('.', 1)[1], secret_key)

                file.close()

                command_line = "let_it_sno.task -db %s -h arn-ibmc.in2p3.fr -algorithms /opt/exp_soft/grand-est/ARN/pyrna_algorithms_2 -python /opt/exp_soft/grand-est/ARN/epd -e 2"%secret_key
                outputDir = os.getenv("HOME")+"/tmp/jobs_let_it_sno.task_on_%s"%secret_key
                
                shutil.os.mkdir(outputDir)

                fh = open("NUL",'w') #to disable console output
                p = Popen(shlex.split(command_line), stdout = fh, stderr = fh)

                with open(os.getenv("HOME")+"/tmp/jobs_let_it_sno.task_on_%s/info"%secret_key,'a+') as info:
                    info.write("process ID: %i\n"%p.pid)

                return jsonify({"secret_key": secret_key, "annotations_count": mongodb[secret_key]['annotations'].count(), 'genomic_sequences_count': mongodb[secret_key]['genomes'].count()} )

        except:
            return jsonify({'status': 'error'})

    def import_genomic_sequences(genomic_data, organism, suffix, db_name):
        db = mongodb[db_name]
        genomic_sequences = []

        if suffix in ['embl', 'gb']:
            genomic_sequence, features =  parsers.parse_embl(genomic_data) if suffix == "embl" else parsers.parse_genbank(genomic_data)
            
            genomic_sequences.append(genomic_sequence)
            
            for (index, row) in features.iterrows():
                annotation = {
                    '_id': str(ObjectId()),
                    'class': row['type'],
                    'source': genomic_sequence.source,
                    'organism': organism,
                    'genomicPositions': row['genomicPositions'],
                    'genomicStrand': row['genomicStrand'],
                    'genome': "%s@genomes"%genomic_sequence._id,
                    'genomeName': genomic_sequence.name,
                }

                for key in row.keys():
                    if not key in ['type', 'genomicPositions', 'genomicStrand'] and not isnull(row[key]):
                        annotation[key] = row[key]

                db['annotations'].insert(annotation)

        elif suffix in ['fasta', 'fna']:
            genomic_sequences += parsers.parse_fasta(genomic_data, type='DNA')

        for genomic_sequence in genomic_sequences:
            
            description = {
                "name": genomic_sequence.name,
                "_id": genomic_sequence._id,
                 "organism":organism,
                "sequence": genomic_sequence.sequence,
                "source": genomic_sequence.source
            }           
            
            db['genomes'].insert(description)

    #@app.route('/letitsno/results')
    def results():
        log = {
            'path': request.path,
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }
        logs_db['websites'].insert(log)
        import os.path
        import scipy.stats

        if os.path.isfile(os.getenv("HOME")+'/tmp/jobs_let_it_sno.task_on_%s/info'%request.args['id']):
            db = mongodb[request.args['id']]
            all_scores = {'snoscan': [], 'snogps': [], 'snoreport': [], 'rnamotif': [], 'cmsearch': []}
            
            for ncRNA in db['ncRNAs'].find({'source': {'$ne': 'tool:letitsno:NA'}}):
                if ncRNA['source'] != 'tool:cmsearch:NA':
                    if ncRNA.has_key('score'): #some hits found by snoscan have not the key 'score'
                        all_scores[ncRNA['source'].split(':')[1]].append(ncRNA['score'])
                else:
                    all_scores['cmsearch'].append(1/ncRNA['score']) #ncRNA['score'] is an e-value

            #### WE CONSTRUCT THE DATAFRAME FROM ALL NCRNA CLUSTERS #########
            snoclusters = []

            genome_dict = {}
            for genome in db["genomes"].find():
                genome_dict[genome['_id']+"@genomes"] = genome['name']

            cd_cmsearch = 0
            cd_rnamotif = 0
            cd_snoscan = 0
            cd_snoreport = 0
            cd_snoreport_cmsearch = 0
            cd_snoreport_rnamotif = 0
            cd_snoreport_snoscan = 0
            cd_cmsearch_rnamotif = 0
            cd_cmsearch_snoscan = 0
            cd_rnamotif_snoscan = 0
            cd_snoreport_cmsearch_rnamotif = 0
            cd_snoreport_cmsearch_snoscan = 0
            cd_snoreport_rnamotif_snoscan = 0
            cd_cmsearch_rnamotif_snoscan = 0
            cd_all = 0

            haca_cmsearch = 0
            haca_rnamotif = 0
            haca_snogps = 0
            haca_snoreport = 0
            haca_snoreport_cmsearch = 0
            haca_snoreport_rnamotif = 0
            haca_snoreport_snogps = 0
            haca_cmsearch_rnamotif = 0
            haca_cmsearch_snogps = 0
            haca_rnamotif_snogps = 0
            haca_snoreport_cmsearch_rnamotif = 0
            haca_snoreport_cmsearch_snogps = 0
            haca_snoreport_rnamotif_snogps = 0
            haca_cmsearch_rnamotif_snogps = 0
            haca_all = 0

            #now we want to make two thinks:
            # 1 - we want to describe each cluster as a JSON object. For the description, we will compute the percentile of the best scores for each tool that participated to the cluster.
            # 2 - we want to count the number of clusters for each category highlighted in the upcoming venn diagrams.

            for cluster in db['ncRNAs'].find({'source': "tool:letitsno:NA"}): #in the ncRNAs table, the ncRNAs predictions have the original tool as source (snoscan, snoGPS), and the clusters have letitsno as source.
                snocluster = {
                    '_id': cluster['_id'],
                    'class': cluster['class'].split(', ')[-1][:-4],
                    'genome': genome_dict[cluster['genome']],
                    'strand': cluster['genomicStrand'],
                    'start': cluster['genomicPositions'][0],
                    'end': cluster['genomicPositions'][1]
                }
                cluster_names = {} #this will be used to construct the name of the cluster: snoRNA name_1 (list of tools that produced a prediction with this name), snoRNA name_2 (list of tools that produced a prediction with this name),... 
                cluster_scores = {} #This will be used to compute the best score and its percentile for each tool. If a cluster adds a key corresponding to the name of a tool, then this cluster contains at least one annotation produced from this tool. 
                                    #It should be noted that the array of scores can be empty if the tool has produced no score for its annotations (like snoscan for example, for which some scores are so weird that its wrapper in db.computations return no score.)
                
                for _id in cluster['ids']:
                    ncRNA = db['ncRNAs'].find_one({'source': {'$ne': "tool:letitsno:NA"}, '_id': _id[:-7]}) #_id = id+"@ncRNAs"
                    source = ncRNA['source'].split(':')[1]
                    if not cluster_scores.has_key(source):
                        cluster_scores[source] = []
                    #a same name can be found by several tools, and differents names can be found by the same tool!
                    if ncRNA.has_key('name'): #key not always present
                        if not cluster_names.has_key(ncRNA['name']): #no redundancy between names
                            cluster_names[ncRNA['name']] = [source]
                        else:
                            if source not in cluster_names[ncRNA['name']]: #no redundancy between tools for the same name
                                cluster_names[ncRNA['name']].append(source)
                    if ncRNA.has_key('score'): #some hits found by snoscan have not the key 'score'
                        cluster_scores[source].append(ncRNA['score'])

                snocluster['tools'] = cluster_scores.keys()
                
                all_names = "Undefined"
                
                #construction of the name of the cluster based on the different names produced by the tools that participated to the cluster
                if len(cluster_names):
                    all_names = ""
                    for item in cluster_names.iteritems():
                        all_names += "%s(%s"%(item[0], item[1][0])
                        if len(item[1]) > 1:
                            for tool in item[1][1:]:
                                all_names += ",%s"%tool
                        all_names += "),"
                    snocluster['name'] = all_names[:-1]

                snocluster['name'] = all_names

                #computation of the best score for each tool and of the percentile of this best score. Percentile X for score A means that X% of the scores are below score A.
                for item in cluster_scores.iteritems():
                    if len(item[1]):
                        best_score = max(item[1]) if item[0] != 'cmsearch' else 1/min(item[1]) #score given by cmsearch is an e-value X, meaning that the lowest value is the best. We need to homogenize with the other systems. So we compute 1/X
                        snocluster[item[0]+"_best_score"] = best_score
                        percentile = scipy.stats.percentileofscore(all_scores[item[0]], best_score)
                        snocluster[item[0]+"_percentile"] = percentile
                    else: #if no score, then we decide that the percentile is 0% (meaning that 0% of the scores are below these "no-score", which is logically acceptable)
                        snocluster[item[0]+"_best_score"] = -1
                        snocluster[item[0]+"_percentile"] = 0

                #now we check to which categories of the upcoming venn diagram this cluster will fit in. If the dict cluster_scores has as key corresponding to the name of a tool, then this cluster contains at least one ncRNA prediction from this tool (even if the array of scores is empty) 
                if cluster_scores.has_key('snoscan'):
                    if cluster_scores.has_key('snoreport'):
                        if cluster_scores.has_key('rnamotif'):
                            if cluster_scores.has_key('cmsearch'):
                                cd_all  +=1     
                            else:
                                cd_snoreport_rnamotif_snoscan += 1
                        elif cluster_scores.has_key('cmsearch'):
                                cd_snoreport_cmsearch  +=1        
                        else:
                            cd_snoreport_snoscan += 1    

                    elif cluster_scores.has_key('rnamotif'):
                            if cluster_scores.has_key('cmsearch'):
                                cd_cmsearch_rnamotif_snoscan += 1
                            else: 
                                cd_rnamotif_snoscan +=1 
                    elif cluster_scores.has_key('cmsearch'):
                        cd_cmsearch_snoscan += 1    
                    else:
                        cd_snoscan += 1
                elif cluster_scores.has_key('snoreport'):
                    if cluster_scores.has_key('rnamotif'):
                        if cluster_scores.has_key('cmsearch'):
                            cd_snoreport_cmsearch_rnamotif +=1     
                        else:
                            cd_snoreport_rnamotif += 1
                    elif cluster_scores.has_key('cmsearch'):
                            cd_snoreport_cmsearch +=1    
                    else:
                        cd_snoreport += 1
                elif cluster_scores.has_key('rnamotif'):
                    if cluster_scores.has_key('cmsearch'):
                        cd_cmsearch_rnamotif +=1     
                    else:
                        cd_rnamotif += 1 
                elif cluster_scores.has_key('cmsearch'):
                    cd_cmsearch += 1        

                if cluster_scores.has_key('snogps'):
                    if cluster_scores.has_key('snoreport'):
                        if cluster_scores.has_key('rnamotif'):
                            if cluster_scores.has_key('cmsearch'):
                                haca_all  +=1     
                            else:
                                haca_snoreport_rnamotif_snogps += 1
                        elif cluster_scores.has_key('cmsearch'):
                                haca_snoreport_cmsearch  +=1        
                        else:
                            haca_snoreport_snogps += 1    

                    elif cluster_scores.has_key('rnamotif'):
                            if cluster_scores.has_key('cmsearch'):
                                haca_cmsearch_rnamotif_snogps += 1
                            else: 
                                haca_rnamotif_snogps +=1 
                    elif cluster_scores.has_key('cmsearch'):
                        haca_cmsearch_snogps += 1    
                    else:
                        haca_snogps += 1
                elif cluster_scores.has_key('snoreport'):
                    if cluster_scores.has_key('rnamotif'):
                        if cluster_scores.has_key('cmsearch'):
                            haca_snoreport_cmsearch_rnamotif +=1     
                        else:
                            haca_snoreport_rnamotif += 1
                    elif cluster_scores.has_key('cmsearch'):
                            haca_snoreport_cmsearch +=1    
                    else:
                        haca_snoreport += 1
                elif cluster_scores.has_key('rnamotif'):
                    if cluster_scores.has_key('cmsearch'):
                        haca_cmsearch_rnamotif +=1     
                    else:
                        haca_rnamotif += 1 
                elif cluster_scores.has_key('cmsearch'):
                    haca_cmsearch += 1  
                
                snoclusters.append(snocluster)
            
            return jsonify({
                'clusters': snoclusters,
                'cd_cmsearch': cd_cmsearch,
                'cd_rnamotif': cd_rnamotif,
                'cd_snoscan': cd_snoscan,
                'cd_snoreport': cd_snoreport,
                'cd_snoreport_cmsearch': cd_snoreport_cmsearch,
                'cd_snoreport_rnamotif': cd_snoreport_rnamotif,
                'cd_snoreport_snoscan': cd_snoreport_snoscan,
                'cd_cmsearch_rnamotif': cd_cmsearch_rnamotif,
                'cd_cmsearch_snoscan': cd_cmsearch_snoscan,
                'cd_rnamotif_snoscan': cd_rnamotif_snoscan,
                'cd_snoreport_cmsearch_rnamotif': cd_snoreport_cmsearch_rnamotif,
                'cd_snoreport_cmsearch_snoscan': cd_snoreport_cmsearch_snoscan,
                'cd_snoreport_rnamotif_snoscan': cd_snoreport_rnamotif_snoscan,
                'cd_cmsearch_rnamotif_snoscan': cd_cmsearch_rnamotif_snoscan,
                'cd_all':cd_all,
                'haca_cmsearch': haca_cmsearch,
                'haca_rnamotif': haca_rnamotif,
                'haca_snogps': haca_snogps,
                'haca_snoreport': haca_snoreport,
                'haca_snoreport_cmsearch': haca_snoreport_cmsearch,
                'haca_snoreport_rnamotif': haca_snoreport_rnamotif,
                'haca_snoreport_snogps': haca_snoreport_snogps,
                'haca_cmsearch_rnamotif': haca_cmsearch_rnamotif,
                'haca_cmsearch_snogps': haca_cmsearch_snogps,
                'haca_rnamotif_snogps': haca_rnamotif_snogps,
                'haca_snoreport_cmsearch_rnamotif': haca_snoreport_cmsearch_rnamotif,
                'haca_snoreport_cmsearch_snogps': haca_snoreport_cmsearch_snogps,
                'haca_snoreport_rnamotif_snogps': haca_snoreport_rnamotif_snogps,
                'haca_cmsearch_rnamotif_snogps': haca_cmsearch_rnamotif_snogps,
                'haca_all':haca_all
            })        

    #@app.route('/letitsno/check')
    def check():
        log = {
            'path': request.path,
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }
        logs_db['websites'].insert(log)
        #we need to check if a submission file is available
        import os.path
        total_jobs_to_submit = 0
        if os.path.isfile(os.getenv("HOME")+'/tmp/jobs_let_it_sno.task_on_%s/info'%request.args['id']):
            with open(os.getenv("HOME")+'/tmp/jobs_let_it_sno.task_on_%s/info'%request.args['id']) as h:
                info_content = h.read()

            groups = re.findall("total jobs: (\d+)", info_content)

            if len(groups):
                total_jobs_to_submit = int(groups[0])
            else:
                return jsonify({
                    'total_jobs_to_submit': 0, 
                    'snoscan_jobs_done' : -1, 
                    'snogps_jobs_done': -1, 
                    'snoreport_jobs_done': -1, 
                    'cmsearch_jobs_done': -1, 
                    'rnamotif_jobs_done': -1, 
                    'clusters_refinement_done': -1
                })        

        else:
            return jsonify({
                'total_jobs_to_submit': -1, 
                'snoscan_jobs_done' : -1, 
                'snogps_jobs_done': -1, 
                'snoreport_jobs_done': -1, 
                'cmsearch_jobs_done': -1, 
                'rnamotif_jobs_done': -1, 
                'clusters_refinement_done': -1
            })
        
        db = mongodb[request.args['id']]
        snoscan_jobs_done = db['computations'].find({"source": {'$regex':'tool:snoscan:.+'} }).count()
        snogps_jobs_done = db['computations'].find({"source": {'$regex':'tool:snogps:.+'} }).count()
        snoreport_jobs_done = db['computations'].find({"source": {'$regex':'tool:snoreport:.+'} }).count()
        cmsearch_jobs_done = db['computations'].find({"source": {'$regex':'tool:cmsearch:.+'} }).count()
        rnamotif_jobs_done = db['computations'].find({"source": {'$regex':'tool:rnamotif:.+'} }).count()
        clusters_refinement_done = 0 
        if db['computations'].find({"source": {'$regex':'tool:check_letitsno_task:.+'} }).count() == 1:
            clusters_refinement_done = 100
        return jsonify({
            'total_jobs_to_submit': total_jobs_to_submit, 
            'snoscan_jobs_done' : int(float(snoscan_jobs_done)/float(total_jobs_to_submit)*100), 
            'snogps_jobs_done': int(float(snogps_jobs_done)/float(total_jobs_to_submit)*100), 
            'snoreport_jobs_done': int(float(snoreport_jobs_done)/float(total_jobs_to_submit)*100), 
            'cmsearch_jobs_done': int(float(cmsearch_jobs_done)/float(total_jobs_to_submit)*100), 
            'rnamotif_jobs_done': int(float(rnamotif_jobs_done)/float(total_jobs_to_submit)*100), 
            'clusters_refinement_done': clusters_refinement_done
            })

    #@app.route('/letitsno/documentation')
    def documentation():
        log = {
            'path': request.path,
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }
        logs_db['websites'].insert(log)
        return render_template('letitsno/documentation.html')

    #@app.route('/letitsno/stats')
    def stats():
        log = {
            'path': request.path,
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }
        logs_db['websites'].insert(log)
        return render_template('letitsno/stats.html')

    #@app.route('/letitsno/example')
    def example():
        log = {
            'path': request.path,
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }
        logs_db['websites'].insert(log)
        return render_template('letitsno/example.html')

class WebservicesHandler(tornado.web.RequestHandler):

    #@app.route("/webservices")
    #@app.route('/api/')
    def webservices():
        return render_template('webservices.html')

    #here starts the low-level webservices: webservices to avoid to install RNA algorithms on the client computer to be able to use PyRNA

    #@app.route('/api/get_key', methods=['GET'])
    def get_key():
        secret_key = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(10))
        ip = request.remote_addr
        record = logs_db['user_keys'].find_one({'ip': ip})
        if record:
            record['key'] = secret_key
            logs_db['user_keys'].save(record)
        else:
            logs_db['user_keys'].insert({
                'ip': ip,
                'date': str(datetime.datetime.now()),
                'key': secret_key
                })
        return secret_key

    def is_registered_user(secret_key):
        return logs_db['user_keys'].find_one({'key': secret_key}) != None

    #@app.route('/api/computations/rnafold', methods=['POST'])
    def rnafold():
        log = {
            '_id': str(ObjectId()),
            'path': request.path,
            'tool': 'rnafold',
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }

        logs_db['webservices'].insert(log)

        name = None
        sequence = None
        constraints = None
        api_key = None
        if 'name' in request.form: 
            name = request.form['name']
        if 'sequence' in request.form:
            sequence = request.form['sequence']
        if 'constraints' in request.form:
            constraints = request.form['constraints']
        if 'api_key' in request.form:
            api_key = request.form['api_key']

        if not is_registered_user(api_key) or not 'rnafold' in enabled_algorithms:
            return abort(401)

        return Rnafold().fold(RNA(name=name, sequence=sequence), constraints, raw_output = True)

    #@app.route('/api/computations/rnaplot', methods=['POST'])
    def rnaplot():
        log = {
            '_id': str(ObjectId()),
            'path': request.path,
            'tool': 'rnaplot',
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }

        logs_db['webservices'].insert(log)

        secondary_structure = None
        api_key = None
        if 'secondary_structure' in request.form: 
            secondary_structure = request.form['secondary_structure']
        if 'api_key' in request.form:
            api_key = request.form['api_key']

        if not is_registered_user(api_key) or not 'rnaplot' in enabled_algorithms:
            return abort(401)

        rnas, secondary_structures = parse_vienna(secondary_structure)

        return Rnaplot().plot(secondary_structures[0], rnas[0], raw_output = True)

    #@app.route('/api/computations/contrafold', methods=['POST'])
    def contrafold():
        log = {
            '_id': str(ObjectId()),
            'path': request.path,
            'tool': 'contrafold',
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }

        logs_db['webservices'].insert(log)

        name = None
        sequence = None
        constraints = None
        api_key = None
        if 'name' in request.form: 
            name = request.form['name']
        if 'sequence' in request.form:
            sequence = request.form['sequence']
        if 'api_key' in request.form:
            api_key = request.form['api_key']

        if not is_registered_user(api_key) or not 'contrafold' in enabled_algorithms:
            return abort(401)

        return Contrafold().fold(RNA(name=name, sequence=sequence), raw_output = True)

    #@app.route('/api/computations/rnaview', methods=['POST'])
    def rnaview():
        log = {
            '_id': str(ObjectId()),
            'path': request.path,
            'tool': 'rnaview',
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }

        logs_db['webservices'].insert(log)

        tertiary_structure = None
        canonical_only = False
        api_key = None
        
        if '3d' in request.form:
            tertiary_structure = request.form['3d']
        if 'canonical_only' in request.form:
            canonical_only = request.form['canonical_only'] == 'true'
        if 'api_key' in request.form:
            api_key = request.form['api_key']

        if not is_registered_user(api_key) or not 'rnaview' in enabled_algorithms:
            return abort(401)

        return Rnaview().annotate(pdb_content = tertiary_structure, canonical_only = canonical_only, raw_output = True)

    #here starts the high-level webservices

    #@app.route('/api/compute/2d', methods=['GET', 'POST'])
    def compute_2d():
        data = None
        tool = None
        version = 1
        pdbid = None
        output = None
        if request.method == 'POST':
            if 'data' in request.form:
                data = request.form['data']
            if 'tool' in request.form:
                tool = request.form['tool']
            if 'output' in request.form:
                output = request.form['output']
            if 'version' in request.form:
                version = request.form['version']
            if 'pdbid' in request.form:
                pdbid = request.form['pdbid'] 
        else:
            if 'data' in request.args:
                data = request.args.get('data', None)
            if 'tool' in request.args:
                tool = request.args.get('tool', None)
            if 'output' in request.args:
                output = request.args.get('output', None)
            if 'version' in request.args:
                version = request.args.get('version', 1)
            if 'pdbid' in request.args:
                pdbid = request.args.get('pdbid', None)
        result = None

        log = {
            'path': request.path,
            'tool': tool,
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }

        logs_db['webservices'].insert(log)

        if data and data.startswith('>'): #2D prediction
            rnas = parse_fasta(data)
            result = []
            if len(rnas) == 1: #single molecule prediction (MFE,...)
                rna = rnas[0]
                secondary_structures = []
                if tool == 'rnafold':
                    secondary_structures.append(base_pairs_to_secondary_structure(rna, Rnafold().fold(rna)))
                elif tool == 'contrafold':
                    secondary_structures.append(base_pairs_to_secondary_structure(rna, Contrafold().fold(rna)))
                elif tool == 'rnasubopt':
                    random_sample = 20
                    if request.method == 'POST':
                        random_sample = int(request.form['random_sample'])
                    else:
                        random_sample = int(request.args.get('random_sample', 20))    
                    for _result in Rnasubopt().fold(rna, random_sample = random_sample):
                        secondary_structures.append(base_pairs_to_secondary_structure(rna, _result))
                for ss in secondary_structures:
                    _result = {
                        '_id': ss._id,
                        'name': ss.name,
                        'source': ss.source,
                        'rna': {
                            'name': ss.rna.name,
                            'sequence': ss.rna.sequence,
                            'source': ss.rna.source,
                            '_id': ss.rna._id
                        }    
                    }

                    helices_descr = []
                    for helix in ss.helices:
                        helix_desc = {
                            'name': helix['name'],
                            'location': {'ends': helix['location']} if version == 1 else helix['location']
                        }
                        if helix.has_key('interactions'):
                            interactions_descr = []
                            for interaction in helix['interactions']:
                                interactions_descr.append({
                                    'orientation': interaction['orientation'],
                                    'edge1': interaction['edge1'],
                                    'edge2': interaction['edge2'],
                                    'location': {'ends': interaction['location']} if version == 1 else interaction['location']
                                })
                            helix_desc['interactions'] = interactions_descr   
                        
                        helices_descr.append(helix_desc)

                    _result['helices'] = helices_descr

                    single_strands_descr = []
                    for single_strand in ss.single_strands:
                        single_strands_descr.append({
                            'name': single_strand['name'],
                            'location': {'ends': single_strand['location']} if version == 1 else single_strand['location']
                        })

                    _result['singleStrands'] = single_strands_descr 

                    tertiary_interactions_descr = []
                    for tertiary_interaction in ss.tertiary_interactions:
                        tertiary_interactions_descr.append({
                            'orientation': tertiary_interaction['orientation'],
                            'edge1': tertiary_interaction['edge1'],
                            'edge2': tertiary_interaction['edge2'],
                            'location': {'ends': tertiary_interaction['location']} if version == 1 else tertiary_interaction['location']
                        })

                    _result['tertiaryInteractions'] = tertiary_interactions_descr
                    result.append(_result) 
            elif len(rnas) >= 2: #structural alignment
                if tool == 'mlocarna':
                    aligned_molecules, consensus2D = Mlocarna().align(rnas)
                    return Response(to_clustalw(consensus2D, aligned_molecules), mimetype='application/txt')

            if tool == 'rnafold' or tool == 'contrafold':
                return Response(ujson.dumps(result[0]), mimetype='application/json')
            else:
                return Response(ujson.dumps(result), mimetype='application/json')
        elif tool == 'rnalifold' and data and data.startswith('CLUSTAL'): #computation of consensus structure from sequence alignment
            return Response(RnaAlifold().align(data), mimetype='application/txt')
        elif tool == 'rnaview': #3D annotation
            rnaview = Rnaview()

            if output == 'rnaml':
                if pdbid:
                    return Response(rnaview.annotate(pdb_content = PDB().get_entry(pdbid), raw_output = True), mimetype='application/txt')
                elif data:
                    return Response(rnaview.annotate(pdb_content = data, raw_output = True), mimetype='application/txt')
                
            else:
                if pdbid:
                    tertiary_structures = parse_pdb(PDB().get_entry(pdbid))
                elif data:
                    tertiary_structures = parse_pdb(data)

                result = []

                for ts in tertiary_structures:

                    (ss, ts) = rnaview.annotate(ts, canonical_only = False)

                    ss.find_junctions()

                    _2D_descr = {
                        '_id': ss._id,
                        'name': ss.name,
                        'source': ss.source,
                        'rna': {
                            'name': ss.rna.name,
                            'sequence': ss.rna.sequence,
                            'source': ss.rna.source,
                            '_id': ss.rna._id
                        }    
                    }

                    helices_descr = []
                    for helix in ss.helices:
                        helix_desc = {
                            'name': helix['name'],
                            'location': {'ends': helix['location']} if version == 1 else helix['location']
                        }
                        if helix.has_key('interactions'):
                            interactions_descr = []
                            for interaction in helix['interactions']:
                                interactions_descr.append({
                                    'orientation': interaction['orientation'],
                                    'edge1': interaction['edge1'],
                                    'edge2': interaction['edge2'],
                                    'location': {'ends': interaction['location']} if version == 1 else interaction['location']
                                })
                            helix_desc['interactions'] = interactions_descr   
                        
                        helices_descr.append(helix_desc)

                    _2D_descr['helices'] = helices_descr

                    single_strands_descr = []
                    for single_strand in ss.single_strands:
                        single_strands_descr.append({
                            'name': single_strand['name'],
                            'location': {'ends': single_strand['location']} if version == 1 else single_strand['location']
                        })

                    _2D_descr['singleStrands'] = single_strands_descr

                    tertiary_interactions_descr = []
                    for tertiary_interaction in ss.tertiary_interactions:
                        tertiary_interactions_descr.append({
                            'orientation': tertiary_interaction['orientation'],
                            'edge1': tertiary_interaction['edge1'],
                            'edge2': tertiary_interaction['edge2'],
                            'location': {'ends': tertiary_interaction['location']} if version == 1 else tertiary_interaction['location']
                        })

                    _2D_descr['tertiaryInteractions'] = tertiary_interactions_descr

                    junctions_descr = []

                    for junction in ss.junctions:
                        junctions_descr.append({
                            'description': junction['description'],
                            'location': junction['location']
                        })


                    _2D_descr['junctions'] = junctions_descr 

                    _3D_descr = {
                        '_id': ts._id,
                        'name': ts.name,
                        'source': ts.source,
                        'rna': {
                            'name': ts.rna.name,
                            'sequence': ts.rna.sequence,
                            'source': ts.rna.source,
                            '_id': ts.rna._id
                        }     
                    }

                    residues_descr = {}
                    keys=[]
                    for k in ts.residues:
                        keys.append(k)

                    keys.sort() #the absolute position are sorted

                    for key in keys:
                        atoms = ts.residues[key]['atoms']
                        
                        atoms_descr = []

                        for atom in atoms:
                            atoms_descr.append({
                                'name': atom['name'],
                                'coords': atom['coords']
                            })
                        residues_descr[str(key)] = {
                            'atoms': atoms_descr
                        }

                    _3D_descr['residues'] = residues_descr


                    result.append({"2D": _2D_descr, "3D": _3D_descr})

                return Response(ujson.dumps(result), mimetype='application/json')

    #@app.route('/api/compute/2dplot', methods=['GET', 'POST'])
    def plot_2d():
        if request.method == 'POST':
            if 'data' in request.form:
                vienna_data = request.form['data']
        else:
            if 'data' in request.args:
                vienna_data = request.args.get('data', None)

        log = {
            'path': request.path,
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }

        logs_db['webservices'].insert(log)

        rnas, base_pairs = parse_vienna(vienna_data)
        rnaplot = Rnaplot()
        plot =  rnaplot.plot(base_pairs[0], rnas[0])
        coords = []
        for (index, row) in plot.iterrows():
            coords.append([row['x'], row['y']])
        return Response(ujson.dumps(coords), mimetype='application/json')

    #@app.route('/api/pdb', methods=['GET', 'POST'])
    def pdb():
        result = None
        collection = None
        query = None
        count = False
        id = None
        if request.method == 'POST':
            if 'coll' in request.form:
                collection = request.form['coll']
            if 'query' in request.form:
                query = request.form['query']
            if 'id' in request.form:
                id = request.form['id']
            if 'count' in request.form:
                count = True
        else:
            if 'coll' in request.args:
                collection = request.args.get('coll', None)
            if 'query' in request.args:
                query = request.args.get('query', None)
            if 'id' in request.args:
                id = request.args.get('id', None)
            if 'count' in request.args:
                count = True
        
        log = {
            'path': request.path,
            'collection': collection,
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }

        logs_db['webservices'].insert(log)

        db = mongodb['PDB']

        if collection and query:
            import simplejson
            result = list(db[collection].find(simplejson.loads(str(query))))
            if count:
                result = '{"count":%i}'%len(result)
        elif collection and id:
            result = db[collection].find_one({'_id':id})
        elif collection:
            result = list(db[collection].find())
            if count:
                result = '{"count":%i}'%len(result)
        return Response(ujson.dumps(result), mimetype='application/json')

    #@app.route('/api/rna3dhub', methods=['GET', 'POST'])
    def rna3dhub():
        result = None
        collection = None
        query = None
        count = False
        id = None
        if request.method == 'POST':
            if 'coll' in request.form:
                collection = request.form['coll']
            if 'query' in request.form:
                query = request.form['query']
            if 'id' in request.form:
                id = request.form['id']
            if 'count' in request.form:
                count = True
        else:
            if 'coll' in request.args:
                collection = request.args.get('coll', None)
            if 'query' in request.args:
                query = request.args.get('query', None)
            if 'id' in request.args:
                id = request.args.get('id', None)
            if 'count' in request.args:
                count = True
        
        log = {
            '_id': str(ObjectId()),
            'path': request.path,
            'collection': collection,
            'ip': request.remote_addr,
            'method': request.method,
            'date': str(datetime.datetime.now())
        }

        logs_db['webservices'].insert(log)

        db = mongodb['RNA3DHub']

        if collection and query:
            import simplejson
            result = list(db[collection].find(simplejson.loads(str(query))))
            if count:
                result = '{"count":%i}'%len(result)
        elif collection and id:
            result = db[collection].find_one({'_id':id})
        elif collection:
            result = list(db[collection].find())
            if count:
                result = '{"count":%i}'%len(result)
        return Response(ujson.dumps(result), mimetype='application/json')

class WebSocketHandler(tornado.websocket.WebSocketHandler):
    clients = {}

    def open(self, *args):
        self.id = uuid.uuid4()
        self.clients[self.id] = {'id':self.id}
        self.clients[self.id]['current_coverages_+'] = []
        self.clients[self.id]['current_coverages_-'] = []
        self.clients[self.id]['current_cufflinks_models'] = []
        self.clients[self.id]['current_gmorse_models'] = []
        print "New client connected %s"%self.request.remote_ip

    def on_message(self, message):
        import json
        message = json.loads(message)

        if message['header'] == 'get genome':
            self.clients[self.id]['species'] = message['species']
            answer = {'header': 'got genome'}
            answer['center'] = message['center']
            answer['annotation_id'] = message['annotation_id']
            self.write_message(answer, binary = False)

        elif message['header'] == 'get motifs':
            genome = mongodb[self.clients[self.id]['species']]['genomes'].find_one({'_id': message['genome_id']})
            genomicStart = message['genomicStart']
            genomicEnd = message['genomicEnd']
            answer = {'header': 'got genomic sequence to display'}
            answer['sequence'] = genome['sequence'][genomicStart-1:genomicEnd]
            answer['genomicStart'] = genomicStart
            answer['genomicEnd'] = genomicEnd
            self.write_message(answer, binary = False)
            for motif in mongodb[self.clients[self.id]['species']]['motifs'].find({'genome': message['genome_id']+"@genomes"}):
                if motif['genomicPositions'][0] >= genomicStart and motif['genomicPositions'][0] <= genomicEnd or motif['genomicPositions'][1] >= genomicStart and motif['genomicPositions'][1] <= genomicEnd:
                    answer = {'header': 'got new motif to display'}
                    answer['motif'] = motif
                    self.write_message(answer, binary = False)

        elif message['header'] == 'display coverage +':
            genomic_range = message['genomic_range']
            genomic_step = message['genomic_step']
            coverage_name = message['name']
            chromosome_name = mongodb[self.clients[self.id]['species']]['genomes'].find_one({'_id':message['genome_id']}, {'name':1})['name']
            self.clients[self.id]['current_coverages_+'].append(coverage_name)
            if self.clients[self.id]['species'] == 'Candida_glabrata_CBS_138':
                answer = {
                    'header': 'got coverage + to display',
                    'name': coverage_name,
                    'genomicStart': genomic_range[0],
                    'genomicEnd': genomic_range[1]
                }
                coverages = []
                for i in xrange(0,121):
                    coverage = {}
                    coverage['rank'] = i
                    coverage['genomicStart'] = int(genomic_range[0]+i*genomic_step*10)
                    coverage['genomicEnd'] = int(genomic_range[0]+(i+1)*genomic_step*10-1)
                    coverage['count'] = sum(counts["%s/%s_+.pickle"%(coverage_name, chromosome_name)][coverage['genomicStart']-1:coverage['genomicEnd']])
                    coverages.append(coverage)
                answer['coverages'] = coverages
                self.write_message(answer, binary = False)

        elif message['header'] == 'hide coverage +':
            self.clients[self.id]['current_coverages_+'].remove(message['name'])

        elif message['header'] == 'display coverage -':
            genomic_range = message['genomic_range']
            genomic_step =  message['genomic_step']
            coverage_name = message['name']
            chromosome_name = mongodb[self.clients[self.id]['species']]['genomes'].find_one({'_id':message['genome_id']}, {'name':1})['name']
            self.clients[self.id]['current_coverages_-'].append(coverage_name)
            if self.clients[self.id]['species'] == 'Candida_glabrata_CBS_138':
                answer = {
                    'header': 'got coverage - to display',
                    'name': coverage_name,
                    'genomicStart': genomic_range[0],
                    'genomicEnd': genomic_range[1]
                }
                coverages = []
                for i in xrange(0,121):
                    coverage = {}
                    coverage['rank'] = i
                    coverage['genomicStart'] = int(genomic_range[0]+i*genomic_step*10)
                    coverage['genomicEnd'] = int(genomic_range[0]+(i+1)*genomic_step*10-1)
                    coverage['count'] = sum(counts["%s/%s_-.pickle"%(coverage_name, chromosome_name)][coverage['genomicStart']-1:coverage['genomicEnd']])
                    coverages.append(coverage)
                answer['coverages'] = coverages
                self.write_message(answer, binary = False)
        
        elif message['header'] == 'hide coverage -':
            self.clients[self.id]['current_coverages_-'].remove(message['name'])

        elif message['header'] == 'display cufflinks models':
            genomic_range = message['genomic_range']
            sample_name = message['name']
            self.clients[self.id]['current_cufflinks_models'].append(sample_name)
            chromosome_name = mongodb[self.clients[self.id]['species']]['genomes'].find_one({'_id':message['genome_id']}, {'name':1})['name']
            if self.clients[self.id]['species'] == 'Candida_glabrata_CBS_138':
                for sample_name in self.clients[self.id].get('current_cufflinks_models', []):
                    for gene_model in mongodb[self.clients[self.id]['species']]['gene_models'].find({'genome': message['genome_id']+"@genomes", 'source':'tool:cufflinks:NA', 'sample':sample_name}):
                        if int(gene_model['genomicPositions'][0]) >= genomic_range[0] and int(gene_model['genomicPositions'][0]) <= genomic_range[1] or int(gene_model['genomicPositions'][1]) >= genomic_range[0] and int(gene_model['genomicPositions'][1]) <= genomic_range[1]:
                            answer = {'header': 'got new cufflinks model to display'}
                            answer['gene model'] = gene_model
                            self.write_message(answer, binary = False)
        
        elif message['header'] == 'hide cufflinks models':
            self.clients[self.id]['current_cufflinks_models'].remove(message['name'])

        elif message['header'] == 'display gmorse models':
            genomic_range = message['genomic_range']
            sample_name = message['name']
            self.clients[self.id]['current_gmorse_models'].append(sample_name)
            chromosome_name = mongodb[self.clients[self.id]['species']]['genomes'].find_one({'_id':message['genome_id']}, {'name':1})['name']
            if self.clients[self.id]['species'] == 'Candida_glabrata_CBS_138':
                for sample_name in self.clients[self.id].get('current_gmorse_models', []):
                    for gene_model in mongodb[self.clients[self.id]['species']]['gene_models'].find({'genome': message['genome_id']+"@genomes", 'source':'tool:gmorse:NA', 'sample':sample_name}):
                        if int(gene_model['genomicPositions'][0]) >= genomic_range[0] and int(gene_model['genomicPositions'][0]) <= genomic_range[1] or int(gene_model['genomicPositions'][1]) >= genomic_range[0] and int(gene_model['genomicPositions'][1]) <= genomic_range[1]:
                            answer = {'header': 'got new gmorse model to display'}
                            answer['gene model'] = gene_model
                            self.write_message(answer, binary = False)
        
        elif message['header'] == 'hide gmorse models':
            self.clients[self.id]['current_gmorse_models'].remove(message['name'])

        elif message['header'] == 'genome browser dragged':
            current_ids = message['current_ids']
            genomic_range = message['genomic_range']
            genomic_step = message['genomic_step']
            chromosome_name = mongodb[self.clients[self.id]['species']]['genomes'].find_one({'_id':message['genome_id']}, {'name':1})['name']
            answer = {'header': 'got new genomic range'}
            self.write_message(answer, binary = False) #just to send a mesage to update the start and end fields in browser.html
            if message['source'] == 'annotationsBrowser':
                for annotation in mongodb[self.clients[self.id]['species']]['annotations'].find({'genome': message['genome_id']+"@genomes"}):
                    if not annotation['_id'] in current_ids and annotation['genomicPositions'][0] >= genomic_range[0] and annotation['genomicPositions'][0] <= genomic_range[1] or annotation['genomicPositions'][1] >= genomic_range[0] and annotation['genomicPositions'][1] <= genomic_range[1]:
                        answer = {'header': 'got new annotation to display'}
                        answer['annotation'] = annotation
                        self.write_message(answer, binary = False)
            elif message['source'] == 'ncrnasBrowser':
                for ncRNA in mongodb[self.clients[self.id]['species']]['ncRNAs'].find({'genome': message['genome_id']+"@genomes"}):
                    if not ncRNA['_id'] in current_ids and ncRNA['genomicPositions'][0] >= genomic_range[0] and ncRNA['genomicPositions'][0] <= genomic_range[1] or ncRNA['genomicPositions'][1] >= genomic_range[0] and ncRNA['genomicPositions'][1] <= genomic_range[1]:
                        answer = {'header': 'got new ncrna to display'}
                        answer['ncrna'] = ncRNA
                        self.write_message(answer, binary = False)
            elif message['source'] == 'cufflinksBrowser':
                for sample_name in self.clients[self.id].get('current_cufflinks_models', []):
                    for gene_model in mongodb[self.clients[self.id]['species']]['gene_models'].find({'genome': message['genome_id']+"@genomes", 'source':'tool:cufflinks:NA', 'sample':sample_name}):
                        if not gene_model['_id'] in current_ids and int(gene_model['genomicPositions'][0]) >= genomic_range[0] and int(gene_model['genomicPositions'][0]) <= genomic_range[1] or int(gene_model['genomicPositions'][1]) >= genomic_range[0] and int(gene_model['genomicPositions'][1]) <= genomic_range[1]:
                            answer = {'header': 'got new cufflinks model to display'}
                            answer['gene model'] = gene_model
                            self.write_message(answer, binary = False)
            elif message['source'] == 'gmorseBrowser':
                for sample_name in self.clients[self.id].get('current_gmorse_models', []):
                    for gene_model in mongodb[self.clients[self.id]['species']]['gene_models'].find({'genome': message['genome_id']+"@genomes", 'source':'tool:gmorse:NA', 'sample':sample_name}):
                        if not gene_model['_id'] in current_ids and int(gene_model['genomicPositions'][0]) >= genomic_range[0] and int(gene_model['genomicPositions'][0]) <= genomic_range[1] or int(gene_model['genomicPositions'][1]) >= genomic_range[0] and int(gene_model['genomicPositions'][1]) <= genomic_range[1]:
                            answer = {'header': 'got new gmorse model to display'}
                            answer['gene model'] = gene_model
                            self.write_message(answer, binary = False)
            elif message['source'] == 'plusReadsBrowser':
                if self.clients[self.id]['species'] == 'Candida_glabrata_CBS_138':
                    for current_coverage in self.clients[self.id].get('current_coverages_+', []):
                        answer = {
                            'header': 'got coverage + to display',
                            'name': current_coverage,
                            'genomicStart': genomic_range[0],
                            'genomicEnd': genomic_range[1]
                        }
                        coverages = []
                        for i in xrange(0,121):
                            coverage = {}
                            coverage['rank'] = i
                            coverage['genomicStart'] = int(genomic_range[0]+i*genomic_step*10)
                            coverage['genomicEnd'] = int(genomic_range[0]+(i+1)*genomic_step*10-1)
                            coverage['count'] = sum(counts["%s/%s_+.pickle"%(current_coverage, chromosome_name)][coverage['genomicStart']-1:coverage['genomicEnd']])
                            coverages.append(coverage)
                        answer['coverages'] = coverages
                        self.write_message(answer, binary = False)
            elif message['source'] == 'minusReadsBrowser':
                if self.clients[self.id]['species'] == 'Candida_glabrata_CBS_138':
                    for current_coverage in self.clients[self.id].get('current_coverages_-', []):
                        answer = {
                            'header': 'got coverage - to display',
                            'name': current_coverage,
                            'genomicStart': genomic_range[0],
                            'genomicEnd': genomic_range[1]
                        }
                        coverages = []
                        for i in xrange(0,121):
                            coverage = {}
                            coverage['rank'] = i
                            coverage['genomicStart'] = int(genomic_range[0]+i*genomic_step*10)
                            coverage['genomicEnd'] = int(genomic_range[0]+(i+1)*genomic_step*10-1)
                            coverage['count'] = sum(counts["%s/%s_-.pickle"%(current_coverage, chromosome_name)][coverage['genomicStart']-1:coverage['genomicEnd']])
                            coverages.append(coverage)
                        answer['coverages'] = coverages
                        self.write_message(answer, binary = False)

            elif message['source'] == 'motifsBrowser':
                for motif in mongodb[self.clients[self.id]['species']]['motifs'].find({'genome': message['genome_id']+"@genomes"}):
                    if not motif['_id'] in current_ids and motif['genomicPositions'][0] >= genomic_range[0] and motif['genomicPositions'][0] <= genomic_range[1] or motif['genomicPositions'][1] >= genomic_range[0] and motif['genomicPositions'][1] <= genomic_range[1]:
                        answer = {'header': 'got new motif to display'}
                        answer['motif'] = motif
                        self.write_message(answer, binary = False)

        elif message['header'] == 'get genome list':
            self.clients[self.id]['species'] = message['species']
            genome_list = mongodb[message['species']]['genomes'].find({}, {"_id":1, "name":1})
            self.clients[self.id]['current_coverages_+'] = []
            self.clients[self.id]['current_coverages_-'] = []
            
            answer = {'header': 'got genome list'}
            answer['genomes'] = genome_list
            first_genome_id = genome_list[0]['_id']
            answer['start'] = 1
            answer['end'] = len(mongodb[message['species']]['genomes'].find_one({"_id": first_genome_id}, {"sequence":1})['sequence'])

            self.write_message(dumps(answer), binary = False)
        
        elif message['header'] == 'get genome attributes':
            answer = {'header': 'got genome attributes'}
            answer['start'] = 1
            answer['end'] = len(mongodb[message['species']]['genomes'].find_one({"_id": message['genome_id']}, {"sequence":1})['sequence'])
            self.write_message(answer, binary = False)

        elif message['header'] == 'get 2d':
            answer = {'header': 'got 2d'}
            answer['2d'], orthologs, aligned_sequences = get_ortholog_2d(self.clients[self.id]['booquet']['alignment'], message['ortholog'])
            self.write_message(answer, binary = False)

        elif message['header'] == 'init booquet':
            junction_diameter = 20
            aligned_sequences = None
            bn = None
            ss_json = None
            consensus_ss_json = None
            orthologs = None
            ncRNA_name = None
            structural_alignment = None
            if message['alignment'].startswith('RF'):
                ncRNA_name = message['alignment']
                stockholm_content = rfam.get_entry(rfam_id = ncRNA_name, format='stockholm')
                (aligned_molecules, organisms, consensus2D) = parse_stockholm(stockholm_content)
                structural_alignment = to_clustalw(consensus2D, aligned_molecules)
            else: #client own alignment
                structural_alignment = message['alignment']
            ss_json, orthologs, aligned_sequences = get_ortholog_2d(structural_alignment, None)
            consensus_ss_json, bn = get_consensus_2d(structural_alignment, junction_diameter)
            self.clients[self.id]['booquet'] = {'alignment':structural_alignment}
            answer = {'header': 'init booquet'}
            answer['secondary_structure'] = ss_json
            answer['consensus_secondary_structure'] = consensus_ss_json
            answer['junction_diameter'] = junction_diameter
            answer['aligned_sequences'] = aligned_sequences
            answer['orthologs'] = orthologs
            answer['bn'] = bn

            self.write_message(answer, binary = False)

        
    def on_close(self):
        self.clients.pop(self.id, None)
        print "Client disconnected %s"%self.request.remote_ip

class Application(tornado.web.Application):
    def __init__(self):

        handlers = [
            (r'/', IndexHandler),
            (r'/help', HelpHandler),
            (r'/browser', GenomeBrowserHandler),
            (r'/data', GetDataHandler),
            (r'/families', FamiliesHandler),
            (r'/species', SpeciesHandler),
            (r'/synthenies', SyntheniesHandler),
            (r'/consensus_structure', ConsensusStructureHandler),
            (r'/s2s', S2SHandler),
            (r'/websocket', WebSocketHandler),
            (r'/booquet', BooquetHandler),
            (r'/booquet/sample', BooquetSampleHandler)
        ]

        settings = {
            'template_path': 'templates',
            'static_path': 'static',
            'debug': True
        }

        tornado.web.Application.__init__(self, handlers, **settings)

if __name__ == '__main__':
    webserver_host = "localhost"
    webserver_port = 8080
    mongodb_host = "localhost"
    mongodb_port = 27017

    if "-wh" in sys.argv:
        webserver_host = sys.argv[sys.argv.index("-wh")+1]
    if "-wp" in sys.argv:
        webserver_port = int(sys.argv[sys.argv.index("-wp")+1])
    if "-mh" in sys.argv:
        mongodb_host = sys.argv[sys.argv.index("-mh")+1]
    if "-mp" in sys.argv:
        mongodb_port = int(sys.argv[sys.argv.index("-mp")+1])
    if "-conf" in sys.argv:
        with open(sys.argv[sys.argv.index("-conf")+1]) as h:
            json_data = h.read()
        enabled_algorithms = ujson.loads(json_data)["rest_server"]["enable-algorithms"]

    try :
        mongodb = MongoClient(mongodb_host, mongodb_port)
        logs_db = mongodb['logs']
    except Exception, e:
        print 'Cannot connect any Mongodb instance hosted at %s:%i'%(mongodb_host, mongodb_port)
        print 'Usage: ./server.py [-wh webserver_host (default: %s)] [-wp webserver_port (default: %i)] [-mh mongodb_host (default: %s)] [-mp mongodb_port (default: %i)] [-conf configuration_file] '%(webserver_host, webserver_port, mongodb_host, mongodb_port)
        sys.exit(-1)

    # sample_directories = [f for f in os.listdir('static/data/RNA-seq/Candida_glabrata_CBS_138/')]
    # for sample_directory in sample_directories:
    #     if sample_directory in ['FJT-05_5_S1', 'FJT-10_10_S2']:
    #         if os.path.isdir('static/data/RNA-seq/Candida_glabrata_CBS_138/%s'%sample_directory):
    #             filenames = [f for f in os.listdir('static/data/RNA-seq/Candida_glabrata_CBS_138/%s'%sample_directory) if re.match(r'.+.pickle', f)]
    #             for filename in filenames:
    #                 print 'load %s'%('static/data/RNA-seq/Candida_glabrata_CBS_138/%s/%s'%(sample_directory, filename))
    #                 with open('static/data/RNA-seq/Candida_glabrata_CBS_138/%s/%s'%(sample_directory, filename), 'rb') as h:
    #                     counts["%s/%s"%(sample_directory,filename)] = pickle.load(h)

    # print "Pickle files for RNA-seq reads parsed efficiently..."

    tornado.options.parse_command_line()
    app = Application()
    server = tornado.httpserver.HTTPServer(app)
    server.listen(webserver_port)
    main_loop = tornado.ioloop.IOLoop.instance()
    main_loop.start()
