#!/usr/bin/env python

from flask import Flask, request, Response, render_template, redirect, url_for, jsonify, abort
from pyrna.computations import Rnafold, Contrafold, Rnaplot, Rnaview, Mlocarna, Rnasubopt, RnaAlifold
from pyrna.db import PDB
from pyrna import parsers
from pyrna.parsers import parse_vienna, parse_fasta, base_pairs_to_secondary_structure, parse_pdb, to_clustalw
from pyrna.features import RNA
import ujson, sys, datetime, shutil, re, os, random, string, shlex
from pymongo import MongoClient
from bson.objectid import ObjectId
from pandas import DataFrame, isnull
from werkzeug import secure_filename
from subprocess import Popen

app = Flask(__name__)

mongodb = None
logs_db = None
enabled_algorithms = []

################# TEAM WEBSITE #########################

@app.route("/")
def index():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)

    return render_template('index.html')

@app.route("/team")
def team():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('team.html')

@app.route("/methods")
def methods():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('methods.html')

@app.route("/publications")
def publications():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('publications.html')

@app.route("/projects")
def projects():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('projects.html')

@app.route("/tools")
def tools():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('tools.html')

@app.route("/contacts")
def contacts():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('contacts.html')

@app.route("/teaching")
def teaching():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('teaching.html')

################# WEBSERVICES #########################

@app.route("/webservices")
@app.route('/api/')
def webservices():
    return render_template('webservices.html')

#here starts the low-level webservices: webservices to avoid to install RNA algorithms on the client computer to be able to use PyRNA

@app.route('/api/get_key', methods=['GET'])
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

@app.route('/api/computations/rnafold', methods=['POST'])
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

@app.route('/api/computations/rnaplot', methods=['POST'])
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

@app.route('/api/computations/contrafold', methods=['POST'])
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

@app.route('/api/computations/rnaview', methods=['POST'])
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

@app.route('/api/compute/2d', methods=['GET', 'POST'])
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
                        'crown': junction['crown']
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
                for k in ts.residues.keys():
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

@app.route('/api/compute/2dplot', methods=['GET', 'POST'])

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

@app.route('/api/pdb', methods=['GET', 'POST'])
@app.route('/api/ark', methods=['GET', 'POST'])
@app.route('/api/charn', methods=['GET', 'POST'])
def ark():
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

    if request.path.endswith('pdb') or request.path.endswith('ark'):
        db = mongodb['PDB']
    elif request.path.endswith('charn'):
        db = mongodb['charndb']

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

################# CHARNDB #########################

@app.route("/charndb")
def charndb():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('charndb.html')

@app.route("/charndb/ncRNAs", methods=['GET', 'POST'])
def ncRNAs():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    species_names = []
    if request.method == 'POST':
        species_names = request.form['species']
    else:
        species_names = request.args.get('species', [])       
    print species_names
    return render_template('/charndb/ncRNAs.html', species_names = species_names)

@app.route("/charndb/species")
def species():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('/charndb/species.html')

@app.route("/charndb/orthoparalogs")
def orthoparalogs():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    ncRNA_name = request.args.get('name', None)
    return render_template('/charndb/orthoparalogs.html', ncRNA_name = ncRNA_name)

@app.route("/charndb/downloads")
def downloads():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('/charndb/downloads.html')

################# LETITSNO #########################

ALLOWED_EXTENSIONS = set(['txt', 'fasta', 'fna', 'embl', 'gb'])

@app.route("/letitsno")
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

@app.route('/letitsno/upload', methods=['POST'])
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

            info = open(os.getenv("HOME")+"/tmp/jobs_let_it_sno.task_on_%s/info"%secret_key,'a+')
            info.write("process ID: %i\n"%p.pid)
            info.close()

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

@app.route('/letitsno/results')
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

@app.route('/letitsno/check')
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
        h = open(os.getenv("HOME")+'/tmp/jobs_let_it_sno.task_on_%s/info'%request.args['id'], 'r')
        info_content = h.read()
        h.close()

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

@app.route('/letitsno/documentation')
def documentation():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('letitsno/documentation.html')

@app.route('/letitsno/stats')
def stats():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('letitsno/stats.html')

@app.route('/letitsno/example')
def example():
    log = {
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }
    logs_db['websites'].insert(log)
    return render_template('letitsno/example.html')

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
        h = open(sys.argv[sys.argv.index("-conf")+1], 'r')
        json_data = h.read()
        h.close()
        enabled_algorithms = ujson.loads(json_data)["rest_server"]["enable-algorithms"]

    try :
        mongodb = MongoClient(mongodb_host, mongodb_port)
        logs_db = mongodb['logs']
    except Exception, e:
        print 'Cannot connect any Mongodb instance hosted at %s:%i'%(mongodb_host, mongodb_port)
        print 'Usage: ./server.py [-wh webserver_host (default: %s)] [-wp webserver_port (default: %i)] [-mh mongodb_host (default: %s)] [-mp mongodb_port (default: %i)] [-conf configuration_file] '%(webserver_host, webserver_port, mongodb_host, mongodb_port)
        sys.exit(-1)

    app.debug = True

    #app.run(host = webserver_host, port = webserver_port)
    from tornado.wsgi import WSGIContainer
    from tornado.httpserver import HTTPServer
    from tornado.ioloop import IOLoop

    http_server = HTTPServer(WSGIContainer(app))
    http_server.listen(webserver_port)
    IOLoop.instance().start()