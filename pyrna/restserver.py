#!/usr/bin/env python

from flask import Flask, request, Response ,render_template, abort
from pyrna.computations import Rnafold, Contrafold, Rnaplot, Rnaview, Mlocarna, Rnasubopt
from pyrna.db import PDB
from pyrna.parsers import parse_vienna, parse_fasta, base_pairs_to_secondary_structure, parse_pdb, to_clustalw
from pyrna.features import RNA
import ujson, sys, datetime, random, string, json
from pymongo import MongoClient
from bson.objectid import ObjectId

app = Flask(__name__)

mongodb = None
logs_db = None
enabled_algorithms = []

################# WEBSERVICES #########################

@app.route("/")
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
        '_id': str(ObjectId()),
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
        '_id': str(ObjectId()),  
        'path': request.path,
        'ip': request.remote_addr,
        'method': request.method,
        'date': str(datetime.datetime.now())
    }

    logs_db['webservices'].insert(log)

    rnas, base_pairs = parse_vienna(vienna_data)
    rnaplot = Rnaplot()
    plot =  rnaplot.plot(base_pairs[0], rna[0])
    coords = []
    for (index, row) in plot.iterrows():
        coords.append([row['x'], row['y']])
    return Response(ujson.dumps(coords), mimetype='application/json')

@app.route('/api/pdb', methods=['GET', 'POST'])
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
        '_id': str(ObjectId()),
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

@app.route('/api/rna3dhub', methods=['GET', 'POST'])
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
        print 'Usage: ./restserver.py [-wh webserver_host (default: %s)] [-wp webserver_port (default: %i)] [-mh mongodb_host (default: %s)] [-mp mongodb_port (default: %i)] [-conf configuration_file] '%(webserver_host, webserver_port, mongodb_host, mongodb_port)
        sys.exit(-1)

    app.debug = True

    #app.run(host = webserver_host, port = webserver_port)
    from tornado.wsgi import WSGIContainer
    from tornado.httpserver import HTTPServer
    from tornado.ioloop import IOLoop

    http_server = HTTPServer(WSGIContainer(app))
    http_server.listen(webserver_port)
    IOLoop.instance().start()