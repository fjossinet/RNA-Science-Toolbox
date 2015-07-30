#!/usr/bin/env python

import ujson, sys, datetime, os, random, string, json

from pyrna.features import RNA
from pyrna.computations import Rnafold, Contrafold, Rnaplot, Rnaview, Mlocarna, Rnasubopt, RnaAlifold
from pyrna.db import PDB
from pyrna import parsers
from pyrna.parsers import parse_vienna, parse_fasta, base_pairs_to_secondary_structure, parse_pdb, to_clustalw
from pymongo import MongoClient
from bson.objectid import ObjectId
from subprocess import Popen

import tornado.httpserver
import tornado.ioloop
import tornado.options
import tornado.web
import tornado.websocket
from tornado.escape import json_encode, native_str

static_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../website')
pages_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../website/pages')

app = None
websockets = []
external_tools_2_websockets = {}
mongodb = None
logs_db = None
webserver_db = None
enabled_algorithms = ['rnafold', 'rnaplot', 'contrafold', 'rnaview']

def is_registered_user(secret_key):
    return logs_db['user_keys'].find_one({'key': secret_key}) != None

#to write messsages to websockets outside of the WebSocketHandler class
def wsSend(message):
    for ws in websockets:
        if not ws.ws_connection.stream.socket:
            print "Web socket does not exist anymore!!!"
            websockets.remove(ws)
        else:
            ws.write_message(message)

class IndexHandler(tornado.web.RequestHandler):
    def get(self):
        if not os.path.exists(static_dir):
            self.write("RNA WebServices running...")
        else:
            self.render('index.html', hostname = hostname)

class ServerHandler(tornado.web.RequestHandler):
    def get(self):
        self.render('server.html', hostname = hostname)

class PyrnaDocHandler(tornado.web.RequestHandler):
    def get(self):
        self.render('pyrna.html', hostname = hostname)

class WebservicesDocHandler(tornado.web.RequestHandler):
    def get(self):
        self.render('webservices.html', hostname = hostname)

class RnapediaHandler(tornado.web.RequestHandler):
    def get(self):
        self.render('rnapedia.html', hostname = hostname)

class AccountHandler(tornado.web.RequestHandler):
    def get(self):
        if not self.get_secure_cookie("login"):
            self.redirect("/login")
            return
        user = webserver_db['users'].find_one({'login': self.get_secure_cookie("login")})
        self.render('account.html', hostname = hostname, user = user, tool_name = self.get_argument('tool_name', default = None), tool_id = self.get_argument('tool_id', default = None), tools_online = external_tools_2_websockets.keys())

class RegisterUserHandler(tornado.web.RequestHandler):
    def post(self):
        webserver_db['users'].insert({
            'login': self.get_argument("login"),
            'password': self.get_argument("password"),
            'external tools linked': []
        })
        self.set_secure_cookie("login", self.get_argument("login"))
        self.redirect("/account")

class LoginHandler(tornado.web.RequestHandler):
    def get(self):
        self.render('login.html', hostname = hostname)

    def post(self):
        user = webserver_db['users'].find_one({'login': self.get_argument("login"), 'password': self.get_argument("password")})
        if not user:
            self.redirect("/login")
        else:
            self.set_secure_cookie("login", user['login'])
            self.redirect("/account")

class LogoutHandler(AccountHandler):
    def get(self):
        user = webserver_db['users'].find_one({'login': self.get_secure_cookie("login")})
        for tool in user['external tools linked']:
            print external_tools_2_websockets.has_key(tool['id'])
            print external_tools_2_websockets
            print tool['id']
            if external_tools_2_websockets.has_key(tool['id']):
                message = {}
                message['header'] = "super"
                external_tools_2_websockets[tool['id']].write_message(message)
        user['last_connection'] = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')
        webserver_db['users'].save(user)
        self.clear_cookie("login")
        self.redirect("/login")

class LinkToolHandler(tornado.web.RequestHandler):
    def post(self):
        if len(self.get_argument("tool_name")) !=0 and len(self.get_argument("tool_id")) != 0:
            user = webserver_db['users'].find_one({'login': self.get_secure_cookie("login")})
            tools = user.get('external tools linked', [])
            tools.append({
                'name': self.get_argument("tool_name"),
                'id': self.get_argument("tool_id")
            })
            user['external tools linked'] = tools
            webserver_db['users'].save(user)
        self.redirect("/account")

class UnLinkToolHandler(tornado.web.RequestHandler):
    def post(self):
        user = webserver_db['users'].find_one({'login': self.get_secure_cookie("login")})
        tools = user.get('external tools linked', [])
        tools = [tool for tool in tools if tool.get('id') != self.get_argument("tool_id")]
        user['external tools linked'] = tools
        webserver_db['users'].save(user)
        self.redirect("/account")

##########################################################
# Here starts the low-level webservices:
# webservices to avoid to install RNA algorithms on the
# client computer to be able to use PyRNA.
##########################################################

#webservice to generate and register an api key (needed to call webservices from python scripts or ipython sessions)
class APIKeyHandler(tornado.web.RequestHandler):
    def get(self):
        secret_key = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(10))
        remote_ip = self.request.remote_ip
        record = logs_db['user_keys'].find_one({'ip': remote_ip})
        if record:
            record['key'] = secret_key
            logs_db['user_keys'].save(record)
        else:
            logs_db['user_keys'].insert({
                '_id': str(ObjectId()),
                'ip': remote_ip,
                'date': datetime.datetime.now(),
                'key': secret_key
                })
        self.write(secret_key)

#webservice to run RNAfold
class RNAfoldHandler(tornado.web.RequestHandler):
    def post(self):
        log = {
            '_id': str(ObjectId()),
            'path': self.request.uri,
            'tool': 'rnafold',
            'ip': self.request.remote_ip,
            'method': self.request.method,
            'date': datetime.datetime.now()
        }

        logs_db['webservices'].insert(log)

        api_key = self.get_argument('api_key', default = None)
        function = self.get_argument('function', default = None)

        if not is_registered_user(api_key) or not 'rnafold' in enabled_algorithms:
            self.send_error(status_code=401)
        else:
            rnafold = Rnafold()
            if function == 'fold':
                name = self.get_argument('name', default = None)
                sequence = self.get_argument('sequence', default = None)
                constraints = self.get_argument('constraints', default = None)
                bp_probabilities = self.get_argument('bp_probabilities', default = "False") == "True"
                self.write(rnafold.fold(RNA(name=name, sequence=sequence), bp_probabilities = bp_probabilities, constraints = constraints, raw_output = True))
            else:
                self.send_error(status_code=401)

#webservice to run RNAplot
class RNAplotHandler(tornado.web.RequestHandler):
    def post(self):
        log = {
            '_id': str(ObjectId()),
            'path': self.request.uri,
            'tool': 'rnaplot',
            'ip': self.request.remote_ip,
            'method': self.request.method,
            'date': datetime.datetime.now()
        }

        logs_db['webservices'].insert(log)

        secondary_structure = self.get_argument('secondary_structure', default = None)
        api_key = self.get_argument('api_key', default = None)

        if not is_registered_user(api_key) or not 'rnaplot' in enabled_algorithms:
            self.send_error(status_code=401)
        else:
            rnas, secondary_structures = parse_vienna(secondary_structure)
            self.write(Rnaplot().plot(secondary_structures[0], rnas[0], raw_output = True))

#webservice to run Contrafold
class ContrafoldHandler(tornado.web.RequestHandler):
    def post(self):
        log = {
            '_id': str(ObjectId()),
            'path': self.request.uri,
            'tool': 'contrafold',
            'ip': self.request.remote_ip,
            'method': self.request.method,
            'date': datetime.datetime.now()
        }

        logs_db['webservices'].insert(log)

        name = self.get_argument('name', default = None)
        sequence = self.get_argument('sequence', default = None)
        constraints = self.get_argument('constraints', default = None)
        api_key = self.get_argument('api_key', default = None)

        if not is_registered_user(api_key) or not 'contrafold' in enabled_algorithms:
            self.send_error(status_code=401)
        else:
            self.write(Contrafold().fold(RNA(name=name, sequence=sequence), raw_output = True))

#webservice to run RNAVIEW
class RnaviewHandler(tornado.web.RequestHandler):
    def post(self):
        log = {
            '_id': str(ObjectId()),
            'path': self.request.uri,
            'tool': 'rnaview',
            'ip': self.request.remote_ip,
            'method': self.request.method,
            'date': datetime.datetime.now()
        }

        logs_db['webservices'].insert(log)

        tertiary_structure = self.get_argument('3d', default = None)
        canonical_only = self.get_argument('canonical_only', default = None)
        api_key = self.get_argument('api_key', default = None)

        if not is_registered_user(api_key) or not 'rnaview' in enabled_algorithms:
            self.send_error(status_code=401)
        else:
            self.write(Rnaview().annotate(pdb_content = tertiary_structure, canonical_only = canonical_only, raw_output = True))

##########################################################
# Here starts the high-level webservices...
##########################################################

class Compute2dHandler(tornado.web.RequestHandler):

    def get(self):
        self.post(self)

    def post(self):
        data = self.get_argument('data', default = None)
        tool = self.get_argument('tool', default = None)
        version = self.get_argument('version', default = 1)
        pdbid =  self.get_argument('pdbid', default = None)
        output = None
        result = None

        log = {
            'path': self.request.uri,
            'tool': tool,
            'ip': self.request.remote_ip,
            'method': self.request.method,
            'date': datetime.datetime.now()
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
                    random_sample = int(self.get_argument('random_sample', default = 20))
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
                if tool == 'rnafold' or tool == 'contrafold':
                    self.write(json_encode(result[0]))
                else:
                    self.write(json_encode(result))
            elif len(rnas) >= 2: #structural alignment
                if tool == 'mlocarna':
                    aligned_molecules, consensus2D = Mlocarna().align(rnas)
                    self.write(to_clustalw(consensus2D, aligned_molecules))
        elif tool == 'rnalifold' and data and data.startswith('CLUSTAL'): #computation of consensus structure from sequence alignment
            self.write(RnaAlifold().align(data))
        elif tool == 'rnaview': #3D annotation
            rnaview = Rnaview()

            if output == 'rnaml':
                if pdbid:
                    self.write(rnaview.annotate(pdb_content = PDB().get_entry(pdbid), raw_output = True))
                elif data:
                    self.write(rnaview.annotate(pdb_content = data, raw_output = True))

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

                self.write(json_encode(result))

class Compute2dplotHandler(tornado.web.RequestHandler):

    def get(self):
        self.post(self)

    def post(self):
        vienna_data = self.get_argument('data', default = None)

        log = {
            'path': self.request.uri,
            'ip': self.request.remote_ip,
            'tool': 'rnaplot',
            'method': self.request.method,
            'date': datetime.datetime.now()
        }

        logs_db['webservices'].insert(log)

        rnas, base_pairs = parse_vienna(vienna_data)
        rnaplot = Rnaplot()
        plot =  rnaplot.plot(base_pairs[0], rnas[0])
        coords = []
        for (index, row) in plot.iterrows():
            coords.append([row['x'], row['y']])
        self.write(json_encode(coords))

class PDBHandler (tornado.web.RequestHandler):

    def get(self):
        self.post(self)

    def post(self):
        result = None
        collection = self.get_argument('coll', default = None)
        query = self.get_argument('query', default = None)
        id = self.get_argument('id', default = None)
        count = self.get_argument('count', default = False)

        log = {
            'path': self.request.uri,
            'collection': collection,
            'ip': self.request.remote_ip,
            'method': self.request.method,
            'date': datetime.datetime.now()
        }

        logs_db['webservices'].insert(log)

        db = mongodb['PDB']

        if collection and query:
            result = list(db[collection].find(ujson.loads(str(query))))
            if count:
                result = '{"count":%i}'%len(result)
        elif collection and id:
            result = db[collection].find_one({'_id':id})
        elif collection:
            result = list(db[collection].find())
            if count:
                result = '{"count":%i}'%len(result)
        self.write(json_encode(result))

class RNA3DHubHandler (tornado.web.RequestHandler):

    def get(self):
        self.post(self)

    def post(self):
        result = None
        collection = self.get_argument('coll', default = None)
        query = self.get_argument('query', default = None)
        id = self.get_argument('id', default = None)
        count = self.get_argument('count', default = False)

        log = {
            'path': self.request.uri,
            'collection': collection,
            'ip': self.request.remote_ip,
            'method': self.request.method,
            'date': datetime.datetime.now()
        }

        logs_db['webservices'].insert(log)

        db = mongodb['RNA3DHub']

        if collection and query:
            result = list(db[collection].find(ujson.loads(str(query))))
            if count:
                result = '{"count":%i}'%len(result)
        elif collection and id:
            result = db[collection].find_one({'_id':id})
        elif collection:
            result = list(db[collection].find())
            if count:
                result = '{"count":%i}'%len(result)

        self.write(json_encode(result))

class WebSocketHandler(tornado.websocket.WebSocketHandler):

    def open(self, *args):
        if self not in websockets:
            websockets.append(self)
        print "New client connected"

    def on_message(self, message):
        message = ujson.loads(message)
        if message['header'] == 'new external tool':
            for websocket in websockets:
                if self == websocket:
                    external_tools_2_websockets[message['id']] = websocket
        elif message['header'] == 'webservices usage':
            db = mongodb['logs']
            now = datetime.datetime.now()
            data = []
            for i in xrange(1, 24):
                counts = {'y': "-%ih"%
                i}
                time_range = {
                    "$gt": now - datetime.timedelta(hours = i),
                    "$lte": now - datetime.timedelta(hours = i-1)
                    }
                counts['RNAfold'] = db['webservices'].find({
                    "date": time_range,
                    "tool": 'rnafold'
                }).count()

                counts['RNAsubopt'] = db['webservices'].find({
                    "date": time_range,
                    "tool": 'rnasubopt'
                }).count()

                counts['RNAalifold'] = db['webservices'].find({
                    "date": time_range,
                    "tool": 'rnaalifold'
                }).count()

                counts['Contrafold'] = db['webservices'].find({
                    "date": time_range,
                    "tool": 'contrafold'
                }).count()

                counts['RNAplot'] = db['webservices'].find({
                    "date": time_range,
                    "tool": 'rnaplot'
                }).count()

                counts['Mlocarna'] = db['webservices'].find({
                    "date": time_range,
                    "tool": 'mlocarna'
                }).count()

                counts['RNAVIEW'] = db['webservices'].find({
                    "date": time_range,
                    "tool": 'rnaview'
                }).count()

                data.append(counts)
            answer = {
                'header': 'webservices usage',
                'data': data
                }
            self.write_message(answer, binary = False)

    def on_close(self):
        if self in websockets:
            websockets.remove(self)
        for external_tool, ws in external_tools_2_websockets.items():
            if ws == self:
                del external_tools_2_websockets[external_tool]
                break
        print "Client disconnected"

class Application(tornado.web.Application):
    def __init__(self):

        handlers = [
            (r'/', IndexHandler),
            (r'/server', ServerHandler),
            (r'/pyrna', PyrnaDocHandler),
            (r'/rnapedia', RnapediaHandler),
            (r'/account', AccountHandler),
            (r'/register', RegisterUserHandler),
            (r'/login', LoginHandler),
            (r'/logout', LogoutHandler),
            (r'/link', LinkToolHandler),
            (r'/unlink', UnLinkToolHandler),
            (r'/webservices', WebservicesDocHandler),
            (r'/websocket', WebSocketHandler),
            (r'/api/get_key', APIKeyHandler),
            (r'/api/computations/rnafold', RNAfoldHandler),
            (r'/api/computations/rnaplot', RNAplotHandler),
            (r'/api/computations/contrafold', ContrafoldHandler),
            (r'/api/computations/rnaview', RnaviewHandler),
            (r'/api/compute/2d', Compute2dHandler),
            (r'/api/compute/2dplot', Compute2dplotHandler),
            (r'/api/pdb', PDBHandler),
            (r'/api/rna3dhub', RNA3DHubHandler)
        ]

        settings = {
            'template_path': pages_dir,
            'static_path': static_dir,
            'debug': True
        }

        tornado.web.Application.__init__(self, handlers, cookie_secret = str(ObjectId()) , **settings)

if __name__ == '__main__':
    hostname = "127.0.0.1"
    webserver_port = 8080
    mongodb_host = "localhost"
    mongodb_port = 27017

    if "-h" in sys.argv:
        print "Usage: ./server.py [-hostname your_hostname]"
        sys.exit(-1)
    if "-hostname" in sys.argv:
        hostname = sys.argv[sys.argv.index("-hostname")+1]

    with open("%s/../conf/pyrna.json"%os.path.dirname(os.path.realpath(__file__))) as config:
        params = config.read()
        params = json.loads(params)
        webserver_port = params['rest_server']['port']
        mongodb_host = params['mongodb']['host']
        mongodb_port = params['mongodb']['port']
        enabled_algorithms = params['rest_server']['enable-algorithms']

    try :
        mongodb = MongoClient(mongodb_host, mongodb_port)
        logs_db = mongodb['logs']
        webserver_db = mongodb['webserver']
    except Exception, e:
        print '\033[91mI cannot connect any Mongodb instance hosted at %s:%i\033[0m'%(mongodb_host, mongodb_port)
        print 'You can check and change your parameters in the file \033[95mconf/pyrna.json\033[0m'
        sys.exit(-1)

    app = Application()
    server = tornado.httpserver.HTTPServer(app)
    server.listen(webserver_port)
    print "\033[92mYour webserver is now accessible at http://%s:%i/\033[0m"%(hostname, webserver_port)
    print "To modify:"
    print "- the hostname: launch the server like: \033[95m./server.py -hostname your_hostname\033[0m"
    print '- all the other parameters: open and edit the file \033[95mconf/pyrna.conf\033[0m\n'
    main_loop = tornado.ioloop.IOLoop.instance()
    main_loop.start()
