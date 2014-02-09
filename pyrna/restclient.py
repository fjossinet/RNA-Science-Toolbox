#!/usr/bin/env python

#To use this script you will need the following libraries installed:
#- requests (https://github.com/kennethreitz/requests)
#- ujson (http://pypi.python.org/pypi/ujson)
import requests
import ujson
from features import SecondaryStructure, StructuralAlignment
import os
from pandas import DataFrame

class RestClient:

	def __init__(self, baseUrl = "http://arn-ibmc.in2p3.fr/api"):
		self.baseUrl = baseUrl

	def list_databases(self):
		return ujson.loads(requests.get(self.baseUrl+"/db/").text)	

	def get_all_ncRNAs(self, db_name):
		return DataFrame(ujson.loads(requests.get(self.baseUrl+"/db/%s/ncRNAs"%db_name).text))

	def get_all_annotations(self, db_name):
		return DataFrame(ujson.loads(requests.get(self.baseUrl+"/db/%s/annotations"%db_name).text))

	def compute_2d(self, seq = None, name = 'rna', tool = 'rnafold'):
		params = {'tool': tool, 'name': name, 'seq': seq}
		return SecondaryStructure(ujson.loads(requests.get(self.baseUrl+"/compute/2d", params = params).text))

	def compute_2ds(self, data, tool = 'mlocarna'):
		params = {'tool': tool, 'data' : data}
		return StructuralAlignment(ujson.loads(requests.post(self.baseUrl+"/compute/2d", data = params).text))

	def annotate_3d(self, tool = 'rnaview', pdb_id = None, data = None):
		secondary_structures = []
		if pdb_id:
			params = {'tool': tool, 'pdbid': pdb_id}
			for annotation in ujson.loads(requests.get(self.baseUrl+"/compute/2d", params = params).text):
				secondary_structures.append(SecondaryStructure(annotation['2D']))	
		elif data:
			params = {'tool': tool, 'data': data}
			for annotation in ujson.loads(requests.post(self.baseUrl+"/compute/2d", data = params).text):
				secondary_structures.append(SecondaryStructure(annotation['2D']))	
		return secondary_structures	

	def plot_2d(self):
		pass

	def query_ark(self):
		pass

if __name__ == '__main__':
	 
	rest_client = RestClient("http://localhost:8000/api")

	#2D predictions
	print "RNAfold prediction\n------------------------------------------------------------------\n"
	print rest_client.compute_2d("AAGGAACCAGATATAGGACCCACCAGAGATAATTGGG").get_helices()

	print "\n\nContrafold prediction\n------------------------------------------------------------------\n"
	print rest_client.compute_2d("AAGGAACCAGATATAGGACCCACCAGAGATAATTGGG", tool = 'contrafold').get_helices()

	print "\n\nMlocarna alignment\n------------------------------------------------------------------\n"
	alignment = rest_client.compute_2ds(data = open(os.path.dirname(os.path.realpath(__file__))+"/../files/telomerase.fasta", 'r').read(), tool = 'mlocarna')
	print alignment.get_aligned_sequences().values
	print alignment.get_consensus_2d()

	print "\n\n3D annotation of a local PDB file\n------------------------------------------------------------------\n"
	for ss in rest_client.annotate_3d(data = open(os.path.dirname(os.path.realpath(__file__))+"/../files/1MFQ.pdb", 'r').read(), tool = 'rnaview'):
		print ss.get_source()
		print ss.get_helices()
		print ss.get_non_canonical_secondary_interactions()
		print ss.get_tertiary_interactions()
		print ss.get_single_strands()
		print ss.get_junctions()
		print ss.get_modified_ribonucleotides()

	print "\n\n3D annotation of files downloaded from the PDB website\n------------------------------------------------------------------\n"
	for pdb_id in ['1EHZ', '1GID']:
		for ss in rest_client.annotate_3d(tool = 'rnaview', pdb_id = pdb_id):
			print ss.get_source()
			print ss.get_helices()
			print ss.get_non_canonical_secondary_interactions()
			print ss.get_tertiary_interactions()
			print ss.get_single_strands()
			print ss.get_junctions()
			print ss.get_modified_ribonucleotides()