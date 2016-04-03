import os, commands, re, shutil, sys, urllib, subprocess, time, fcntl, urllib, urllib2
from string import maketrans
from pandas import DataFrame
import parsers, utils
from features import RNA, SecondaryStructure, TertiaryStructure
from parsers import base_pairs_to_secondary_structure, parse_bn, to_fasta, to_pdb
from distutils.spawn import find_executable

def get_api_key(rest_server):
    response = urllib.urlopen("http://%s/api/get_key"%rest_server)
    api_key = str(response.read())
    return api_key

class Tool:
    def __init__(self, cache_dir = '/tmp', rest_server = None, api_key = None):
        self.cache_dir = cache_dir
        if not os.path.exists(self.cache_dir):
            os.mkdir(self.cache_dir)
        self.rest_server = rest_server
        self.api_key = api_key

    def find_executable(self, executable):
        if not find_executable(executable):
            raise Exception("%s is not available in your PATH"%executable)

    def submit(self, tool_name, parameters):
        parameters['api_key'] = self.api_key
        parameters = urllib.urlencode(parameters)
        req = urllib2.Request("http://%s/api/computations/%s"%(self.rest_server, tool_name), parameters)
        response = urllib2.urlopen(req)
        output = str(response.read())
        response.close()
        return output

    def download_file(self, uri):
        """
        Download a file from a URI (local file or Web address).

        Parameters:
        -----------
        - uri: an URI. A string starting with 'file://', 'http://' or 'https://'"

        Returns:
        --------
        the local path of the file downloaded.
        """
        if not uri.endswith("/"):
            file_name = self.cache_dir+'/'+uri.split('/')[-1]
        else:
            print "Error: your URI should not end with the character '/'"
            sys.exit(1)
        if uri.startswith("https://") or uri.startswith("http://"):
            urllib.urlretrieve(uri, file_name)
        elif uri.startswith("file://"):
            shutil.copy(uri.split('file://')[1], self.cache_dir)
        else:
            print "Error: your URI must start with 'file://', 'http://' or 'https://'" #if the given path start with /Users/, URI must start with file:///Users/
            sys.exit(1)
        return file_name

class Augustus(Tool):

    """
    Application Controller for Augustus.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("augustus")

    def search(self, molecules, species = "saccharomyces_cerevisiae_S288C"):
        """"
        Parameters:
        ---------
        - molecules: the target molecules that will be used to do the search (as a list of Molecule objects, see pyrna.features)
        - species: the identifier of the species that will be used to do the search (for the list of available species, see the augustus readme file)
        """

        fileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'
        f = open(fileName, 'w')
        f.write(parsers.to_fasta(molecules))
        f.close()

        data = commands.getoutput('cd %s ; augustus --species=%s %s'%(self.cache_dir, species, fileName))

        return data, self.parse_output(data)

    def parse_output(self, output):

        """
        Parameters:
        -----------
        - output: the Augustus output as a String.

        Returns:
        -----------
        A pandas DataFrame describing all the Augustus genes. The index stores genes ids. The column are:
        - sequence_name
        - source
        - genomic_positions
        - score
        - strand
        - frame
        - gene_id
        - protein_sequence
        """

        genes = []
        i = 0
        lines = output.split('\n')
        for line in lines:
            if line.startswith('# start gene'):
                tokens = re.split('\t', lines[i+1])
                gene = {
                    "sequence_name": tokens[0],
                    "source": tokens[1],
                    "genomic_positions": [tokens[3], tokens[4]],
                    "score": tokens[5],
                    "strand": tokens[6],
                    "frame": tokens[7],
                    "gene_id": tokens[8]
                }
                j = 1
                flag = False
                while flag == False:
                    if lines[i+j].startswith('# protein sequence'):
                        protein_sequence = lines[i+j].split('[')[1].split(']')[0]
                        if not ']' in lines[i+j]:
                            k = 1
                            while not ']' in lines[i+j+k]:
                                protein_sequence += (lines[i+j+k].split('# ')[1].split(']')[0])
                                k += 1
                        gene['protein_sequence'] = protein_sequence
                        flag = True
                    j += 1
                genes.append(gene)
            i += 1
        return DataFrame(genes)

class Bcheck(Tool):
    """
    Application Controller for Bcheck.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("Bcheck")

    def check(self, target_molecules):
        """
        Scan the target molecules for RnaseP hits

        Parameters:
        -----------
        - target_molecules: the target molecules (as a list of Molecule objects, see pyrna.features)

        Returns:
        --------
        A pandas DataFrame describing all the hits. The columns are:
        - e_value
        - score
        - target_strand ('+' or '-')
        - target_positions
        - target_name
        - sequence (the hit primary sequence)
        - organism
        """
        fileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'

        with open(fileName, 'w') as f:
            f.write(parsers.to_fasta(target_molecules))

        hits = []
        output = commands.getoutput("Bcheck -f %s"%fileName)

        #a hit is like NW_139495/31959-32288/-1/nuc_fugA, Score = 150.79,  E = 4.491e-49,  GC =  47

        for line in output.split('\n'):
            if re.findall('.+Score = .+', line):
                tokens = re.split('\s+', line)
                sub_tokens = tokens[0].split('/')
                start, end = map(lambda pos: int(pos),sub_tokens[1].split('-'))
                hit = {
                    'e_value': float(tokens[tokens.index('E')+2][:-1]), #[:-1] is to remove the trailing ',' in 'E = 4.491e-49,'
                    'score': float(tokens[tokens.index('Score')+2][:-1]), #[:-1] is to remove the trailing ',' in 'Score = 150.79,'
                    'target_positions': [start, end],
                    'target_strand': '-' if sub_tokens[2] == '-1' else '+'
                }
                for target_molecule in target_molecules:
                    if target_molecule.name == sub_tokens[0]:
                        hit['target_name'] = target_molecule.name
                        hit['organism'] = target_molecule.organism
                        if tokens[-1] == '+':
                            hit['sequence'] = target_molecule.sequence[start-1:end]
                        else:
                            hit['sequence'] = target_molecule.get_complement()[start-1:end][::-1]
                hits.append(hit)

        return DataFrame(hits)

class Blast(Tool):
    def __init__(self, target_molecules, cache_dir = "/tmp", rest_server = None, api_key = None):
        """"
        Parameters:
        ---------
        - target_molecules:  the target molecules that will be used to make the blast database and to do the search (as a list of Molecule objects, see pyrna.features)
        """
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("formatdb")
            self.find_executable("blastall")
        self.target_molecules = target_molecules

    def format_db(self, is_nucleotide=True):
        """
        Dump the target_molecules into a fasta file and format them into a Blast database.

        Parameters:
        ---------
        is_nucleotide (default: True): state if the target molecules are nucleotides or proteins
        """
        with open("%s/input.fasta"%self.cache_dir, 'w+b') as fasta_file:
            fasta_file.write(parsers.to_fasta(self.target_molecules))

        commands.getoutput("cd %s ; formatdb -i %s -p %s -o"%(self.cache_dir, fasta_file.name, "F" if is_nucleotide else "T"))
        self.formatted_db = fasta_file.name

    def parse_output(self, output):
        """
        Parse a blast output

        Parameters:
        ---------
        - output_file: the blast output content as a String

        Returns:
        --------
        A pandas DataFrame describing all the blast hits. The index stores hit ids. The columns are:
        - e_value
        - target_strand ('+' or '-')
        - target_positions
        - target_name
        - query_strand ('+' or '-')
        - query_positions
        - sequence (the hit primary sequence)
        - source
        - organism
        - name
        """
        lines = output.split('\n')
        hits=[]
        i = 0
        query_plus_strand = True
        subject_plus_strand = True
        evalue = None
        source = None
        query_positions = []
        subject_positions = []
        sequence_name = None
        query_name = None
        while i < len(lines):
            line = lines[i].strip()
            #print line
            if line.startswith("BLAST"):
                source = "tool:blast:%s"%line.lower()
            elif line.startswith("Query= "):
                query_name = line.split("Query= ")[1]
            elif line.startswith("Score ="):
                if len(query_positions) and len(subject_positions) : # we have a previous hit to store
                    if not subject_plus_strand:
                        subject_positions = subject_positions[::-1]
                    for m in self.target_molecules:
                        if m.name == sequence_name:
                            if subject_plus_strand:
                                sequence = m.sequence[subject_positions[0][0]-1:subject_positions[-1][1]]
                            else:
                                sequence = m.get_complement()[subject_positions[0][0]-1:subject_positions[-1][1]][::-1]
                            hits.append({
                                "name": query_name,
                                "target_name":sequence_name,
                                "target_positions":[(subject_positions[0][0],subject_positions[-1][1])],
                                "target_strand":'+' if subject_plus_strand else '-',
                                "query_positions":[(query_positions[0][0],query_positions[-1][1])],
                                "query_strand":'+' if query_plus_strand else '-',
                                "e_value":float(evalue),
                                "sequence":sequence,
                                "source": source,
                                "organism": m.organism
                            })
                query_positions = []
                subject_positions = []
                evalue = line.split("Expect =")[1].split(',')[0].strip() #the instruction .split(',') is to be compatible with the blastR output (which is a blastP output)
                if evalue.startswith('e'): #sometimes blast output evalues like 'e-104'
                    evalue = '1'+evalue
                if len(lines[i-3].strip()): #otherwise its a hit for the same sequence
                    sequence_name = lines[i-3][1:].strip()
            elif line.startswith("Strand ="):
                strand_orientations = line.split("Strand = ")[1].strip().split("/")
                query_plus_strand = strand_orientations[0].strip() == "Plus"
                subject_plus_strand = strand_orientations[1].strip() == "Plus"
            elif line.startswith("Query:"):
                tokens = line.split()
                fragment_start = int(tokens[1])
                fragment_end = int(tokens[-1])
                if fragment_start < fragment_end:
                    query_positions.append((fragment_start,fragment_end))
                else:
                    query_positions.append((fragment_end,fragment_start))
            elif line.startswith("Sbjct:"):
                tokens = line.split()
                fragment_start = int(tokens[1])
                fragment_end = int(tokens[-1])
                if fragment_start < fragment_end:
                    subject_positions.append((fragment_start,fragment_end))
                else:
                    subject_positions.append((fragment_end,fragment_start))
            i += 1
        if len(query_positions) and len(subject_positions) : # the last hit to store
            if not subject_plus_strand:
                subject_positions = subject_positions[::-1]
            for m in self.target_molecules:
                if m.name == sequence_name:
                    if subject_plus_strand:
                        sequence = m.sequence[subject_positions[0][0]-1:subject_positions[-1][1]]
                    else:
                        sequence = m.get_complement()[subject_positions[0][0]-1:subject_positions[-1][1]][::-1]
                    hits.append({
                        "name": query_name,
                        "target_name":sequence_name,
                        "target_positions":[(subject_positions[0][0],subject_positions[-1][1])],
                        "target_strand":'+' if subject_plus_strand else '-',
                        "query_positions":[(query_positions[0][0],query_positions[-1][1])],
                        "query_strand":'+' if query_plus_strand else '-',
                        "e_value":float(evalue),
                        "sequence":sequence,
                        "source": source,
                        "organism": m.organism
                    })
                    break
        return DataFrame(hits)

    def blastn(self, query_molecule):
        """
        Blast a query against the formated target molecules

        Parameters:
        -----------
        - query_molecule: a Molecule object (see pyrna.features).

        Returns:
        --------
        A pandas DataFrame describing all the blast hits. The index stores hit ids. The columns are:
        - e_value
        - target_strand ('+' or '-')
        - target_positions
        - target_name
        - query_strand ('+' or '-')
        - query_positions
        - sequence (the hit primary sequence)
        - source
        - organism
        """
        tmp_dir = os.path.dirname(self.formatted_db)
        with open("%s/query.fasta"%tmp_dir, 'w+b') as query_file:
            query_file.write(query_molecule.to_fasta())

        return self.parse_output(commands.getoutput("cd %s ; blastall -p blastn -d %s -i %s"%(self.cache_dir, self.formatted_db, query_file.name)))

    def rpsblast(self):
        pass

class Blastr(Blast):

    """
    Application Controller for blastR.
    """
    def __init__(self, target_molecules, cache_dir = "/tmp", rest_server = None, api_key = None):
        """"
        Parameters:
        ---------
        - target_molecules:  the target molecules that will be used to make the blast database and to do the search (as a list of Molecule objects, see pyrna.features)
        """
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("formatdbR.pl")
            self.find_executable("blastallR.pl")
        self.target_molecules = target_molecules

    def format_db(self):
        """
        Dump the target_molecules into a fasta file and format them into a Blast database.
        """
        with open("%s/input.fasta"%self.cache_dir, 'w+b') as fasta_file:
            fasta_file.write(parsers.to_fasta(self.target_molecules))
        commands.getoutput("cd %s ; formatdbR.pl -i %s"%(self.cache_dir, fasta_file.name))
        self.formatted_db = fasta_file.name

    def blastallr(self, query_molecule):
        """
        Blast a query against the formated target molecules

        Parameters:
        -----------
        - query_molecule: a Molecule object (see pyrna.features).

        Returns:
        --------
        For now, this method returns the raw output of blastr as a String.
        """
        tmp_dir = os.path.dirname(self.formatted_db)
        with open("%s/query.fasta"%tmp_dir, 'w+b') as query_file:
            query_file.write(query_molecule.to_fasta())
        return self.parse_output(commands.getoutput("cd %s ; blastallR.pl -p blastr -i %s -d %s"%(self.cache_dir, query_file.name, self.formatted_db)))

class Bowtie2(Tool):
    """
    Application Controller for Bowtie2.
    """
    def __init__(self, cache_dir = "/tmp", index_path = None, rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("bowtie2-build")
            self.find_executable("bowtie2")
        self.index_path = index_path

    def parse_sam(self, sam_file, target_molecules):
        """
        Parse a SAM file and returns aligned reads.

        Returns:
        --------
        A pandas DataFrame describing the reads. The columns are:
        - genomicStart (an int)
        - genomicEnd (an int)
        - genomicStrand ('+' or '-')
        - genomeName (a String)
        """

        aligned_reads, total_reads, tids  = parsers.parse_sam(sam_file)
        reads = []
        total_aligned_reads = 0
        for aligned_reads_per_genome in aligned_reads:
            total_aligned_reads += len(aligned_reads_per_genome)
            for target_molecule in target_molecules:
                if target_molecule.name == tids[aligned_reads_per_genome[0]['tid']]:
                    current_molecule = target_molecule
            genomeName = tids[aligned_reads_per_genome[0]['tid']]

            for aligned_read in aligned_reads_per_genome :
                reads.append({
                    'genomicStart': aligned_read['genomicStart'],
                    'genomicEnd': aligned_read['genomicEnd'],
                    'genomeName': genomeName,
                    'genomicStrand': aligned_read['genomicStrand']
                })
        print "%i reads found, %i reads aligned..."%(total_reads, total_aligned_reads)
        return DataFrame(reads)

    def align(self, target_molecules, fastq_file, no_parsing = False, user_defined_options=[]):
        """
        Align reads against target molecules.

        Parameters:
        -----------
        - target_molecules: the genomic sequences to be used for the alignment (an array of Molecule objects, see pyrna.features)
        - fastq_file: the full path for the fastq file containing the reads (as a String)
        - no_parsing (default: False): if True, the function returns the full path of the SAM file without parsing it

        Returns:
        --------
        The full path of the SAM file or a pandas DataFrame describing the reads. The columns are:
        - genomicStart (an int)
        - genomicEnd (an int)
        - genomicStrand ('+' or '-')
        - genomeName (a String)

        Once done, the full path to the directory containing the index plus the prefix of the index files is recorded in an attribute named index_path:
        bowtie2 = Bowtie2()
        bowtie2.align(....)
        print bowtie2.index_path
        """

        result_file = self.cache_dir+'/'+os.path.basename(fastq_file)+'.sam'

        if not self.index_path:
            self.index_path = self.build_index(target_molecules)

        elif not os.path.exists(self.index_path+".1.bt2"):
            self.build_index(target_molecules)

        print "bowtie2 %s -x %s -q \"%s\" -S %s"%(' '.join(user_defined_options), self.index_path, fastq_file, result_file)
        commands.getoutput("bowtie2 %s -x %s -q \"%s\" -S %s"%(' '.join(user_defined_options), self.index_path, fastq_file, result_file))
        print "SAM file %s produced successfully!!"%result_file

        if no_parsing:
            return result_file

        return self.parse_sam(result_file, target_molecules)

    def build_index(self, target_molecules):
        """
        Build a new index for the target molecules

        Parameters:
        -----------
        - target_molecules: the genomic sequences to be used for the index (an array of Molecule objects, see pyrna.features)

        Returns:
        --------
        The full path to the directory containing the index plus the prefix of the index files
        """
        if not self.index_path:
            random_name = utils.generate_random_name(7)
            self.index_path = self.cache_dir+"/"+random_name
            fasta_file_name = self.cache_dir+'/'+random_name+'.fa'
            fasta_file = open(fasta_file_name, 'w')
            fasta_file.write(to_fasta(target_molecules))
            fasta_file.close()
            print "bowtie2-build %s %s"%(fasta_file_name, self.index_path)
            commands.getoutput("bowtie2-build %s %s"%(fasta_file_name, self.index_path))
            print "Index files produced successfully!!"
        else:
            random_name = self.index_path.split('/')[-1]
            fasta_file_name = self.index_path+".fa"
            fasta_file = open(fasta_file_name, 'w')
            fasta_file.write(to_fasta(target_molecules))
            fasta_file.close()
            print "bowtie2-build %s %s"%(fasta_file_name, self.index_path)
            commands.getoutput("bowtie2-build %s %s"%(fasta_file_name, self.index_path))
            print "Index files produced successfully!!"

        return self.index_path


class Clustalw(Tool):

    """
    Application Controller for clustalw.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("clustalw2")

    def align(self, molecules):
        """
        Parameters:
        ---------
        - molecules: the molecules to align (as a list of Molecule objects, see pyrna.features)

        Returns:
        --------
        the clustalw output as a String and an array of aligned molecules
        """
        fileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'

        with open(fileName, 'w') as f:
            f.write(parsers.to_fasta(molecules))

        commands.getoutput("clustalw2 -infile=%s"%fileName)

        with open(fileName.split('.fasta')[0]+".aln") as output_file:
            data = output_file.read()

        return data, self.parse_output(data, molecules)

    def parse_output(self, output, molecules):
        """
        Return a list of aligned molecules

        Parameters:
        -----------
        - output: the clustalw raw output as a String
        - molecules: the list of molecules that have been used to produce the output (as a list of Molecule objects, see pyrna.features)
        """
        aligned_molecules = {}

        for molecule in molecules:
            aligned_molecules[molecule.name] = ""

        for line in output.split('\n'):
            tokens = re.split('\s+',line)
            if len(tokens) == 2 and aligned_molecules.has_key(tokens[0]):
                aligned_molecules[tokens[0]] += tokens[1]

        rnas = []
        for k,v in aligned_molecules.iteritems():
            rnas.append(RNA(name=k, sequence=v))

        return rnas

class Cmalign(Tool):
    """
    Application Controller for cmalign.
    """
    def __init__(self, cache_dir = "/tmp", local_mode = True, rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("cmalign")
        self.local_mode = local_mode

    def align(self, molecules, rfam_id = None, rfam = None, stockholm_content = None, cm_content = None):
        """
        Align new ncRNA candidates (the molecules object argument) to the RFAM family (defined by the rfam_id argument) using the CM model stored in cm_file.

        Parameters:
        -----------
        - molecules: a list of molecules to align
        - rfam_id: the id of the RFAM family to use for the alignment
        - rfam: an Rfam object

        Returns:
        --------
        a tuple like: (list of all the aligned molecules, dict of organism names (keys) and accession numbers/start-end (values), Dataframe of the consensus 2D)
        """
        with open("%s/input.fasta"%self.cache_dir, 'w') as fasta_file:
            fasta_file.write(parsers.to_fasta(molecules))
        if not stockholm_content:
            try:
                stockholm_content = rfam.get_entry(rfam_id, format='stockholm')
            except Exception, e:
                raise e
        if rfam_id:
            stockholm_file = open("%s/%s.stk"%(self.cache_dir,rfam_id), 'w')
        else:
            stockholm_file = open("%s/%s.stk"%(self.cache_dir, utils.generate_random_name(7)), 'w')

        stockholm_file.write(stockholm_content)
        stockholm_file.close()

        output = None

        cm_file = None

        if cm_content:
            cm_file = self.cache_dir+'/'+utils.generate_random_name(7)+'.cm'

            with open(cm_file, 'w') as f:
                f.write(cm_content)

        if self.local_mode:
            output = commands.getoutput("cmalign -l --withali %s %s %s"%(stockholm_file.name, cm_file if cm_file else rfam.cache_dir+'/CMs/'+rfam_id+".cm", fasta_file.name))
        else:
            output = commands.getoutput("cmalign --withali %s %s %s"%(stockholm_file.name, cm_file if cm_file else rfam.cache_dir+'/CMs/'+rfam_id+".cm", fasta_file.name))
        #shutil.rmtree(self.cache_dir)
        if rfam_id:
            output = "#=GF AC "+rfam_id+"\n"+output
        return parsers.parse_stockholm(output)

class Cmbuild(Tool):
    """
    Application Controller for Cmbuild.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("cmbuild")

    def build(self, stockholm_content):
        """
        Parameters:
        ---------
        -stockholm_content:

        Returns:
        --------
        the content of a covariance model (as a String)
        """
        name = utils.generate_random_name(7)
        stockholm_file = self.cache_dir+'/'+name+'.sto'
        cm_file = self.cache_dir+'/'+name+'.cm'

        with open(stockholm_file, 'w') as f:
            f.write(stockholm_content)

        commands.getoutput("cmbuild %s %s "%(cm_file, stockholm_file))

        cm_content = None
        with open(cm_file) as f:
            cm_content = f.read()

        return cm_content

class Cmcalibrate(Tool):
    """
    Application Controller for Cmcalibrate.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("cmcalibrate")

    def calibrate(self, cm_content):
        """
        Parameters:
        ---------
        -cm_content:

        Returns:
        --------
        the calibrated content of a covariance model (as a String)
        """
        name = utils.generate_random_name(7)
        cm_file = self.cache_dir+'/'+name+'.cm'

        with open(cm_file, 'w') as f:
            f.write(cm_content)

        commands.getoutput("cmcalibrate %s"%cm_file)

        with open(cm_file) as f:
            cm_content = f.read()

        return cm_content

class Cmsearch(Tool):
    """
    Application Controller for Cmsearch.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("cmsearch")

    def search(self, molecules, rfam_id = None, rfam = None, cm_content = None,  gathering_threshold = True):
        """
        Launch a search with cmsearch

        Parameters:
        -----------
        - molecules: the molecules used to do the search
        - rfam_id : to id of the RFAM family
        - rfam: an Rfam object (see pyrna.db)
        - cm_content: the content of a CM file as a String (default is None).

        Returns:
        --------
        A pandas DataFrame describing all the cmsearch hits. The index stores hit ids. The column are:
        - e_value
        - p_value
        - score
        - target_strand ('+' or '-')
        - sequence (the hit primary sequence)
        - target_positions
        - target_name
        - query_positions
        - cm_file (corresponds to the query sequence, but its a covariance models for cmsearch)
        - RFAM_family
        - source
        - organism
        """
        if cm_content:
            fileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'

            with open(fileName, 'w') as f:
                f.write(parsers.to_fasta(molecules))

            cmFile = self.cache_dir+'/'+utils.generate_random_name(7)+'.cm'

            with open(cmFile, 'w') as f:
                f.write(cm_content)

            if not gathering_threshold:
                return self.parse_output(commands.getoutput("cmsearch "+cmFile+" "+fileName), molecules, False)
            else:
                return self.parse_output(commands.getoutput("cmsearch --ga "+cmFile+" "+fileName), molecules, True)
        else:
            #write molecules as FASTA
            fileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'

            with open(fileName, 'w') as f:
                f.write(parsers.to_fasta(molecules))

            if not gathering_threshold:
                return self.parse_output(commands.getoutput("cmsearch "+rfam.cache_dir+"/CMs/"+rfam_id+".cm "+fileName), molecules, False)
            else:
                return self.parse_output(commands.getoutput("cmsearch --ga "+rfam.cache_dir+"/CMs/"+rfam_id+".cm "+fileName), molecules, True)

    def parse_output(self, output, molecules, gathering_threshold = True):
        """
        Parse the output of the cmsearch tool.

        Parameters:
        -----------
        - output: the cmsearch output as a String
        - molecules: the molecules used to do the search

        Returns:
        --------
        A pandas DataFrame describing the cmsearch hits. The index stores hit ids. The column are:
        - e_value
        - p_value
        - score
        - target_strand ('+' or '-')
        - sequence (the hit primary sequence)
        - target_positions
        - target_name
        - query_positions
        - cm_file (corresponds to the query sequence, but its a covariance models for cmsearch)
        - RFAM_family
        - source
        - organism
        """
        flag = False
        hits=[]
        i = 0
        plus_strand = True
        lines = output.split('\n');
        source = None
        while i < len(lines):
            if lines[i].startswith("# INFERNAL"):
                source = lines[i][1:].strip()
                flag = True
            elif lines[i].startswith("# command:"):
                cm_file = lines[i].split("cmsearch --ga")[1].split()[0] if gathering_threshold else lines[i].split("cmsearch")[1].split()[0]
            elif lines[i].strip().startswith("CM:"):
                rfam_family_id = lines[i].split("CM:")[1].strip()
            elif lines[i].strip().startswith("Plus strand results:"):
                plus_strand = True
                if len(lines[i-2].strip()) != 0: #otherwise its Plus Strand results after Minus Strand results for the same sequence
                    sequence_name = lines[i-2][1:].strip()
            elif lines[i].strip().startswith("Minus strand results:"):
                plus_strand = False
                if len(lines[i-2].strip()) != 0: #otherwise its Minus Strand results after Plus Strand results for the same sequence
                    sequence_name = lines[i-2][1:].strip()
            elif lines[i].strip().startswith("Query ="):
                tokens = lines[i].split(", Target = ")[0].strip().split("Query = ")[1].split(" - ")
                query_start = int(tokens[0])
                query_end = int(tokens[1])
                tokens = lines[i].split(", Target = ")[1].strip().split(" - ")
                target_start = int(tokens[0])
                target_end = int(tokens[1])
                i+=1
                scores = lines[i].split(", ")
                score = None
                p_value = None
                e_value = None
                for s in scores:
                    if s.startswith("Score"):
                        score = float(s.strip().split(" = ")[1].strip())
                    elif s.startswith("E"):
                        e_value = float(s.strip().split(" = ")[1].strip())
                    elif s.startswith("P"):
                        p_value = float(s.strip().split(" = ")[1].strip())
                i+=3
                query_positions = []
                target_positions = []
                current_query_start = query_start
                current_target_start = target_start
                query_sequence = ""
                target_sequence = ""
                target_turn = False #we can have situation where the target and query lines have the same start position
                while True:
                    if not target_turn:
                        tokens = lines[i].strip().split()
                        seq = ''.join(tokens[1:-1])
                        query_sequence += seq
                        target_turn = True
                        i+=2
                    elif target_turn:
                        tokens = lines[i].strip().split()
                        seq = ''.join(tokens[1:-1])
                        target_sequence += seq
                        target_turn = False
                        if lines[i].strip().endswith("%i"%target_end) and lines[i-2].strip().endswith("%i"%query_end)  or len(lines[i+2].strip()) == 0:
                            break
                        i += 3
                current_query_start = query_start
                current_target_start = target_start
                subsequences = query_sequence.split("*")
                if len(subsequences) == 1:
                    residues = len(subsequences[0])-subsequences[0].count('.')-subsequences[0].count('-')
                    query_positions.append((current_query_start,current_query_start + residues-1))
                else:
                    for subsequence in subsequences:
                        if subsequence.startswith("["):
                            current_query_start += int(subsequence[1:-1].strip())
                        else:
                            residues = len(subsequence)-subsequence.count('.')-subsequence.count('-')
                            query_positions.append((current_query_start,current_query_start+residues-1))
                            current_query_start += residues
                subsequences = target_sequence.split("*")
                if len(subsequences) == 1:
                    residues = len(subsequences[0])-subsequences[0].count('.')-subsequences[0].count('-')
                    if plus_strand:
                        target_positions.append((current_target_start,current_target_start + residues-1))
                    else:
                        target_positions.append((current_target_start - residues+1, current_target_start))
                else:
                    for subsequence in subsequences:
                        if subsequence.startswith("["):
                            if plus_strand:
                                current_target_start += int(subsequence[1:-1].strip())
                            else:
                                current_target_start -= int(subsequence[1:-1].strip())
                        else:
                            residues = len(subsequence)-subsequence.count('.')-subsequence.count('-')
                            if plus_strand:
                                target_positions.append((current_target_start,current_target_start+residues-1))
                                current_target_start += residues
                            else:
                                target_positions.append((current_target_start-residues+1,current_target_start))
                                current_target_start -= residues
                if e_value != None:
                    for m in molecules:
                        if m.name == sequence_name:
                            if plus_strand:
                                sequence = m.sequence[target_positions[0][0]-1:target_positions[-1][1]]
                            else:
                                target_positions = target_positions[::-1]
                                sequence = m.get_complement()[target_positions[0][0]-1:target_positions[-1][1]][::-1]
                            hit = {
                                "cm_file": cm_file,
                                "RFAM_family": rfam_family_id,
                                "target_name": sequence_name,
                                "target_positions": target_positions,
                                "sequence": sequence,
                                "query_positions": query_positions,
                                "target_strand":'+' if plus_strand else '-',
                                "source": source,
                                "organism": m.organism
                            }

                            if score != None:
                                hit["score"] = score
                            if e_value != None:
                                hit["e_value"] = e_value
                            if p_value != None:
                                hit["p_value"] = p_value
                            hits.append(hit)
                            break
                else:
                    raise Exception("No e-value found. The covariance model %s was probably not calibrated"%cm_file)
            i+=1
        index= []
        for i in range(1,len(hits)+1):
            index.append("%s_candidate_%i"%(rfam_family_id,i))
        if not flag:
            raise Exception("No Cmsearch output")
        return DataFrame(hits, index)

class Contrafold(Tool):
    """
    Application Controller for Contrafold.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("contrafold")

    def fold(self, molecule, raw_output = False):
        """
        Parameters:
        ---------
        - molecule: a pyrna.features.Molecule object (RNA or DNA)
        - raw_output (default: False): if True, the function returns the raw output instead of a pandas DataFrame

        Returns:
        -------
        the secondary structure as a list of base-pairs in a pandas DataFrame
        """
        output = None
        if self.rest_server:
            values = {
                'name' : molecule.name,
                'sequence': molecule.sequence,
                'api_key': self.api_key
            }
            data = urllib.urlencode(values)
            req = urllib2.Request("http://%s/api/computations/contrafold"%self.rest_server, data)
            response = urllib2.urlopen(req)
            output = str(response.read())
            response.close()
        else:
            fileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'
            with open(fileName, 'w') as fasta_file:
                fasta_file.write(parsers.to_fasta([molecule]))
            output = commands.getoutput("cd %s ; contrafold predict %s"%(self.cache_dir, fileName)).strip()
        filtered_lines = []
        for line in output.split('\n'):
            if not line.startswith('>structure'):
                filtered_lines.append(line)
        output = '\n'.join(filtered_lines)
        if raw_output:
            return output
        else:
            rnas, base_pairs = parsers.parse_vienna(output)
            return base_pairs[0]

class Cufflinks(Tool):
    """
    Application Controller for Cufflinks.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("cufflinks")

    def __cmp(self, dict1, dict2):
        """
        Comparison function for the sort() method of a list.

        first sorting: by growing genome name
        second sorting: by growing start
        third sorting: by decreasing end
        fourth and fifth sorting: by class: 1-gene, 2-mRNA or other, 3-CDS
        """
        if dict1['genomeName'] > dict2['genomeName']:
            return 1
        elif dict1['genomeName'] < dict2['genomeName']:
            return -1
        else:
            if dict1['start'] > dict2['start']:
                return 1
            elif dict1['start'] < dict2['start']:
                return -1
            else:
                if dict1['end'] > dict2['end']:
                    return -1
                elif dict1['end'] < dict2['end']:
                    return 1
                else:
                    if dict1['class'] == 'gene' and dict2['class'] != 'gene':
                        return -1
                    elif dict1['class'] != 'gene' and dict2['class'] == 'gene':
                        return 1
                    else:
                        if dict1['class'] != 'CDS' and dict2['class'] == 'CDS':
                            return -1
                        elif dict1['class'] == 'CDS' and dict2['class'] != 'CDS':
                            return 1
                        else:
                            return 0

    def predict_genes(self, db_name, bam_file, db_host = "localhost", db_port = 27017, annotations = None):
        """
        Method that runs the Cufflinks program with RABT assembly option

        Parameters:
        ----------
        - db_name: the name of the database containing the gene annotations to be used as a reference (needed to characterize the new genes for example)
        - bam_file: the full path of a BAM sorted file
        - db_host : the host for the MongoDB containing the gene annotations
        - db_port : the port for the MongoDB containing the gene annotations
        - annotations : an array of annotations used to be dumped as a GFF file
        """

        if not annotations:
            data_str = self.__mongo_to_gff3(db_name, db_host, db_port)
        else:
            data_str = ""
            for annotation in annotations:
                data_str += "%s\t.\t%s\t%i\t%i\t0.0\t%s\t.\tid %s\n"%(annotation['genomeName'], annotation['class'], annotation['genomicPositions'][0], annotation['genomicPositions'][1], annotation['genomoicStrand'], annotation['_id'])

        gff_file = utils.generate_random_name(7)+'.gff'
        fh = open('%s/%s'%(self.cache_dir, gff_file), 'w')
        fh.write(data_str)
        fh.close()

        commands.getoutput("cd %s ; cufflinks -g %s %s"%(self.cache_dir, gff_file, bam_file))
        df = self.__parse_output()
        return df

    def __mongo_to_gff3(self, db_name, db_host = "localhost", db_port = 27017):
        """
        Converts annotations from a Mongo database into gff3 data as Cufflinks like to see them. More precisely, this means that Cufflinks needs:
        - parent-child relation gene, mRNA, CDS, exons, introns
        - if introns are in phase 0, 1 or 2

        Parameters:
        ----------
        - db_name: database name
        - db_host: database host
        - db_port: database port

        Returns:
        -------
        The gff3 data as a String.
        """
        from pyrna.db import charnDB
        charnDB = charnDB(db_host, db_port)
        db = charnDB.get_database(db_name)

        ##RETRIEVING annotations from the database
        annotation_list = []
        intron_list = []
        for annotation in db['annotations'].find({'source': {'$regex': '^db:ncbi:.+$'}}): #on the server 'charn'
            annotation_dict = {'genomeName': annotation['genomeName'],
                       'source': annotation['source'],
                       'class': annotation['class'],
                       'start': annotation['genomicPositions'][0],#Integer
                       'end': annotation['genomicPositions'][-1],#Integer
                       'score': str(annotation.get('score', '.')),#String
                       'genomicStrand': annotation['genomicStrand'],
                       'name': annotation.get('product', annotation.get('locus_tag', 'unknown'))
                       }
            if annotation['class'] != 'intron':
                annotation_list.append(annotation_dict)
            else:
                intron_list.append(annotation_dict)
        charnDB.disconnect()

        ##SORTING annotations
        annotation_list.sort(self.__cmp)

        if intron_list:
            print "Total number of introns in the database %s: %i"%(db_name, len(intron_list))
            intron_list.sort(self.__cmp)
            intron_counter = 0
        else:
            print "No intron in the database %s"%db_name

        ##GFF3 formatting
        gene = -1
        rna = -1
        exon = -1
        cds = -1
        data = "##gff-version 3\n"
        for annotation in annotation_list:
            chain = "%s\t%s\t%s\t%i\t%i\t%s\t%s\t.\t"%(annotation['genomeName'], annotation['source'], annotation['class'], annotation['start'], annotation['end'], annotation['score'], annotation['genomicStrand'])
            if annotation['class'] == 'gene':
                gene += 1
                gene_parent = "gene%s"%gene
                data += chain+"ID=%s;Name=%s\n"%(gene_parent, annotation['name'])
            elif annotation['class'] in ['mRNA','tRNA','rRNA','ncRNA']:
                rna += 1
                rna_parent = "rna%s"%rna
                data += chain+"ID=%s;Name=%s;Parent=%s\n"%(rna_parent, annotation['name'], gene_parent)
                retained_introns = []
                if intron_list:#list that discharges progressively
                    introns_copy = intron_list[:]
                    j = 0
                    for i in range (0, len(introns_copy)):
                        intron = introns_copy[i]
                        if intron['genomeName'] == annotation['genomeName'] and intron['genomicStrand'] == annotation['genomicStrand']:
                            if intron['start'] > annotation['start'] and intron['end'] < annotation['end']:
                                intron_counter += 1
                                retained_introns.append(intron)
                                intron_list.pop(i-j)
                                j += 1
                                if not intron_list:
                                    print "Processing of all intron data (%i introns): done"%intron_counter
                if retained_introns:
                    exon_coord = []
                    for i in range(0, len(retained_introns)):
                        if i == 0:
                            exon_coord.append((annotation['start'], retained_introns[i]['start'] -1))
                        if i > 0:
                            exon_coord.append((retained_introns[i-1]['end'] +1, retained_introns[i]['start'] -1))
                        if i == len(retained_introns)-1:
                            exon_coord.append((retained_introns[i]['end'] +1, annotation['end']))
                    if annotation['genomicStrand'] == '-':
                        exon_coord.reverse()
                    for coordinates in exon_coord:
                        exon += 1
                        data += "%s\t%s\texon\t%i\t%i\t%s\t%s\t.\tID=id%s;Parent=%s\n"%(annotation['genomeName'], annotation['source'], coordinates[0], coordinates[1], annotation['score'], annotation['genomicStrand'], exon, rna_parent)
                    cds += 1
                    for i in range(0, len(exon_coord)):
                        coordinates = exon_coord[i]
                        exon_len = (coordinates[1] - coordinates[0])+1
                        if i == 0:
                            exon1_len = exon_len
                            phase = 0
                        elif 0 < i <= len(exon_coord)-1:
                            exon2_len = exon_len
                            if exon1_len%3 == 0:
                                phase = 0
                            elif exon1_len%3 == 1:
                                phase = 2
                            elif exon1_len%3 == 2:
                                phase = 1
                            exon1_len = exon1_len + exon2_len
                        if i == len(exon_coord)-1:
                            if exon1_len%3 != 0:
                                print "Warning: the RNA on genome %s at coordinates [%i,%i] has a length (%i) that is not a multiple of 3"%(annotation['genomeName'],annotation['start'],annotation['end'],exon1_len)
                                sys.exit()
                        data += "%s\t%s\tCDS\t%i\t%i\t%s\t%s\t%i\tID=cds%s;Parent=%s\n"%(annotation['genomeName'], annotation['source'], coordinates[0], coordinates[1], annotation['score'], annotation['genomicStrand'], phase, cds, rna_parent)
                else:
                    exon += 1
                    data += "%s\t%s\texon\t%i\t%i\t%s\t%s\t.\tID=id%s;Parent=%s\n"%(annotation['genomeName'], annotation['source'], annotation['start'], annotation['end'], annotation['score'], annotation['genomicStrand'], exon, rna_parent)
                    cds += 1
                    data += "%s\t%s\tCDS\t%i\t%i\t%s\t%s\t0\tID=cds%s;Parent=%s\n"%(annotation['genomeName'], annotation['source'], annotation['start'], annotation['end'], annotation['score'], annotation['genomicStrand'], cds, rna_parent)

        return data

    def __parse_output(self):
        """
        Method that parses the Cufflinks output file 'isoforms.fpkm_tracking'

        Returns:
        --------
        A pandas DataFrame describing all the Cufflinks genes. The columns are:
        - gene_names: list of gene names or '-'
        - genome
        - genomicPositions (Remark: coordinate start in the isoforms.fpkm_tracking file is incorrect (-1))
        - coverage
        - score: FPKM (Fragments Per Kilobase of exon per Million fragments mapped)
        - status: four different gene types are possible: 'expressed' or 'unexpressed' or 'isoform' or 'novel'
        """
        fh = open("%s/isoforms.fpkm_tracking"%self.cache_dir, 'r')
        lines = fh.readlines()
        transcripts = []
        print "Total number of genes in the Cufflinks output (file isoforms.fpkm_tracking): %i"%len(lines[1:])
        for line in lines[1:]:#first line contains the description of columns
            tokens = line[:-1].split('\t')
            locus = tokens[6].split(':')
            coordinates = locus[1].split('-')
            transcript = {
                    'tracking_id': tokens[0],#provisional key to be deleted
                    'gene_id': tokens[3],#provisional key to be deleted
                    'gene_names': tokens[4],
                    'genome': locus[0],
                    'genomicPositions': [int(coordinates[0])+1,int(coordinates[1])],
                    'coverage': float(tokens[8]),
                    'score': float(tokens[9])
                    }
            if not transcripts:
                transcripts.append([transcript])
            else:
                if transcripts[-1][0]['gene_id'] == transcript['gene_id']:
                    transcripts[-1].append(transcript)
                else:
                    transcripts.append([transcript])

        fh.close()
        total_genes = []#list that will contain total number of annotated genes (expressed or not), isoforms and novel genes
        for liste in transcripts:
            if len(liste) > 1:#isoform and its corresponding gene(s) that are referenced with the same CUFF gene_id (transcript which can overlap one or several genes)
                names = []
                iso = []
                for transcript in liste:
                    if transcript['tracking_id'].startswith('CUFF.'):
                        iso.append(transcript)
                    else:
                        names.append(transcript['gene_names'])
                        if transcript['coverage'] != 0:
                            transcript['class'] = 'expressed'
                        else:
                            transcript['class'] = 'unexpressed'
                        del transcript['gene_id']
                        del transcript['tracking_id']
                        total_genes.append(transcript)
                if iso:
                    for transcript in iso:
                        transcript['gene_names'] = ','.join(names)
                        transcript['class'] = 'isoform'
                        del transcript['gene_id']
                        del transcript['tracking_id']
                        total_genes.append(transcript)
            else:
                transcript = liste[0]
                if transcript['tracking_id'].startswith('CUFF.'):#novel gene
                    transcript['class'] = 'novel'
                else:#other genes expressed or not (neither isoform nor novel gene)
                    if transcript['coverage'] != 0:
                        transcript['class'] = 'expressed'
                    else:
                        transcript['class'] = 'unexpressed'
                del transcript['gene_id']
                del transcript['tracking_id']
                total_genes.append(transcript)

        return DataFrame(total_genes)

class Gmorse(Tool):
    """
    Application Controller for Gmorse.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("extractNonMappedReads")
            self.find_executable("coverage")
            self.find_executable("build_covtigs")
            self.find_executable("gmorse")

    def extract_unmapped_reads(self, aligned_reads_file, fastq_file):
        print "extractNonMappedReads %s %s"%(aligned_reads_file, fastq_file)
        output = commands.getoutput("extractNonMappedReads %s %s"%(aligned_reads_file, fastq_file))
        unmapped_reads_file = self.cache_dir+'/unmapped_reads.fa'
        with open(unmapped_reads_file, 'w') as f:
            f.write (output)
        return unmapped_reads_file

    def calculate_depth_coverage(self, aligned_reads_file):
        print "coverage %s"%(aligned_reads_file)
        output = commands.getoutput("coverage %s"%(aligned_reads_file))
        depth_coverage_file = self.cache_dir+'/depth_coverage'
        with open(depth_coverage_file, 'w') as f:
            f.write (output)
        return depth_coverage_file

    def build_covtigs(self, depth_coverage_file, depth_treshold):
        print "build_covtigs %s %i"%(depth_coverage_file, depth_treshold)
        output = commands.getoutput("build_covtigs %s %i"%(depth_coverage_file, depth_treshold))
        covtigs_file = self.cache_dir+'/covtigs'
        with open(covtigs_file, 'w') as f:
            f.write (output)
        return covtigs_file

    def make_model(self, unmapped_reads_file, covtigs_file, target_molecules):
        scaffolds_file = self.cache_dir+'/scaffolds.fa'
        output_gene_models = self.cache_dir+'/gene_model'
        extended_covtigs = self.cache_dir+'/extended_covtigs'
        validated_junctions = self.cache_dir+'/validated_junctions'
        with open(scaffolds_file, 'w') as f:
            f.write (to_fasta(target_molecules))
        print "gmorse -r %s -c %s -f %s -G %s -C %s -J %s"%(unmapped_reads_file, covtigs_file, scaffolds_file, output_gene_models, extended_covtigs, validated_junctions)
        commands.getoutput("gmorse -r %s -c %s -f %s -G %s -C %s -J %s"%(unmapped_reads_file, covtigs_file, scaffolds_file, output_gene_models, extended_covtigs, validated_junctions))
        o = open(output_gene_models)
        model = o.read()
        o.close()
        return self.parse_model(model)

    def parse_model(self, model):
        """
        Parameters:
        -----------
        - output: the Gmorse output as a String.

        Returns:
        -----------
        A pandas DataFrame describing the Gmorse gene models. The columns are:
        - the genome name
        - the class of the annotaton (UTR, exon, CDS, mRNA)
        - the genomic strand ('+' or '-')
        - the genomic positions
        """
        hits = []
        lines = model.split('\n')
        lines = lines[1:-1]
        for line in lines:
            tokens = re.split('\t', line)
            hit = {
                "genome": tokens[0],
                "class": tokens[2],
                "genomicPositions": [int(tokens[3]), int(tokens[4])],
                "genomicStrand": tokens[6]
            }
            hits.append(hit)
        return DataFrame(hits)

class Gotohscan(Tool):
    """
    Application Controller for GotohScan.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("GotohScan2a")

    def scan(self, query_molecule, target_molecules, evalue = 1e-3):
        """
        Scan the target molecules for hits of query molecule

        Parameters:
        -----------
        - query_molecule: the query molecule (as a Molecule object, see pyrna.features)
        - target_molecules: the target molecules (as a list of Molecule objects, see pyrna.features)
        - evalue (default: 1e-3): only hits with a evalue better or equal to this value will be returned

        Returns:
        --------
        A pandas DataFrame describing all the hits. The columns are:
        - e_value
        - score
        - target_strand ('+' or '-')
        - target_positions
        - target_name
        - sequence (the hit primary sequence)
        - organism
        """

        hits = []
        queryFileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'
        with open(queryFileName, 'w') as query_file:
            query_file.write(parsers.to_fasta([query_molecule]))

        targetFileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'
        with open(targetFileName, 'w') as target_file:
            target_file.write(parsers.to_fasta(target_molecules))

        output = commands.getoutput("GotohScan2a -e %e -d %s -q %s"%(evalue, targetFileName, queryFileName))

        for line in output.split('\n'):
            if line.startswith(query_molecule.name):
                tokens = re.split('\s+',line)
                start = int(tokens[-5])
                end = int(tokens[-4])
                hit = {
                    'e_value': float(tokens[-3]),
                    'score' : float(tokens[-2]),
                    'target_positions': [start, end] if tokens[-1] == '+' else [end, start],
                    'target_strand': tokens[-1]
                }
                for target_molecule in target_molecules:
                    if target_molecule.name == tokens[1]:
                        hit['target_name'] = target_molecule.name
                        hit['organism'] = target_molecule.organism
                        if tokens[-1] == '+':
                            hit['sequence'] = target_molecule.sequence[start-1:end]
                        else:
                            hit['sequence'] = target_molecule.get_complement()[end-1:start][::-1]
                hits.append(hit)

        return DataFrame(hits)

class Mlocarna(Tool):

    """
    Application Controller for mlocarna.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("mlocarna")

    def align(self, molecules):
        """
        Returns:
        --------
        a tuple like (list of aligned molecules, secondary structure computed as a list of base-pairs in a pandas DataFrame)
        """
        fileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'
        with open(fileName, 'w') as fasta_file:
            fasta_file.write(parsers.to_fasta(molecules))
        output = commands.getoutput("mlocarna %s"%fileName)

        aligned_molecules = {}
        consensus2D = None

        for molecule in molecules:
            aligned_molecules[molecule.name] = ""

        for line in output.split('\n'):
            tokens = re.split('\s+', line)
            if len(tokens) == 2 and aligned_molecules.has_key(tokens[0]):
                aligned_molecules[tokens[0]] += tokens[1]
            elif tokens[0] == 'alifold':
                consensus2D = parsers.parse_bn(tokens[1])

        rnas = []
        for k,v in aligned_molecules.iteritems():
            rnas.append(RNA(name=k, sequence=v))

        return (rnas, consensus2D)

class RnaAlifold(Tool):

    """
    Application Controller for RNAalifold.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("RNAalifold")

    def align(self, alignment):
        """
        Returns:
        --------
        the bracket notation of the MFE structure as a String
        """
        fileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.aln'
        with open(fileName, 'w') as aln_file:
            aln_file.write(alignment)
        return commands.getoutput("cd %s ; RNAalifold < %s"%(self.cache_dir, fileName)).strip().split('\n')[-1].split(' ')[0]

class Rnafold(Tool):

    """
    Application Controller for RNAfold.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("RNAfold")

    def fold(self, molecule, constraints = None, bp_probabilities = False, raw_output = False):
        """
        Parameters:
        -----------
        - molecule: a pyrna.features.Molecule object (RNA or DNA)
        - constraints: a string defining the constraints. See the RNAfold documentation for the notation to be used.
        - bp_probabilities (default: False): if True, the method returns the base-pair probabilities
        - raw_output (default: False): if True, the method returns the raw output instead of the pandas Dataframe.

        Returns:
        --------
        - base-pair probabilities in a pandas DataFrame (if parameter bp_probabilities is True)
        - a secondary structure as a list of base-pairs in a pandas DataFrame (if parameter bp_probabilities is False)
        """
        if self.rest_server:
            parameters = {
                'name' : molecule.name,
                'sequence': molecule.sequence,
                'function': 'fold',
                'bp_probabilities': bp_probabilities
            }
            if constraints:
                parameters['constraints'] = constraints
            output = self.submit('rnafold', parameters)
        else:
            fileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'
            with open(fileName, 'w') as fasta_file:
                fasta_file.write(parsers.to_fasta([molecule], single_line=True))
                if constraints:
                    fasta_file.write("\n"+constraints)

            if constraints:
                output = commands.getoutput("cd %s ; RNAfold -C < %s"%(self.cache_dir, fileName)).strip()
            elif bp_probabilities:
                commands.getoutput("cd %s ; RNAfold -p < %s"%(self.cache_dir, fileName)).strip()
                with open("%s/%s_dp.ps"%(self.cache_dir, molecule.name), 'r') as ps_file:
                    output = ps_file.read()
            else:
                output = commands.getoutput("cd %s ; RNAfold < %s"%(self.cache_dir, fileName)).strip()
        if raw_output:
            return output
        else:
            if bp_probabilities:
                probabilities = []
                for i in xrange(0, len(molecule.sequence)):
                    probabilities.append(0)
                for line in output.split('\n'):
                    tokens = line.split(' ')
                    if len(tokens) == 4 and tokens[-1] == 'ubox':
                        probabilities[int(tokens[0])-1] = probabilities[int(tokens[0])-1]+float(tokens[2])*float(tokens[2])
                        probabilities[int(tokens[1])-1] = probabilities[int(tokens[1])-1]+float(tokens[2])*float(tokens[2])
                result = []
                for i in xrange(0, len(molecule.sequence)):
                    result.append({'position':i+1, 'residue': molecule.sequence[i], 'probability': probabilities[i] })
                return DataFrame(result, columns = ['position', 'residue', 'probability'])
            else:
                vienna_data = ""
                for line in output.split('\n'):
                    if re.match('^[\.()]+', line):
                        vienna_data += line.split(' ')[0] #we remove the stability value
                    elif not line.startswith("Warning"): #sometimes, the output ends with a line like "Warning from traverse_loop. Loop 1 has crossed regions"
                        vienna_data += line+"\n"
                rnas, base_pairs = parsers.parse_vienna(vienna_data)
                return base_pairs[0]



class Rnainverse(Tool):
    """
    Application Controller for RNAinverse.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("RNAinverse")

    def compute_sequences(self, secondary_structure, molecule, repeats=1):
        """
        Computes sequences folding into a predefined structure.

        Parameters:
        -----------
        - secondary_structure: a secondary structure described as a list of base pairs in a pandas DataFrame
        - molecule: the starting molecule. Any Characters other than "AUGC" will be treated as wild cards and replaced by a random character
        - repeats:

        Returns:
        --------
        an array of RNA molecules

        """

        fileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.bn'
        with open(fileName, 'w') as input_file:
            input_file.write(parsers.to_bn(secondary_structure, len(molecule))+'\n')
            input_file.write('N'*len(molecule))

        rnas = []
        i = 0
        output = commands.getoutput("RNAinverse -R %i < %s"%(repeats, fileName))
        for line in output.split('\n'):
            i+=1
            name = "%s_%i"%(molecule.name, i)
            rnas.append(RNA(name = name, sequence = line.split()[0].strip()))
        return rnas

class RNAMotif(Tool):
    """
    Application Controller for RNAMotif.
    """
    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("rnamotif")

    def search(self, target_molecules, descriptor_file, rna_class = None):
        """
        Launch a search with RNAMotif

        Parameters:
        -----------
        - target_molecules:
        - descriptor_file:
        - rna_class (default: None): optional argument

        Returns:
        --------
        A pandas DataFrame reporting all the SnoReport hits. The column are:
        - source
        - score
        - class (= RNA class)
        - name (= class since unknown)
        - target_positions
        - target_name
        - target_strand ('+' or '-')
        - sequence (the hit primary sequence)
        - bracket_notation as a String
        """
        flag = False
        fasta_file_name = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'
        with open(fasta_file_name, 'w') as fasta_file:
            fasta_file.write(parsers.to_fasta(target_molecules))
        descr_file_name = self.download_file(descriptor_file)
        output = commands.getoutput("rnamotif -descr %s %s | rmprune"%(descr_file_name, fasta_file_name))
        hits = []
        target_molecule = None
        lines = output.split('\n')
        for i in range(0, len(lines)) :
            line = lines[i]
            if line.startswith(self.cache_dir+'/'):
                flag = True
            if line.startswith('>'):
                defline = line[1:-1]
                if lines[i+1].startswith(defline):
                    tokens = re.split('\s+', lines[i+1])
                    if len(tokens) < 6: #the command could produce a segmentation fault. So we need to be sure that the line describing a hit is complete. We need at least 6 tokens to have a "full" line describing the hit. Limitation: the sequence for this hit cound be truncated.
                        break
                    hit = {
                        'source': 'tool:rnamotif:NA',
                        'score': float(tokens[1])
                    }
                    hit['class'] = rna_class
                    hit['name'] = hit['class']
                    for m in target_molecules:
                        if defline == m.name:
                            hit['target_name'] = m.name
                            target_molecule = m
                            break
                    if tokens[2] == '0':
                        hit['target_strand'] = '+'
                        hit['target_positions'] = [int(tokens[3]), int(tokens[3])+int(tokens[4])-1]
                        if target_molecule:
                            hit['sequence'] = target_molecule[hit['target_positions'][0]-1:hit['target_positions'][1]]
                    elif tokens[2] == '1':
                        hit['target_strand'] = '-'
                        hit['target_positions'] = [int(tokens[3])-(int(tokens[4])-1), int(tokens[3])]
                        if target_molecule:
                            hit['sequence'] = target_molecule.get_complement()[hit['target_positions'][0]-1:hit['target_positions'][1]][::-1]
                    length_list = []
                    for sequence in tokens[5:]:
                        length_list.append(len(sequence))
                    bracket_list = []
                    with open(descr_file_name) as infh:
                        line = infh.readline()
                        while line:
                            if line.lstrip().startswith('h5'):
                                bracket_list.append('(')
                            elif line.lstrip().startswith('ss'):
                                bracket_list.append('.')
                            elif line.lstrip().startswith('h3'):
                                bracket_list.append(')')
                            line = infh.readline()
                    chain = ''
                    if len(length_list) == len(bracket_list):
                        j = 0
                        while j != len(bracket_list):
                            chain += bracket_list[j]*length_list[j]
                            j += 1
                    hit['bracket_notation'] = chain
                    hits.append(hit)
                    target_molecule = None
        if not flag:
            raise Exception("No RNAMotif output")

        if len(hits):
            return DataFrame(hits, columns = ['source', 'score', 'class', 'name', 'target_positions', 'target_name', 'target_strand', 'sequence', 'bracket_notation'])
        else:
            return DataFrame()

class Rnaplfold(Tool):
    """
    Application Controller for RNAplfold.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("RNAplfold")

    def scan(self, molecule, winsize, span, width, free_energies = False, raw_output = False):
        """
        Compute the mean probability that regions (a stretch of x consecutive nucleotides) of length 1 to width are unpaired for a sequence. This allows to "scan" sequences for short stable RNA structures.

        Molecule: XXXXXXXXXXXXXXXXXXXXXXXXXXXXX <- the molecule
        winsize:       |-------|                <- the sliding window used to compute the pairing probabilities for each position
        width:            |--|                  <- the sliding window used to compute a score (free energy or probability) to be unpaired.

        Parameters:
        -----------
        - molecule: a Molecule object to scan (see pyrna.features)
        - winsize: window size
        - span: allow only pairs (i,j) with j-i<=span
        - width: length of the regions
        - free_energies (default: False): if True, the method returns free energies, else switch to probabilities
        - raw_output (default: False): if True, the method returns the raw output instead of the pandas Dataframe.

        Returns:
        --------
        a pandas DataFrame reporting the mean probability that regions of length 1 to width are unpaired for a sequence. The index of the Dataframe lists the Molecule positions (from A to molecule's length) and the columns list the length of the region (from 1 to width).
        """

        fileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'
        with open(fileName, 'w') as fasta_file:
            fasta_file.write(parsers.to_fasta([molecule], single_line=True))

        if free_energies:
            commands.getoutput("cd %s ; RNAplfold -W %i -L %i -u %i -O < %s"%(self.cache_dir, winsize, span, width, fileName)).strip()

            h = open('%s/test_openen'%self.cache_dir)
            output = h.read()
            h.close()
        else:
            commands.getoutput("cd %s ; RNAplfold -W %i -L %i -u %i < %s"%(self.cache_dir, winsize, span, width, fileName)).strip()

            with open('%s/test_lunp'%self.cache_dir) as h:
                output = h.read()

        import numpy as np

        if raw_output:
            return output
        else:
            index = []
            columns = None
            data = []
            for line in output.split('\n'):
                line = line.strip()
                if len(line) and not line.startswith("#"):
                    tokens = line.split('\t')
                    index.append(int(tokens[0]))
                    row = {}
                    tokens = tokens[1:]
                    tokens = [float(x) if x != "NA" else np.nan for x in tokens]
                    for pos in range(0,len(tokens)):
                        row[columns[pos]] = tokens[pos]
                    data.append(row)
                elif len(line) and line.startswith("#i$"):
                    columns = [int(x) for x in line.split("l=")[1].split('\t')]
            return DataFrame(data, index = index, columns = columns)

class Rnaplot(Tool):
    """
    Application Controller for RNAplot.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("RNAplot")

    def plot(self, secondary_structure, rna, raw_output = False):
        """
        Parameters:
        -----------
        - secondary_structure: a secondary structure described as a list of base pairs in a pandas DataFrame
        - rna: an RNA object (see pyrna.features)
        - raw_output (default: False): if True, the method returns the raw output instead of the pandas Dataframe.

        Returns:
        --------
        a pandas DataFrame reporting the x and y coordinates for each RNA residue
        """
        output = None
        if self.rest_server:
            _rna = RNA(name="rna", sequence=rna.sequence)
            values = {
                'secondary_structure' : parsers.to_vienna([secondary_structure], [_rna], single_line=True),
                'api_key': self.api_key
            }
            data = urllib.urlencode(values)
            req = urllib2.Request("http://%s/api/computations/rnaplot"%self.rest_server, data)
            response = urllib2.urlopen(req)
            output = str(response.read())
            response.close()
        else:
            _rna = RNA(name="rna", sequence=rna.sequence) #the name of the rna object should not contains any / character.
            vienna_file_name = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'
            with open(vienna_file_name, 'w') as f:
                f.write(parsers.to_vienna([secondary_structure], [_rna], single_line=True))

            commands.getoutput("cd %s ; RNAplot -o svg < %s"%(self.cache_dir, vienna_file_name))

            for f in os.listdir(self.cache_dir):
                if f.endswith('.svg'):
                    with open("%s/%s"%(self.cache_dir, f)) as svg_file:
                        output = svg_file.read()

        if raw_output:
            return output

        import xml.etree.ElementTree as ET

        svg = ET.fromstringlist(output.split('\n'))

        coords = []
        pos = 1
        for text in svg.getchildren()[-1].getchildren()[-1].getchildren():
            coords.append({
                'x': float(text.attrib['x']),
                'y': float(text.attrib['y'])
            })
            pos += 1

        df = DataFrame(coords, index = range(1, len(rna)+1))
        return df

class Rnasubopt(Tool):
    """
    Application Controller for RNAsubopt.
    """

    def __init__(self, cache_dir = "/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, rest_server = rest_server, api_key = api_key)
        if not self.rest_server:
            self.find_executable("RNAsubopt")

    def fold(self, molecule, range = None, random_sample = None):
        """
        Parameters:
        ---------
        - molecule: a Molecule object (see pyrna.features)
        - range (default: None): calculate suboptimal structures within range kcal/mol of the mfe.
        - random_sample (default: None): instead of producing all suboptimals in an energy range, produce a random sample of n suboptimal structures.

        Returns:
        --------
        all the suboptimal secondary structures as a list of pandas DataFrames. Each pandas Dataframe contains a list of base-pairs.
        """
        fileName = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'

        with open(fileName, 'w') as fasta_file:
            fasta_file.write(parsers.to_fasta([molecule], single_line=True))

        output = commands.getoutput("cd %s ; RNAsubopt %s %s < %s"%(self.cache_dir, "-e %i"%range if range else "" ,  "-p %i"%random_sample if random_sample else "", fileName)).strip()
        secondary_structures = []
        for line in output.split('\n'):
            tokens = line.split(' ')
            if not line.startswith('>') and re.match("^[.()]+$", tokens[0]):
                secondary_structures.append(parse_bn(tokens[0]))
        return secondary_structures

class Rnaview(Tool):
    """
    Application Controller for RNAVIEW.
    """
    def __init__(self, cache_dir="/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, api_key = api_key, rest_server = rest_server)
        if not self.rest_server:
            self.find_executable("rnaview")

    def annotate(self, tertiary_structure = None, pdb_content = None, canonical_only = False, raw_output = False):
        """
        Parameters:
        -----------
        - tertiary_structure (default: None): a TertiaryStructure object (see pyrna.features)
        - pdb_content (default: None): the content of a PDB file
        - canonical_only (default: False): if set to True, the helices will be made exclusively with canonical base-pairs: AU c(), GC c() or GU c().
        - raw_output (default: False): if set to True, the method returns the raw RNAML result produced with RNAVIEW .
        """
        if not tertiary_structure and not pdb_content:
            raise Exception("No data provided")
        if self.rest_server:
            values = {
                '3d' : pdb_content if pdb_content else to_pdb(tertiary_structure, export_numbering_system = True),
                'canonical_only': "true" if canonical_only else "false",
                'api_key': self.api_key
            }
            data = urllib.urlencode(values)
            req = urllib2.Request("http://%s/api/computations/rnaview"%self.rest_server, data)
            response = urllib2.urlopen(req)
            xml_content = str(response.read())
            response.close()
        else:
            pdb_file_name = self.cache_dir+'/'+utils.generate_random_name(7)+'.pdb'
            with open(pdb_file_name, 'w') as pdb_file:
                if pdb_content:
                    pdb_file.write(pdb_content)
                else:
                    pdb_file.write(to_pdb(tertiary_structure, export_numbering_system = True))

            commands.getoutput("rnaview -p %s"%(pdb_file_name))

            xml_file_name = pdb_file_name+".xml"
            xml_content = ""
            if os.path.exists(xml_file_name):
                with open(xml_file_name) as xml_file:
                    xml_content = xml_file.read()
            else:
                raise Exception("No file %s"%xml_file_name)
        if raw_output:
            return xml_content
        else:
            import xml.etree.ElementTree as ET
            rnaml_tree = ET.fromstring(xml_content)

            molecule = rnaml_tree.find('molecule')

            rna = RNA(name = tertiary_structure.rna.name, sequence = re.sub('\s+','',molecule.find('sequence').find('seq-data').text))

            secondary_structure = SecondaryStructure(rna)
            secondary_structure.source='tool:rnaview:N.A.'
            new_3D = None

            if len(rna) != len(tertiary_structure.rna): #RNAVIEW can have problems with some residues. Consequently, RNAVIEW produces an RNA molecule with a different sequence. We need to fit the 3D to this molecule.
                new_3D = TertiaryStructure(rna)
                new_3D.source = "tool:rnaview:N.A."
                numbering_system = re.sub('\s{2,}', ' ', molecule.find('sequence').find('numbering-table').text).strip().split(' ')
                #the strategy is the following:
                #- the numbering-table in the XML output stores the labels of the 3D residues used by RNAVIEW
                #- for each residue label, we recover its absolute position in the numbering system of the initial 3D
                residue_absPos = 1
                for residue_label in numbering_system:
                    for absPos, label in tertiary_structure.numbering_system.items():
                        if label == residue_label:
                            new_3D.residues[residue_absPos] = tertiary_structure.residues[int(absPos)]
                            break
                    residue_absPos += 1
            else: #no problem, then we can substitute the RNA of the 2D for the RNA of the 3D
                secondary_structure.rna = tertiary_structure.rna

            if not canonical_only:

                for helix in molecule.find('structure').find('model').find('str-annotation').findall('helix'):
                    secondary_structure.add_helix(helix.get('id'), int(helix.find('base-id-5p').find('base-id').find('position').text), int(helix.find('base-id-3p').find('base-id').find('position').text), int(helix.find('length').text));

                #for single_strand in molecule.find('structure').find('model').find('str-annotation').findall('single-strand'):
                #    end5 = int(single_strand.find('segment').find('base-id-5p').find('base-id').find('position').text)
                #    end3 = int(single_strand.find('segment').find('base-id-3p').find('base-id').find('position').text)
                #    secondary_structure.add_single_strand(single_strand.find('segment').find('seg-name').text, end5, end3-end5+1);

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

                    secondary_structure.add_base_pair(base_pair.find('bond-orientation').text.lower(), edge1, edge2, int(base_pair.find('base-id-5p').find('base-id').find('position').text), int(base_pair.find('base-id-3p').find('base-id').find('position').text));

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

                    orientation = base_pair.find('bond-orientation').text.lower()

                    pos1 = int(base_pair.find('base-id-5p').find('base-id').find('position').text)
                    residue1 = secondary_structure.rna.sequence[pos1-1]
                    pos2 = int(base_pair.find('base-id-3p').find('base-id').find('position').text)
                    residue2 = secondary_structure.rna.sequence[pos2-1]

                    canonical_bps.append([orientation, edge1, edge2, pos1, pos2]) if utils.is_canonical(residue1, residue2, orientation, edge1, edge2) else non_canonical_bps.append([orientation, edge1, edge2, pos1, pos2])

                secondary_structure = base_pairs_to_secondary_structure(secondary_structure.rna, DataFrame(canonical_bps, columns=['orientation', 'edge1', 'edge2', 'pos1', 'pos2']))

                for bp in non_canonical_bps: #the non-canonical interactions are tertiary ones
                    secondary_structure.add_tertiary_interaction(bp[0], bp[1], bp[2], bp[3], bp[4])

            secondary_structure.find_single_strands()

            if new_3D:
                return (secondary_structure, new_3D)
            else:
                return (secondary_structure, tertiary_structure)

class Samtools(Tool):
    """
    Application Controller for Samtools.
    """
    def __init__(self, sam_file, cache_dir="/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, api_key = api_key, rest_server = rest_server)
        self.sam_file = sam_file
        if not self.rest_server:
            self.find_executable("samtools")

    def run(self, user_defined_options = []):
        """
        Generic function to run samtools
        """
        return commands.getoutput("samtools %s"%(" ".join(user_defined_options)))

    def sort_and_index(self):
        """
        Sort and index a SAM file.

        Parameters:
        ----------
        sam_file: the full path of the SAM file

        Returns:
        --------
        the full path of the indexed and sorted BAM file
        """
        path = os.path.realpath(self.sam_file).split('.sam')[0]
        if not os.path.exists("%s.bam"%path):
            self.run(["view", "-bS %s"%self.sam_file, "> %s.bam"%path])
            #commands.getoutput("samtools view -bS %s > %s.bam"%(self.sam_file, path))
            print "BAM file done"
        else:
            print "BAM file already done"
        if not os.path.exists("%s.sorted.bam"%(path)):
            self.run(["sort", "%s.bam"%path, "%s.sorted"%path])
            #commands.getoutput("samtools sort %s.bam %s.sorted"%(path, path))
            print "Sorted BAM file done"
        else:
            print "Sorted BAM file already done"
        if not os.path.exists("%s.sorted.bam.bai"%(path)):
            self.run(["index", "%s.sorted.bam"%path])
            #commands.getoutput("samtools index %s.sorted.bam"%path)
            print "Indexed BAM file done"
        else:
            print "Indexed BAM file already done"
        return "%s.sorted.bam"%path

    def count(self, chromosome_name, start, end, restrict_to_plus_strand = False, restrict_to_minus_strand = False):
        path = os.path.realpath(self.sam_file).split('.sam')[0]
        sorted_bam_file = "%s.sorted.bam"%path
        restrict_to = ""
        if restrict_to_plus_strand:
            restrict_to = '-F 16'
        if  restrict_to_minus_strand:
            restrict_to = '-f 16'
        query = "%s"%chromosome_name
        if start and end:
            query = "%s:%i-%i"%(query, start, end)
        return int(self.run(["view", restrict_to, "-c", sorted_bam_file, query]))

class SnoGPS(Tool):
    """
    Application Controller for SnoGPS.
    """
    def __init__(self, cache_dir="/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, api_key = api_key, rest_server = rest_server)
        if not self.rest_server:
            self.find_executable("pseudoU_test")

    def search(self, target_molecules, scores_table_file, targets_file, descriptor_file):
        """
        Launch a search with SnoGPS

        Parameters:
        ----------
        - target_molecules:
        - scores_table_file:
        - targets_file:
        - descriptor_file:

        Returns:
        --------
        A pandas DataFrame describing the SnoGPS hits. The columns are:
        - source
        - class
        - target_name
        - organism
        - score
        - target_positions
        - target_rRNA
        - name (Rfam name of the snoRNA)
        - target_strand ('+' or '-')
        - sequence (the hit primary sequence)
        - H-box (genomicPositions & sequence)
        - ACA-box (genomicPositions & sequence)
        - L-guide (genomicPositions & sequence)
        - R-guide (genomicPositions & sequence)
        - bracket_notation
        """
        fasta_file_name = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'
        with open(fasta_file_name, 'w') as fasta_file:
            fasta_file.write(parsers.to_fasta(target_molecules))
        scores_file_name = self.download_file(scores_table_file)
        targets_file_name = self.download_file(targets_file)
        descr_file_name = self.download_file(descriptor_file)
        #redirection of the stderr in /dev/null to avoid text insertion (corresponding to a standard error) in the output
        output = commands.getoutput("pseudoU_test -q -L %s -T %s -F /dev/null %s %s 2> /dev/null"%(scores_file_name, targets_file_name, fasta_file_name, descr_file_name))

        hits=[]
        target_molecule = None
        lines = output.split('\n')
        if not len(lines):
            raise Exception("No SnoGPS output")

        for j in range(0, len(lines)):
            line = lines[j]
            if line.startswith('>'):
                hit = {
                    'source': 'tool:snogps:NA',
                    'class' : "Gene, snRNA, snoRNA, HACA-box"
                }
                tokens = re.split('\s+', line)
                for i in range(0, len(tokens)):
                    if tokens[i].startswith('>'):
                        for m in target_molecules:
                            if tokens[i][1:].startswith(m.name):
                                hit["target_name"] = m.name
                                hit["organism"]= m.organism
                                target_molecule = m
                                break
                        hit['score'] = float(tokens[i+1])
                        hit['target_positions'] = [int(tokens[i+2].split('-')[0][1:]), int(tokens[i+2].split('-')[1][:-1])]
                    elif tokens[i].startswith('Cmpl:'):
                        hit['target_rRNA'] = tokens[i+1]
                        hit['name'] = tokens[i+2]
                    elif tokens[i].startswith('Pairs:'):
                        if tokens[i+4].strip() == '(W)':
                            hit['target_strand'] = "+"
                            hit['target_positions'] = [hit['target_positions'][0]+1, hit['target_positions'][1]+1]
                            if target_molecule:
                                hit['sequence'] = target_molecule[hit['target_positions'][0]-1:hit['target_positions'][1]]
                        else:
                            hit['target_strand'] = "-"
                            if target_molecule:
                                hit['sequence'] = target_molecule.get_complement()[hit['target_positions'][0]-1:hit['target_positions'][1]][::-1]
            if target_molecule:
                if line.startswith('#stem1'):
                    chain1 = lines[j+2]
                elif line.startswith('#stem2'):
                    chain2 = lines[j+2]
                    chain0 = ''
                    for i in range(0, len(hit['sequence'])-len(chain1+chain2)):
                        chain0 += ' '
                    chain = chain1+chain0+chain2 if hit['target_strand'] == "+" else chain2+chain0+chain1
                    if len(chain) == len(hit['sequence']):
                        chain = chain.replace('C', 'H')
                        ih = hit['target_positions'][0]+chain.find('H') if hit['target_strand'] is "+" else hit['target_positions'][1]-chain.find('H')
                        ic = hit['target_positions'][0]+chain.rfind('H')-5 if hit['target_strand'] is "+" else hit['target_positions'][1]-chain.rfind('H')+5
                        il = hit['target_positions'][0]+chain.find('L') if hit['target_strand'] is "+" else hit['target_positions'][1]-chain.find('L')
                        jl = hit['target_positions'][0]+chain.rfind('L') if hit['target_strand'] is "+" else hit['target_positions'][1]-chain.rfind('L')
                        ir = hit['target_positions'][0]+chain.find('R') if hit['target_strand'] is "+" else hit['target_positions'][1]-chain.find('R')
                        jr = hit['target_positions'][0]+chain.rfind('R') if hit['target_strand'] is "+" else hit['target_positions'][1]-chain.rfind('R')
                        hit['H-box'] = {'genomicPositions': [ih, ih+5], 'sequence': target_molecule[ih-1:ih+5]} if hit['target_strand'] is "+" else {'genomicPositions': [ih-5, ih], 'sequence': target_molecule.get_complement()[ih-6:ih][::-1]}
                        hit['ACA-box'] = {'genomicPositions': [ic, ic+2], 'sequence': target_molecule[ic-1:ic+2]} if hit['target_strand'] is "+" else {'genomicPositions': [ic-2, ic], 'sequence': target_molecule.get_complement()[ic-3:ic][::-1]}
                        hit['L-guide'] = {'genomicPositions': [il, jl], 'sequence': target_molecule[il-1:jl]} if hit['target_strand'] is "+" else {'genomicPositions': [jl, il], 'sequence': target_molecule.get_complement()[jl-1:il][::-1]}
                        hit['R-guide'] = {'genomicPositions': [ir, jr], 'sequence': target_molecule[ir-1:jr]} if hit['target_strand'] is "+" else {'genomicPositions': [jr, ir], 'sequence': target_molecule.get_complement()[jr-1:ir][::-1]}
                        trantab = maketrans("HLR ", "....")
                        chars = chain.translate(trantab)
                        x = 0
                        i = 0
                        chain = ''
                        for char in chars:
                            if char == 'X' and x < 2 and i == 0:
                                x = 1
                                chain += '('
                            elif char == '.' and 1 <= x <= 2 and i == 0:
                                x = 2
                                chain += '.'
                            elif char == 'I' and i < 2 and x == 2:
                                i = 1
                                chain += '('
                            elif char == '.' and 1 <= i <= 2 and x == 2:
                                i = 2
                                chain += '.'
                            elif char == 'I' and 2 <= i <= 3 and x == 2:
                                i = 3
                                chain += ')'
                            elif char == '.' and i == 3 and x == 2:
                                i = 0
                                chain += '.'
                            elif char == 'X' and 2 <= x <= 3 and i == 0:
                                x = 3
                                chain += ')'
                            elif char == '.' and x == 3 and i == 0:
                                x = 0
                                chain += '.'
                            else:
                                chain += '.'
                        hit['bracket_notation'] = chain
                        hits.append(hit)
                        target_molecule = None
                    else:
                        hits.append(hit)
                        target_molecule = None
                        #raise Exception("SnoRNA hit with different lengths of RNA sequence and structure annotation")
        if hits:
            return DataFrame(hits, columns = ['source', 'class', 'target_name', 'organism', 'score',  'target_positions', 'target_rRNA', 'name', 'target_strand', 'sequence', 'H-box', 'ACA-box', 'L-guide', 'R-guide', 'bracket_notation'])
        else:
            return DataFrame()

class Snoreport(Tool):
    """
    Application Controller for SnoReport.
    """
    def __init__(self, cache_dir="/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, api_key = api_key, rest_server = rest_server)
        if not self.rest_server:
            self.find_executable("snoReport")

    def search(self, molecule, reverse_complement = False):
        """
        Launch a search with SnoReport

        Parameters:
        -----------
        - molecule: a molecule object (see pyrna.features). Its sequence will be used to do the search.
        - reverse_complement (defaut: True): state if the search has do be done on the reverse complement strand

        Returns:
        --------
        A pandas DataFrame describing all the SnoReport hits. The columns are:
        - source
        - score (prediction probability)
        - target_strand ('+' or '-')
        - target_name
        - class ("Gene, snRNA, snoRNA, CD-box" or "Gene, snRNA, snoRNA, HACA-box")
        - name (= class since unknown)
        - target_positions
        - sequence (the hit primary sequence)
        - bracket_notation as a String
        - C-box (genomicPositions & sequence)
        - D-box (genomicPositions & sequence)
        - H-box (genomicPositions & sequence)
        - ACA-box (genomicPositions & sequence)
        """
        fasta_file_name = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'
        with open(fasta_file_name, 'w') as fasta_file:
            fasta_file.write(parsers.to_fasta([molecule], single_line=True))
        output = commands.getoutput("snoReport %s < %s"%("-r" if reverse_complement else "", fasta_file_name))
        hits = []
        lines = output.split('\n')
        if not lines[0].startswith('>'):
            raise Exception("No Snoreport output")
        for i in range(0, len(lines)) :
            line = lines[i]
            if line.startswith('CD') or line.startswith('HACA'):
                words = line.split()
                hit = {
                'source': 'tool:snoreport:NA',
                'score': float(words[2]),
                'target_strand': words[-1],
                'target_name': molecule.name, #or words[-2]
                }
                hit['class'] = "Gene, snRNA, snoRNA, CD-box" if words[0] == 'CD' else "Gene, snRNA, snoRNA, HACA-box"
                hit['name'] = hit['class']
                if not reverse_complement:
                    hit['target_positions'] = [int(words[-4].split(':')[-1])+1, int(words[-3].split(':')[-1])+1]
                    hit['sequence'] = molecule[hit['target_positions'][0]-1:hit['target_positions'][1]]
                else:
                    i2 = len(molecule) - (int(words[-3].split(':')[-1])+1)
                    j2 = len(molecule) - int(words[-4].split(':')[-1])
                    hit['target_positions'] = [i2+1, j2]
                    hit['sequence'] = molecule.get_complement()[i2:j2][::-1]
                #hit['bracket_notation'] = parsers.parse_bn(lines[i-1]) #Panda DataFrame object cannot be encoded by pymongo
                hit['bracket_notation'] = lines[i-1]
                cross_notation = lines[i-2]
                i3 = cross_notation.find('x')
                j3 = cross_notation.rfind('x')
                if line.startswith('CD'):
                    hit['C-box'] = {'genomicPositions': [hit['target_positions'][0]+i3, hit['target_positions'][0]+i3+6], 'sequence': molecule[hit['target_positions'][0]+i3-1:hit['target_positions'][0]+i3+6]} if hit['target_strand'] is "+" else {'genomicPositions': [j2-(i3+6), j2-i3], 'sequence': molecule.get_complement()[j2-(i3+7):j2-i3][::-1]}
                    hit['D-box'] = {'genomicPositions': [hit['target_positions'][0]+j3-3, hit['target_positions'][0]+j3], 'sequence': molecule[hit['target_positions'][0]+j3-4:hit['target_positions'][0]+j3]} if hit['target_strand'] is "+" else {'genomicPositions': [j2-j3, j2-(j3-3)], 'sequence': molecule.get_complement()[j2-j3-1:j2-(j3-3)][::-1]}
                else:
                    hit['H-box'] = {'genomicPositions': [hit['target_positions'][0]+i3, hit['target_positions'][0]+i3+5], 'sequence': molecule[hit['target_positions'][0]+i3-1:hit['target_positions'][0]+i3+5]} if hit['target_strand'] is "+" else {'genomicPositions': [j2-(i3+5), j2-i3], 'sequence': molecule.get_complement()[j2-(i3+6):j2-i3][::-1]}
                    hit['ACA-box'] = {'genomicPositions': [hit['target_positions'][0]+j3-2, hit['target_positions'][0]+j3], 'sequence': molecule[hit['target_positions'][0]+j3-3:hit['target_positions'][0]+j3]} if hit['target_strand'] is "+" else {'genomicPositions': [j2-j3, j2-(j3-2)], 'sequence': molecule.get_complement()[j2-j3-1:j2-(j3-2)][::-1]}
                hits.append(hit)
        if len(hits):
            return DataFrame(hits, columns = ['source', 'score',  'target_strand', 'target_name', 'class', 'name', 'target_positions', 'sequence', 'bracket_notation', 'C-box', 'D-box', 'H-box', 'ACA-box'])
        else:
            return DataFrame()

class Snoscan(Tool):
    """
    Application Controller for SnoScan.
    """
    def __init__(self, cache_dir="/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, api_key = api_key, rest_server = rest_server)
        if not self.rest_server:
            self.find_executable("snoscan")

    def search(self, target_molecules, meth_sites_file, r_rna_file):
        """
        Launch a search with SnoScan

        Parameters:
        ---------
        - target_molecules: a list of Molecule objects (see pyrna.features)
        - meth_sites_file: full path of the methylation file
        - r_rna_file: full path of the file describing the target RNA sequence

        Returns:
        --------
        A pandas DataFrame describing all the SnoScan hits. The columns are:
        - source
        - class
        - target_name
        - organism
        - score
        - target_positions
        - target_strand ('+' or '-')
        - sequence (the hit primary sequence)
        - target_rRNA
        - name (Rfam name of the snoRNA)
        - C-box (genomicPositions & sequence)
        - D-box (genomicPositions & sequence)
        - guide_sequence (genomicPositions & sequence)
        """
        flag = False
        fasta_file_name = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'

        with open(fasta_file_name, 'w') as fasta_file:
            fasta_file.write(parsers.to_fasta(target_molecules))

        meth_file_name = self.download_file(meth_sites_file)
        rrna_file_name = self.download_file(r_rna_file)
        output = commands.getoutput("snoscan -m %s %s %s"%(meth_file_name, rrna_file_name, fasta_file_name))
        hits=[]
        target_molecule = None
        for line in output.split('\n'):
            if line.startswith('Snoscan'):
                flag = True
            if line.startswith('>>'):
                tokens = re.split('\s+', line)
                hit = {'source': 'tool:snoscan:NA', 'class' : "Gene, snRNA, snoRNA, CD-box"}
                for m in target_molecules:
                    if tokens[1].startswith(m.name):
                        hit["target_name"] = m.name
                        hit["organism"]= m.organism
                        target_molecule = m
                        break
                if re.match('^\d\d\.\d\d$', tokens[2]): #bug: for some hits, the score is not a 4-digit float number!
                    hit['score'] = float(tokens[2])
                target_positions = [int(tokens[3].split('-')[0][1:]), int(tokens[3].split('-')[1][:-1])]
                if target_positions[0] < target_positions[1]:
                    hit['target_strand'] = "+"
                    hit['target_positions'] = target_positions
                    if target_molecule:
                        hit['sequence'] = target_molecule[hit['target_positions'][0]-1:hit['target_positions'][1]]
                elif target_positions[0] > target_positions[1]:
                    hit['target_strand'] = "-"
                    hit['target_positions'] = target_positions[::-1]
                    if target_molecule:
                        hit['sequence'] = target_molecule.get_complement()[hit['target_positions'][0]-1:hit['target_positions'][1]][::-1]
                else:
                    raise Exception("Hit with incorrect target positions")
                hit['target_rRNA'] = tokens[5]
                if tokens[6][1:-1] != '-' and tokens[6][1:-1] != 'Undet':
                    hit['name'] = tokens[6][1:-1]
                else:
                    hit['name'] = "Undefined"
            elif target_molecule and 'C-D box dist:' in line:
                pattern = re.compile('\((\d+)-(\d+)\)\s+C-D box dist:\s+(\d+) bp')
                match = pattern.search(line)
                if match:
                    i = int(match.group(1))
                    j = int(match.group(2))
                    dist_cd = int(match.group(3))
                    hit['C-box'] = {'genomicPositions': [i, j], 'sequence': target_molecule[i-1:j]} if hit['target_strand'] is "+" else {'genomicPositions': [j, i], 'sequence': target_molecule.get_complement()[j-1:i][::-1]}
                    hit['D-box'] = {'genomicPositions': [j+dist_cd+1, j+dist_cd+4], 'sequence': target_molecule[j+dist_cd:j+dist_cd+4]} if hit['target_strand'] is "+" else {'genomicPositions': [j-dist_cd-4, j-dist_cd-1], 'sequence': target_molecule.get_complement()[j-dist_cd-5:j-dist_cd-1][::-1]}
            elif target_molecule and line.startswith('Qry seq:'):
                pattern = re.compile('\((\d+)-(\d+)\)')
                match = pattern.search(line)
                if match:
                    i = int(match.group(1))
                    j = int(match.group(2))
                    hit['guide_sequence'] = {'genomicPositions': [j, i], 'sequence': target_molecule[j-1:i]} if hit['target_strand'] is "+" else {'genomicPositions': [i, j], 'sequence': target_molecule.get_complement()[i-1:j][::-1]}
                hits.append(hit)
                target_molecule = None
        if not flag:
            raise Exception("No Snoscan output")
        if len(hits):
            return DataFrame(hits, columns = ['source', 'class', 'target_name', 'organism', 'score', 'target_positions', 'target_strand', 'sequence', 'target_rRNA', 'name', 'C-box', 'D-box', 'guide_sequence'])
        else:
            return DataFrame()

class TrnaScanSE(Tool):
    """
    Application Controller for tRNAscan-SE.
    """
    def __init__(self, cache_dir="/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, api_key = api_key, rest_server = rest_server)
        if not self.rest_server:
            self.find_executable("tRNAscan-SE")

    def scan(self, target_molecules):
        """
        Launch a search with tRNAscan-SE

        Args:
        - target_molecules

        Returns:
        A pandas DataFrame describing all the tRNAscan-SE hits. The column are:
        - source
        - score
        - class (= RNA class)
        - name (=class if unknown)
        - target_positions
        - target_name
        - target_strand ('+' or '-')
        - sequence (the hit primary sequence)
        - tRNA type
        - anticodon
        """
        fasta_file_name = self.cache_dir+'/'+utils.generate_random_name(7)+'.fasta'
        fasta_file = open(fasta_file_name, 'w')
        fasta_file.write(to_fasta(target_molecules))
        fasta_file.close()

        output = commands.getoutput("tRNAscan-SE %s"%fasta_file_name)
        hits=[]
        target_molecule = None
        for target_molecule in target_molecules:
            molecule_name = target_molecule.name.split(" ")[0]
            for line in output.split('\n'):
                if line.startswith(molecule_name):
                    tokens = re.split('\s+', line)
                    hit = {
                        'source': 'tool:tRNAscan-SE:NA',
                        'score': float(tokens[8]),
                        'class' : 'Gene, tRNA',
                        'target_name': molecule_name,
                        'tRNA_type': tokens[4],
                        'anticodon': tokens[5]
                    }
                    hit['name'] = hit['class']
                    target_positions = [int(tokens[2]), int(tokens[3])]
                    if target_positions[0] < target_positions[1]:
                        hit['target_strand'] = "+"
                        hit['target_positions'] = target_positions
                        hit['sequence'] = target_molecule.sequence[hit['target_positions'][0]:hit['target_positions'][1]]
                    elif target_positions[0] > target_positions[1]:
                        hit['target_strand'] = "-"
                        hit['target_positions'] = target_positions[::-1]
                        hit['sequence'] = target_molecule.get_complement()[hit['target_positions'][0]:hit['target_positions'][1]][::-1]
                    else:
                        print "Error: hit with incorrect target positions"
                        print "##########\n" + line + "\n##########"
                    hits.append(hit)
        target_molecule = None
        return DataFrame(hits)

class Tophat(Bowtie2):
    """
    Application Controller for Tophat.
    """
    def __init__(self, cache_dir="/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, api_key = api_key, rest_server = rest_server)
        if not self.rest_server:
            self.find_executable("tophat")

    def align(self, target_molecules, fastq_file, bowtie2_index = None, no_convert_bam = False, no_parsing = True, user_defined_options=[]):
        """
        Align reads to target molecules

        Parameters:
        -----------
        - target_molecules: the genomic sequences to be used for the alignment (an array of Molecule objects, see pyrna.features)
        - fastq_file: the full path for the fastq file containing the reads (as a String)
        - bowtie2_index: the full path of the bowtie2 index. If None, a new index will be build before to do the alignment (as a String)

        Returns:
        --------
        The full path of the SAM file or a pandas DataFrame describing the reads. The columns are:
        - genomicStart (an int)
        - genomicEnd (an int)
        - genomicStrand ('+' or '-')
        - genomeName (a String)
        """

        if not bowtie2_index:
            bowtie2_index = Bowtie2(cache_dir = self.cache_dir).build_index(target_molecules)

        print "tophat %s %s -o %s %s %s"%(' '.join(user_defined_options), "--no-convert-bam" if no_convert_bam else "", self.cache_dir, bowtie2_index, fastq_file)
        commands.getoutput("tophat %s %s -o %s %s %s"%(' '.join(user_defined_options), "--no-convert-bam" if no_convert_bam else "", self.cache_dir, bowtie2_index, fastq_file))

        result_file = None
        if no_convert_bam:
            result_file = self.cache_dir+'/accepted_hits.sam'
        else:
            result_file = self.cache_dir+'/accepted_hits.bam'

        if no_parsing:
            return result_file

        return self.parse_sam(result_file, target_molecules)

class Tophat2(Bowtie2):
    """
    Application Controller for Tophat2.
    """
    def __init__(self, cache_dir="/tmp", rest_server = None, api_key = None):
        Tool.__init__(self, cache_dir = cache_dir, api_key = api_key, rest_server = rest_server)
        if not self.rest_server:
            self.find_executable("tophat2")
            self.find_executable("bed_to_juncs")

    def align(self, target_molecules, fastq_file, bowtie2_index = None, no_convert_bam = False, no_parsing = True, user_defined_options=[]):
        """
        Align reads to target molecules

        Parameters:
        -----------
        - target_molecules: the genomic sequences to be used for the alignment (an array of Molecule objects, see pyrna.features)
        - fastq_file: the full path for the fastq file containing the reads (as a String)
        - bowtie2_index: the full path of the bowtie2 index. If None, a new index will be build before to do the alignment (as a String)

        Returns:
        --------
        The full path of the SAM file or a pandas DataFrame describing the reads. The columns are:
        - genomicStart (an int)
        - genomicEnd (an int)
        - genomicStrand ('+' or '-')
        - genomeName (a String)
        """

        if not bowtie2_index:
            bowtie2_index = Bowtie2(cache_dir = self.cache_dir).build_index(target_molecules)
        elif not os.path.exists(bowtie2_index+".1.bt2"):
            bowtie2_index = Bowtie2(cache_dir = self.cache_dir, index_path = bowtie2_index).build_index(target_molecules)

        print "tophat2 %s %s -o %s %s %s"%(' '.join(user_defined_options), "--no-convert-bam" if no_convert_bam else "", self.cache_dir, bowtie2_index, fastq_file)
        commands.getoutput("tophat2 %s %s -o %s %s %s"%(' '.join(user_defined_options), "--no-convert-bam" if no_convert_bam else "", self.cache_dir, bowtie2_index, fastq_file))

        result_file = None
        if no_convert_bam:
            result_file = self.cache_dir+'/accepted_hits.sam'
        else:
            result_file = self.cache_dir+'/accepted_hits.bam'

        if no_parsing:
            return result_file

        return self.parse_sam(result_file, target_molecules)

    def bed_to_juncs(self, junc_file):
        """
        This function converts the file junctions.bed into the file junc_file given as parameter. This junc_file can be provided to tophat2 with the option -j.

        Parameters:
        -----------
        - junc_file: the name of the output file
        """
        lines = open(self.cache_dir+"/junctions.bed").readlines()
        open(self.cache_dir+"/junctions_without_header.bed", 'w').writelines(lines[1:]) #the first line seems to crash the process....
        commands.getoutput("bed_to_juncs < %s > %s"%(self.cache_dir+"/junctions_without_header.bed", junc_file))
        print "bed_to_juncs < %s > %s"%(self.cache_dir+"/junctions_without_header.bed", junc_file)
