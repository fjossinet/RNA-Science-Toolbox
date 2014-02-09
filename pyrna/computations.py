import os, commands, re, shutil, sys, urllib
from string import maketrans
from pandas import DataFrame
import parsers, utils
from features import RNA, SecondaryStructure, TertiaryStructure
from parsers import base_pairs_to_secondary_structure, parse_bracket_notation, to_fasta
from distutils.spawn import find_executable

class Tool:
    def __init__(self, working_dir):
        self.working_dir = working_dir

    def find_executable(self, executable):
        if not find_executable(executable):
            raise Exception("%s is not available in your PATH"%executable)   

    def download_file(self, uri):
        """
        Download a file from a URI (local file or Web address).

        Args:
        - uri: an URI. Must start with 'file://', 'http://' or 'https://'"

        Returns:
        the local path of the file downloaded.
        """
        path = self.working_dir+"/"+utils.generate_random_name(7)
        os.mkdir(path)
        if not uri.endswith("/"):
            file_name = path+'/'+uri.split('/')[-1]
        else:
            print "Error: your URI should not end with the character '/'"
            sys.exit(1)
        if uri.startswith("https://") or uri.startswith("http://"):
            urllib.urlretrieve(uri, file_name)
        elif uri.startswith("file://"):
            shutil.copy(uri.split('file://')[1], self.working_dir)
        else:
            print "Error: your URI must start with 'file://', 'http://' or 'https://'" #if the given path start with /Users/, URI must start with file:///Users/ 
            sys.exit(1)
        return file_name

class Bcheck(Tool):
    """
    Application Controller for Bcheck.
    """

    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("Bcheck")

    def check(self, target_molecules):
        """
        Scan the target molecules for RnaseP hits 

        Args:
        - target_molecules: the target molecules (as a list of Molecule objects, see pyrna.features)

        Returns:
        A pandas DataFrame describing all the hits. The columns are:
        - e_value
        - score 
        - target_strand ('+' or '-')
        - target_positions
        - target_name
        - sequence (the hit primary sequence)
        - organism
        """
        fileName = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'

        f = open(fileName, 'w');
        f.write(parsers.to_fasta(target_molecules))
        f.close()
        
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
    def __init__(self, target_molecules, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("formatdb")
        self.find_executable("blastall")
        self.target_molecules = target_molecules

    def format_db(self, is_nucleotide=True):
        """
        This function saves the target_molecules in a fasta file and format them for blast search.
        It stores the name and location of the fasta file generated.
        """
        path = self.working_dir+"/"+utils.generate_random_name(7)
        os.mkdir(path)
        fasta_file = open("%s/input.fasta"%path, 'w+b')
        fasta_file.write(parsers.to_fasta(self.target_molecules))
        fasta_file.close()
        commands.getoutput("cd %s ; formatdb -i %s -p %s -o"%(self.working_dir, fasta_file.name, "F" if is_nucleotide else "T"))
        self.formatted_db = fasta_file.name

    def parse_output(self, output):
        """
        Parse the output of the blast tool as a String

        Args:
        - output_file: the blast output file

        Returns:
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
        Return all the hits in a pandas DataFrame

        Args:
        - query_molecule: a Molecule object (see pyrna.features).

        Returns:
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
        query_file = open("%s/query.fasta"%tmp_dir, 'w+b')
        query_file.write(query_molecule.to_fasta())
        query_file.close()
        return self.parse_output(commands.getoutput("cd %s ; blastall -p blastn -d %s -i %s"%(self.working_dir, self.formatted_db, query_file.name)))

    def rpsblast(self):
        pass

class Blastr(Blast):

    """
    Application Controller for blastR.
    """

    def __init__(self, target_molecules, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("formatdbR.pl")
        self.find_executable("blastallR.pl")
        self.target_molecules = target_molecules

    def format_db(self):
        path = self.working_dir+"/"+utils.generate_random_name(7)
        os.mkdir(path)
        fasta_file = open("%s/input.fasta"%path, 'w+b')
        fasta_file.write(parsers.to_fasta(self.target_molecules))
        fasta_file.close()
        commands.getoutput("cd %s ; formatdbR.pl -i %s"%(self.working_dir, fasta_file.name))
        self.formatted_db = fasta_file.name

    def blastallr(self, query_molecule):
        tmp_dir = os.path.dirname(self.formatted_db)
        query_file = open("%s/query.fasta"%tmp_dir, 'w+b')
        query_file.write(query_molecule.to_fasta())
        query_file.close()
        return self.parse_output(commands.getoutput("cd %s ; blastallR.pl -p blastr -i %s -d %s"%(self.working_dir, query_file.name, self.formatted_db)))

class Bowtie(Tool):
    """
    Application Controller for Bowtie.
    """
    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("bowtie-build")
        self.find_executable("bowtie")

    def align(self, target_molecules, fastq_file, threshold = 2, fill_cluster_with_genomic_annotations = False):
        """
        Launch an alignment of reads with Bowtie

        Args:
        - target_molecules: the genomic sequences to be used for the alignment (an array of Molecule objects, see pyrna.features) 
        - fastq_file: the absolute path for the fastq file containing the reads (as a String)
        - threshold: this methods will only return the clusters whose size is >= threshold (default: 1)
        - fill_cluster_with_genomic_annotations: if True, all the genomic annotations making the cluster are stored in a list linked to the key 'genomic_annotations' (default: False).

        Returns:
        A pandas DataFrame describing all the cluster of reads hits. The column are:
        - genomicStart (an int)
        - genomicEnd (an int)
        - annotations_count (an int)
        - genomic_annotations (a list. Is empty if the argument fill_cluster_with_genomic_annotations is False)
        - genomeName (a String)
        """

        fastq_file = os.path.abspath(os.path.normpath(fastq_file))        
        fasta_file_name = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'
        fasta_file = open(fasta_file_name, 'w')
        fasta_file.write(to_fasta(target_molecules))
        fasta_file.close()
        
        db_path = self.working_dir+"/bowtie_db_"+utils.generate_random_name(7)
        result_file = self.working_dir+'/'+utils.generate_random_name(7)+'.sam'
        commands.getoutput("bowtie-build %s %s"%(fasta_file_name, db_path))
        commands.getoutput("bowtie %s -q \"%s\" -S %s"%(db_path, fastq_file, result_file))
        print "SAM file %s produced successfully!!"%result_file

        from parsers import parse_sam
        from pyrna.utils import cluster_genomic_annotations

        aligned_reads, total_reads, tids  = parse_sam(result_file)
        all_clustered_reads = []
        total_aligned_reads = 0
        for aligned_reads_per_genome in aligned_reads:
            total_aligned_reads += len(aligned_reads_per_genome)
            if len(aligned_reads_per_genome):
                clustered_reads = cluster_genomic_annotations(aligned_reads_per_genome, threshold, fill_cluster_with_genomic_annotations) #we keep only the clusters with at least "threshold" reads
                genomeName = tids[aligned_reads_per_genome[0]['tid']]
                for clustered_read in clustered_reads:
                    clustered_read['genomeName'] = genomeName      
                all_clustered_reads += clustered_reads
        print "%i reads found, %i reads aligned..."%(total_reads, total_aligned_reads)
        return DataFrame(all_clustered_reads)

class Clustalw(Tool):

    """
    Application Controller for clustalw.
    """

    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("clustalw2")

    def align(self, molecules):
        """
        Return the clustalw output as a String and an array of aligned molecules
        """
        fileName = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'

        f = open(fileName, 'w');
        f.write(parsers.to_fasta(molecules))
        f.close()
        
        commands.getoutput("clustalw2 -infile=%s"%fileName) 

        outputFile = open(fileName.split('.fasta')[0]+".aln", 'r')
        data = outputFile.read()
        outputFile.close()

        return data, self.parse_output(data, molecules)

    def parse_output(self, output, molecules):
        """
        Return an array of aligned molecules

        Args:
        - output: the clustalw output as a String
        - molecules: the list of molecules that have been used to produce the output
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
    def __init__(self, local_mode = True, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("cmalign")
        self.local_mode = local_mode

    def align(self, molecules, rfam_id = None, rfam = None, stockholm_content = None, cm_content = None):
        """
        Align new ncRNA candidates (the molecules object argument) to the RFAM family (defined by the rfam_id argument) using the CM model stored in cm_file.

        Args:
        - molecules: an array of molecules to align
        - rfam_id: the id of the RFAM family to use for the alignment
        - rfam: an Rfam object

        Returns:
        a tuple like: (list of all the aligned molecules, dict of organism names (keys) and accession numbers/start-end (values), Dataframe of the consensus 2D)
        """
        path = self.working_dir+"/"+utils.generate_random_name(7)
        os.mkdir(path)
        fasta_file = open("%s/input.fasta"%path, 'w+b')
        fasta_file.write(parsers.to_fasta(molecules))
        fasta_file.close()
        if not stockholm_content:
            try:
                stockholm_content = rfam.getEntry(rfam_id, format='stockholm')
            except Exception, e:
                raise e
        if rfam_id:
            stockholm_file = open("%s/%s.stk"%(path,rfam_id), 'w+b')
        else:
            stockholm_file = open("%s/%s.stk"%(path, utils.generate_random_name(7)), 'w+b')
        
        stockholm_file.write(stockholm_content)
        stockholm_file.close()
        
        output = None

        cm_file = None

        if cm_content:
            cm_file = self.working_dir+'/'+utils.generate_random_name(7)+'.cm'

            f = open(cm_file, 'w')
            f.write(cm_content)
            f.close()
        
        if self.local_mode:
            output = commands.getoutput("cmalign -l --withali %s %s %s"%(stockholm_file.name, cm_file if cm_file else rfam.rootDir+'/CMs/'+rfam_id+".cm", fasta_file.name))
        else:
            output = commands.getoutput("cmalign --withali %s %s %s"%(stockholm_file.name, cm_file if cm_file else rfam.rootDir+'/CMs/'+rfam_id+".cm", fasta_file.name))
        shutil.rmtree(path)
        return parsers.parse_stockholm("#=GF AC "+rfam_id+"\n"+output) #we add the accession number to the output of cmalign to be able to use the stockholm parser which need this information

class Cmbuild(Tool):
    """
    Application Controller for Cmbuild.
    """

    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("cmbuild")

    def build(self, stockholm_content):
        """
        Return the content of a covariance model (as a String)
        """
        name = utils.generate_random_name(7)
        stockholm_file = self.working_dir+'/'+name+'.sto'
        cm_file = self.working_dir+'/'+name+'.cm'

        f = open(stockholm_file, 'w');
        f.write(stockholm_content)
        f.close()
        
        commands.getoutput("cmbuild %s %s "%(cm_file, stockholm_file))

        cm_content = None
        f = open(cm_file, 'r')
        cm_content = f.read()
        f.close()

        return cm_content

class Cmcalibrate(Tool):
    """
    Application Controller for Cmcalibrate.
    """

    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("cmcalibrate")

    def calibrate(self, cm_content):
        """
        Return the calibrated content of a covariance model (as a String)
        """
        name = utils.generate_random_name(7)
        cm_file = self.working_dir+'/'+name+'.cm'

        f = open(cm_file, 'w')
        f.write(cm_content)
        f.close()

        commands.getoutput("cmcalibrate %s"%cm_file)
        
        f = open(cm_file, 'r')
        cm_content = f.read()
        f.close()

        return cm_content   

class Cmsearch(Tool):
    """
    Application Controller for Cmsearch.
    """

    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("cmsearch")

    def search(self, molecules, rfam_id = None, rfam = None, cm_content = None,  gathering_threshold = True):
        """
        Launch a search with cmsearch

        Args:
        - molecules: the molecules used to do the search
        - rfam_id : to id of the RFAM family
        - rfam: an Rfam object (see pyrna.db)
        - cm_content: the content of a CM file as a String (default is None).

        Returns:
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
            fileName = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'

            f = open(fileName, 'w');
            f.write(parsers.to_fasta(molecules))
            f.close()
            
            cmFile = self.working_dir+'/'+utils.generate_random_name(7)+'.cm'

            f = open(cmFile, 'w');
            f.write(cm_content)
            f.close()

            if not gathering_threshold:
                return self.parse_output(commands.getoutput("cmsearch "+cmFile+" "+fileName), molecules, False)
            else:
                return self.parse_output(commands.getoutput("cmsearch --ga "+cmFile+" "+fileName), molecules, True)

        else:
            #write molecules as FASTA
            fileName = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'

            f = open(fileName, 'w');
            f.write(parsers.to_fasta(molecules))
            f.close()

            if not gathering_threshold:
                return self.parse_output(commands.getoutput("cmsearch "+rfam.rootDir+"/CMs/"+rfam_id+".cm "+fileName), molecules, False)
            else:
                return self.parse_output(commands.getoutput("cmsearch --ga "+rfam.rootDir+"/CMs/"+rfam_id+".cm "+fileName), molecules, True)

    def parse_output(self, output, molecules, gathering_threshold = True):
        """
        Parse the output of the cmsearch tool.

        Args:
        - output: the cmsearch output as a String
        - molecules: the molecules used to do the search

        Returns:
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
                        if lines[i].strip().endswith("%i"%target_end) or len(lines[i+2].strip()) == 0:
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
    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("contrafold")

    def fold(self, molecule):
        """
        Return the secondary structure as a list of base-pairs in a pandas DataFrame
        """
        fileName = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'
        fasta_file = open(fileName, 'w')
        fasta_file.write(parsers.to_fasta([molecule]))
        fasta_file.close()
        output = commands.getoutput("cd %s ; contrafold predict %s"%(self.working_dir, fileName)).strip()
        return parsers.parse_vienna(output)[0][1]


class Gotohscan(Tool):
    """
    Application Controller for GotohScan.
    """

    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("GotohScan2a") 

    def scan(self, query_molecule, target_molecules, evalue = 1e-3):
        """
        Scan the target molecules for hits of query molecule
        Args:
        - query_molecule: the query molecule (as a Molecule object)
        - target_molecules: the target molecules (as a list of Molecule objects)
        - evalue: only hits with a evalue better or equal to this value will be returned (default is 1e-3)

        Returns:
        A pandas DataFrame describing all the  hits. The columns are:
        - e_value
        - score
        - target_strand ('+' or '-')
        - target_positions
        - target_name
        - sequence (the hit primary sequence)
        - organism
        """

        hits = []
        queryFileName = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'
        query_file = open(queryFileName, 'w')
        query_file.write(parsers.to_fasta([query_molecule]))
        query_file.close()

        targetFileName = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'
        target_file = open(targetFileName, 'w')
        target_file.write(parsers.to_fasta(target_molecules))
        target_file.close()

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

    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("mlocarna") 

    def align(self, molecules):
        """
        Returns a tuple like (list of aligned molecules, secondary structure computed as a list of base-pairs in a pandas DataFrame)
        """
        fileName = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'
        fasta_file = open(fileName, 'w')
        fasta_file.write(parsers.to_fasta(molecules))
        fasta_file.close()
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
                consensus2D = parsers.parse_bracket_notation(tokens[1])

        rnas = []
        for k,v in aligned_molecules.iteritems():
            rnas.append(RNA(name=k, sequence=v))

        return (rnas, consensus2D)

class RnaAlifold(Tool):

    """
    Application Controller for RNAalifold.
    """

    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("RNAalifold") 

    def align(self, alignment):
        """
        Returns the mfe structure in bracket notation as a String
        """
        fileName = self.working_dir+'/'+utils.generate_random_name(7)+'.aln'
        aln_file = open(fileName, 'w')
        aln_file.write(alignment)
        aln_file.close()
        return commands.getoutput("cd %s ; RNAalifold < %s"%(self.working_dir, fileName)).strip().split('\n')[-1].split(' ')[0]

class Rnafold(Tool):

    """
    Application Controller for RNAfold.
    """

    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("RNAfold") 

    def fold(self, molecule, constraints = None):
        """
        Returns the secondary structure as a list of base-pairs in a pandas DataFrame
        """
        fileName = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'
        fasta_file = open(fileName, 'w')
        fasta_file.write(parsers.to_fasta([molecule], single_line=True))
        if constraints:
            fasta_file.write("\n"+constraints)
        fasta_file.close()
        if constraints:
            output = commands.getoutput("cd %s ; RNAfold -noLP -C < %s"%(self.working_dir, fileName)).strip()
        else:
            output = commands.getoutput("cd %s ; RNAfold -noLP < %s"%(self.working_dir, fileName)).strip()
        vienna_data = "" 
        for line in output.split('\n'):
            if re.match('^[\.()]+', line):
                vienna_data += line.split(' ')[0] #we remove the stability value
            elif not line.startswith("Warning"): #sometimes, the output ends with a line like "Warning from traverse_loop. Loop 1 has crossed regions"
                vienna_data += line+"\n"
        return parsers.parse_vienna(vienna_data)[0][1] 

class Rnainverse(Tool):
    """
    Application Controller for RNAinverse.
    """
    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("RNAinverse")

    def compute_sequences(self, secondary_structure, molecule, repeats=1):
        """
        Computes sequences folding into a predefined structure.

        Args:
        - secondary_structure: a secondary structure described as a list of base pairs in a pandas data frame
        - molecule: the starting molecule. Any Characters other than "AUGC" will be treated as wild cards and replaced by a random character
        - repeats:

        Returns: an array of RNA molecules

        """

        fileName = self.working_dir+'/'+utils.generate_random_name(7)+'.bn'
        input_file = open(fileName, 'w')
        input_file.write(parsers.to_bn(secondary_structure, len(molecule))+'\n')
        input_file.write('N'*len(molecule))
        input_file.close()

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
    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("rnamotif")

    def search(self, target_molecules, descriptor_file, rna_class = None):
        """
        Launch a search with RNAMotif

        Args:
        - target_molecules 
        - descriptor_file
        - rna_class: optional argument 

        Returns:
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
        fasta_file_name = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'
        fasta_file = open(fasta_file_name, 'w')
        fasta_file.write(parsers.to_fasta(target_molecules))
        fasta_file.close()
        descr_file_name = self.download_file(descriptor_file)
        output = commands.getoutput("rnamotif -descr %s %s | rmprune"%(descr_file_name, fasta_file_name))
        hits = []
        target_molecule = None
        lines = output.split('\n')
        for i in range(0, len(lines)) :
            line = lines[i]
            if line.startswith(self.working_dir+'/'):
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
                    infh = open(descr_file_name, 'r')
                    line = infh.readline()
                    while line:
                        if line.lstrip().startswith('h5'):
                            bracket_list.append('(')
                        elif line.lstrip().startswith('ss'):
                            bracket_list.append('.')
                        elif line.lstrip().startswith('h3'):
                            bracket_list.append(')')
                        line = infh.readline()
                    infh.close()
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
    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("RNAplfold")

    def fold(self, molecule, winsize, span):
        """
        Args:
        - molecule: an Molecule object (see pyrna.features)
        - winsize: window size
        - span: allow only pairs (i,j) with j-i<=span

        Returns:
        a pandas data frame reporting the average pair probabilities over windows of size winsize.
        """

        fileName = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'
        fasta_file = open(fileName, 'w')
        fasta_file.write(parsers.to_fasta([molecule], single_line=True))
        fasta_file.write("\n"+constraints)
        fasta_file.close()
        output = commands.getoutput("cd %s ; RNAplfold -W %i -L %i < %s"%(self.working_dir, winsize, span, fileName)).strip()
        print output

class Rnaplot(Tool):
    """
    Application Controller for RNAplot.
    """
    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("RNAplot")

    def plot(self, secondary_structure, rna):
        """
        Args:
        - secondary_structure: a secondary structure described as a list of base pairs in a pandas data frame
        - rna: an RNA object (see pyrna.features)

        Returns:
        a pandas data frame reporting the x and y coordinates for each rna position
        """
        path=self.working_dir+"/"+utils.generate_random_name(7)
        os.mkdir(path)

        _rna = RNA(name="toto", sequence=rna.sequence) #the name of the rna object should not contains any / character.
        vienna_file_name = path+'/'+utils.generate_random_name(7)+'.fasta'
        f = open(vienna_file_name, 'w')
        f.write(parsers.to_vienna(secondary_structure, [_rna], single_line=True))
        f.close()

        import xml.etree.ElementTree as ET

        print "cd %s ; RNAplot -o svg < %s"%(path, vienna_file_name)

        commands.getoutput("cd %s ; RNAplot -o svg < %s"%(path, vienna_file_name))

        for f in os.listdir(path):
            if f.endswith('.svg'):
                svg_file = open("%s/%s"%(path, f), 'r')
                svg = ET.fromstringlist(svg_file.readlines())
                svg_file.close()

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
    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("RNAsubopt")

    def fold(self, molecule, range = None, random_sample = None):
        """
        Returns all the suboptimal secondary structures as a list of pandas DataFrames. Each pandas Dataframe contains a list of base-pairs.
        """
        fileName = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'
        fasta_file = open(fileName, 'w')
        fasta_file.write(parsers.to_fasta([molecule], single_line=True))
        fasta_file.close()
        output = commands.getoutput("cd %s ; RNAsubopt -noLP %s %s < %s"%(self.working_dir, "-e %i"%range if range else "" ,  "-p %i"%random_sample if random_sample else "", fileName)).strip()
        secondary_structures = []
        for line in output.split('\n'):
            tokens = line.split(' ')
            if not line.startswith('>') and re.match("^[.()]+$", tokens[0]):
                secondary_structures.append(parse_bracket_notation(tokens[0]))
        return secondary_structures

class Rnaview(Tool):
    """
    Application Controller for RNAVIEW.
    """
    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("rnaview")

    def annotate(self, tertiary_structure, canonical_only = False, raw_output = False):
        """
        Args:
        - tertiary_structure: a TertiaryStructure object (see pyrna.features). Could be raw PDB content (if raw_output is True)
        - canonical_only: if True, the helices will be made with canonical base-pairs only: meaning AU c(), GC c() or GU c(). Default value is False.
        - raw_output: if True, the method returns the raw RNAML result produced with RNAVIEW (default is False). 
        """
        if raw_output:
            pdb_file_name = self.working_dir+'/'+utils.generate_random_name(7)+'.pdb'
            pdb_file = open(pdb_file_name, 'w')
            pdb_file.write(tertiary_structure)
            pdb_file.close()

            commands.getoutput("rnaview -p %s"%(pdb_file_name))
            
            xml_file_name = pdb_file_name+".xml"
            xml_content = ""
            if os.path.exists(xml_file_name):
                xml_file = open(xml_file_name, 'r')
                xml_content = xml_file.read()
                xml_file.close()
                return xml_content
            else:
                raise Exception("No file %s"%xml_file_name)
        else:
            pdb_file_name = self.working_dir+'/'+utils.generate_random_name(7)+'.pdb'
            pdb_file = open(pdb_file_name, 'w')
            pdb_file.write(parsers.to_pdb(tertiary_structure, export_numbering_system = True))
            pdb_file.close()

            commands.getoutput("rnaview -p %s"%(pdb_file_name))
            
            xml_file_name = pdb_file_name+".xml"

            if os.path.exists(xml_file_name):
                
                xml_file = open(xml_file_name, 'r')
                xml_content = xml_file.read()
                xml_file.close()

                #print xml_content
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
                    previous_absolute_position = 0
                    missing_residues = 0
                    #the strategy is the following:
                    #- the numbering-table in the XML output stores the labels of the 3D residues used by RNAVIEW
                    #- for each residue label, we recover its absolute position in the numbering system of the initial 3D
                    #- we're keeping track of the number of missing residues
                    #- in the new 3D, the new absolute position = the original absolute position - the missing residues
                    #- in the new 3D, the numbering-system link the residue label to its new absolute position
                    print tertiary_structure.numbering_system
                    for residue_label in numbering_system:
                        absolute_position = int(tertiary_structure.numbering_system[residue_label]) #we get the absolute position in the initial 3D according to the residue label in the numbering table
                        missing_residues += (absolute_position-previous_absolute_position-1)
                        new_3D.residues[absolute_position-missing_residues] = tertiary_structure.residues[absolute_position] #in the new tertiary structure, the new absPos is the previous one minus the missing residues
                        new_3D.numbering_system[residue_label] = absolute_position-missing_residues
                        previous_absolute_position = absolute_position
                else: #no problem, then we can substitute the RNA of the 2D for the RNA of the 3D 
                    secondary_structure.rna = tertiary_structure.rna

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
                
                if new_3D:
                    return (secondary_structure, new_3D)
                else:
                    return (secondary_structure, tertiary_structure)
            else:
                raise Exception("No file %s"%xml_file_name) 

class SnoGPS(Tool):
    """
    Application Controller for SnoGPS.
    """
    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("pseudoU_test")

    def search(self, target_molecules, scores_table_file, targets_file, descriptor_file):
        """
        Launch a search with SnoGPS

        Returns:
        A pandas DataFrame describing all the SnoGPS hits. The columns are:
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
        fasta_file_name = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'
        fasta_file = open(fasta_file_name, 'w')
        fasta_file.write(parsers.to_fasta(target_molecules))
        fasta_file.close()
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
    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("snoReport")

    def search(self, molecule, reverse_complement = False):
        """
        Launch a search with SnoReport

        Args:
        - molecule: a molecule object (see pyrna.features). Its sequence will be used to do the search.
        - reverse_complement: do the search on the reverse complement strand (defaut is True)

        Returns:
        A pandas DataFrame describing all the SnoReport hits. The column are:
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
        fasta_file_name = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'
        fasta_file = open(fasta_file_name, 'w')
        fasta_file.write(parsers.to_fasta([molecule], single_line=True))
        fasta_file.close()
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
                #hit['bracket_notation'] = parsers.parse_bracket_notation(lines[i-1]) #Panda DataFrame object cannot be encoded by pymongo 
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
    def __init__(self, working_dir="/tmp"):
        Tool.__init__(self, working_dir)
        self.find_executable("snoscan")

    def search(self, target_molecules, meth_sites_file, r_rna_file):
        """
        Launch a search with SnoScan

        Returns:
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
        fasta_file_name = self.working_dir+'/'+utils.generate_random_name(7)+'.fasta'
        fasta_file = open(fasta_file_name, 'w')
        fasta_file.write(parsers.to_fasta(target_molecules))
        fasta_file.close()
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