from pandas import DataFrame
import re
from bson.objectid import ObjectId
from itertools import groupby
from operator import itemgetter

class Block:
    """
    A continuous range of molecular positions, with a single start and end point. 
    """
    def __init__(self, start, end):
        if start < end:
            self.start = start
            self.end = end
        else:
            self.start = end
            self.end = start

    def is_before(self, block):
        pass

    def is_beside(self, block):
        pass

    def intersects(self, block):
        pass

    def merge(self, block):
        pass

class Location:
    """
    A Location defines a range of molecular positions, continuous or not. A location is made with Block objects.
    """
    def __init__(self, start = None, end = None, single_positions = None, nested_lists = None):
        """
        To instantiate a Location, you can:
        - set a start and end position: Location(start=34, end=69). The location will contain all the positions between the start end end.
        - list all the single positions (sorted or not) to be contained in the location: Location(single_positions=[34, 56, 57, 58, 67, 68, 69])
        - list the ranges of continuous positions as nested lists: Location(nested_lists=[[34,34], [56,58], [67,69]]) 
        """       
        self.blocks = []
        if start and end:
            self.add_block(Block(start, end))
        elif single_positions:
            single_positions.sort()
            for k, g in groupby(enumerate(single_positions), lambda (i,x):i-x):
                _range = map(itemgetter(1), g)
                self.blocks.append(Block(min(_range), max(_range)))
        elif nested_lists:
            for nested_list in nested_lists:
                self.blocks.append(Block(min(nested_list), max(nested_list)))

    def add_block(self, block):
        blocks_to_remove = []
        
        for _block in self.blocks:
            if block.is_before(_block) and not block.is_beside(_block):
                break
            elif block.intersects(_block) or block.is_beside(_block):
                block.merge(_block)
                blocks_to_remove.append(_block)
                #its necessary to continue to see if the new Block can merge with other blocks
                continue
            elif len(blocks_to_remove):
                break
        
        for block_to_remove in blocks_to_remove:
            self.blocks.remove(block_to_remove)
        
        self.blocks.append(block)
        self.blocks = sorted(self.blocks, key=lambda block: block.start)

    def remove_location(self, location):
        """
        Return a new Location object from the difference between the current Location and the Location given as argument.
        Difference means all the positions not found in the Location given as argument
        """
        single_positions_1 = self.get_single_positions()
        single_positions_2 = location.get_single_positions()

        diff = list(set(single_positions_1) - set(single_positions_2))

        return Location(single_positions = diff)

    def remove_locations(self, locations):
        """
        Return a new Location object from the difference between the current Location with all the Locations given in a list as argument.
        Difference means all the positions not found in the Locations given as argument
        """
        single_positions_1 = self.get_single_positions()
        single_positions_2 = []

        for location in locations:
            single_positions_2 += location.get_single_positions()

        diff = list(set(single_positions_1) - set(single_positions_2))

        return Location(single_positions = diff)


    def get_single_positions(self):
        """
        Returns:
        ------
        all the single positions making this Location as a list.
        """
        single_positions = []
        for block in self.blocks:
            single_positions += xrange(block.start, block.end+1)
        return single_positions

    def has_position(self, position):
        """
        Test if the location encloses a single position. 
        Parameters:
        ---------
        position: an integer
        """
        return position in self.get_single_positions()

    def start(self):
        return self.blocks[0].start

    def end(self):
        return self.blocks[-1].end        


class Molecule:
    def __init__(self, name):
        self._id = str(ObjectId())
        self.modified_residues = []
        self.name = name
        self.family = None
        self.organism = None
        self.lineage = None
        self.source = 'N.A.:N.A.:N.A.'
        self.sequence = ""

    def get_gaps_positions(self):
        positions = []
        i = 0
        for c in list(self.sequence):
            if c == '-':
                positions.append(i)
            i += 1
        return positions

    def to_fasta(self, single_line=False):
        lines = []
        lines.append(">" + self.name)
        if single_line:
            lines.append(self.sequence)
        else:
            c = 0
            while c < len(self.sequence):
                d = min(len(self.sequence), c + 79)
                lines.append(self.sequence[c:d])
                c += 79
        return '\n'.join(lines)

    def _repr_html_(self):
        from pyrna.utils import chunks
        subsequences = [''.join(subsequence) for subsequence in chunks(list(self.sequence),60)] 
        html = "<pre>"
        i = 0
        for subsequence in subsequences:
            html += str(i*60+1)+"\t"+re.sub('U', '<font color="green">U</font>', re.sub('T', '<font color="green">T</font>', re.sub('C', '<font color="orange">C</font>', re.sub('G', '<font color="red">G</font>', re.sub('A', '<font color="blue">A</font>',subsequence)))))+"\n"
            i += 1
        html += "</pre>"
        return html

    def __add__(self, seq):
        if seq.__class__ == str:
            self.sequence = ''.join([self.sequence, seq])

    def __sub__(self, length):
        if length.__class__ == int and length <= len(self.sequence):
            self.sequence = self.sequence[0: len(self.sequence)-length]

    def __len__(self):
        return len(self.sequence)

    def __iter__(self):
        return iter(self.sequence)

    def __getslice__(self, i, j):
        return self.sequence.__getslice__(i, j)

    def __getitem__(self, i):
        return self.sequence.__getitem__(i)


class DNA(Molecule):
    def __init__(self, sequence, name = 'dna'):
        Molecule.__init__(self, name)
        self.sequence = sequence

    def get_complement(self):
        """
        Returns:
        ------
        the complement sequence as a string.
        """
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        letters = list(self.sequence)
        letters = [basecomplement[base] if basecomplement.has_key(base) else base for base in letters]
        return ''.join(letters)


class RNA(Molecule):
    def __init__(self, sequence, name = 'rna'):
        Molecule.__init__(self, name)

        for residue in list(sequence):
            self.add_residue(residue)

    def add_residue(self, residue):
        if modified_ribonucleotides.has_key(residue):
            self.modified_residues.append((residue, len(self.sequence)+1))
            residue = modified_ribonucleotides[residue]
        if residue in ['A', 'U', 'G', 'C']:
            self.sequence = ''.join([self.sequence, residue])
        elif residue in ['.', '_', '-']:
            self.sequence = ''.join([self.sequence, '-'])
        else:
            #print "Unknown residue "+residue
            self.sequence = ''.join([self.sequence, residue])

    def get_complement(self):
        """
        Returns:
        ------
        the complement sequence as a string.
        """
        basecomplement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
        letters = list(self.sequence)
        letters = [basecomplement[base] if basecomplement.has_key(base) else base for base in letters]
        return ''.join(letters)    

class SecondaryStructure:

    def __init__(self, rna):
        self.name = "2D"
        self.rna = rna
        self.helices = []
        self.single_strands = []
        self.tertiary_interactions = []
        self.junctions = []
        self.stem_loops = []
        self.source = "N.A:N.A:N.A"
        self._id = str(ObjectId())
        self.__step = None

    def _repr_html_(self):
        if self.__step:
            return self.draw_as_d3()
        else:
            return "No plot available"

    def __walk(self, helix, x_coords, current_y, verbose = False):
        from numpy import mean
        enclosed_stem_loops = []
        if verbose:
            print "walking helix", helix['location']
        #do we have a >= 3-way junction linked to this helix?
        next_junction = None
        for junction in self.junctions:
            junction_location = sorted(junction['location'])
            if len(junction_location) >= 3 and helix['location'][0][-1] == junction_location[0][0] :
                next_junction = junction
                if verbose:
                    print "linked to >=3 junction",junction_location
                #the occupancy will be the number of residues on the largest side
                junction_occupancy = mean([junction_location[-1][-1]-junction_location[-1][0]+1, junction_location[0][-1]-junction_location[0][0]+1])-1
                for i in range(len(junction_location)-1):
                    for h in self.helices: #next helices in junction
                        if h['location'][0][0] == junction_location[i][-1]:
                            if verbose:
                                print "next helix in junction is helix", h['location']
                            self.__walk(h, x_coords, current_y-(helix['location'][0][-1]-helix['location'][0][0])*self.__residue_occupancy-1.5*self.__junction_diameter, verbose)
                    for stem_loop in self.stem_loops: #this helix will lead to which stem loops?
                        if stem_loop['location'][0][0] >= junction_location[i][-1] and stem_loop['location'][-1][-1] <= junction_location[i+1][0]:
                            enclosed_stem_loops.append(stem_loop)
            elif len(junction_location) == 2 and helix['location'][0][-1] == junction_location[0][0]:
                next_junction = junction
                if verbose:
                    print "linked to 2-way junction", junction_location
                #the occupancy will be the number of residues on the largest side
                junction_occupancy = mean([junction_location[-1][-1]-junction_location[-1][0]+1, junction_location[0][-1]-junction_location[0][0]+1])-1
                for h in self.helices:
                    if h['location'][0][0] == junction_location[0][-1]:
                        if verbose:
                            print "next helix in junction is helix", h['location']
                        self.__walk(h, x_coords, current_y-(helix['location'][0][-1]-helix['location'][0][0])*self.__residue_occupancy-1.5*self.__junction_diameter, verbose)
                for stem_loop in self.stem_loops: #this helix will lead to which stem loops?
                    if stem_loop['location'][0][0] >= junction_location[0][0] and stem_loop['location'][-1][-1] <= junction_location[-1][-1]:
                        enclosed_stem_loops.append(stem_loop)
            elif len(junction_location) == 1 and helix['location'][0][-1] == junction_location[0][0]:
                next_junction = junction
                if verbose:
                    print "linked to apical loop", junction_location
        if not len(enclosed_stem_loops): #there was no junction linked to this helix, so it should be in a stem-loop
            for stem_loop in self.stem_loops:
                if helix['location'][0][0] >= stem_loop['location'][0][0] and helix['location'][-1][-1] <= stem_loop['location'][-1][-1]:
                    enclosed_stem_loops.append(stem_loop)
        _x_coords = []
        for enclosed_stem_loop in enclosed_stem_loops:
            _x_coords.append(x_coords[self.stem_loops.index(enclosed_stem_loop)])
        m = mean(_x_coords)
        helix['coords'] = [[m, current_y], [m, current_y-(helix['location'][0][-1]-helix['location'][0][0])*self.__residue_occupancy]]
        if verbose:
            print "helix", helix['location']
            print "coords", helix['coords']
        next_junction['coords'] = [[m, helix['coords'][-1][-1]-self.__junction_diameter/2-5]]
        if verbose:
            print "junction", next_junction['location']
            print "coords", next_junction['coords']

    def compute_plot(self, step = 25, residue_occupancy = 5, junction_diameter = 15, verbose = False):
        if not self.stem_loops:
            self.find_stem_loops()
        self.helices = sorted(self.helices, key=lambda x: x['location'][0][0])
        if verbose:
            print len(self.helices), "helices"
            print len(self.stem_loops), "stem-loops"
            print len(self.single_strands), "single-strands"
            print len(self.junctions), "junctions"
        self.__step = step
        self.__residue_occupancy = residue_occupancy
        self.__junction_diameter = junction_diameter
        if not len(self.helices):
            raise Exception("Your secondary structure contains no helices!!") 
        x = 0
        if verbose:
            print "\nStem-loops placement\n"
        x_coords = []
        if verbose:
            print "stem loop", self.stem_loops[0]['location']
            print "x:", x
        x_coords.append(x)
        for i in range(0, len(self.stem_loops)-1):
            before = self.stem_loops[i]['location'][-1][-1]
            after = self.stem_loops[i+1]['location'][0][0]
            total_residues = 0
            total_junctions = 0
            for junction in self.junctions:
                if len(junction['location']) >=3:
                    for single_strand_location in sorted(junction['location'])[1:-1]: #we only use the single-strands that are not on the left and right "sides"
                        if before <= single_strand_location[0] and after >= single_strand_location[1]:
                            total_residues += single_strand_location[1]-single_strand_location[0]+1
                            total_junctions += 1
            if verbose:
                print "total residues", total_residues
                print "total junctions", total_junctions
            if total_junctions == 0: #we have two stem-loops linked with no junctions. ___||___||___
                x += (after-before+1)*self.__residue_occupancy+self.__junction_diameter
            else:    
                x += total_junctions*self.__step
            if verbose:
                print "stem loop", self.stem_loops[i+1]['location']
                print "x:", x
            x_coords.append(x)

        if verbose:
            print "\nHelices placement\n"
        helix = self.helices[0]
        currentPos = helix['location'][-1][-1]
        current_y = 200
        self.__walk(helix, x_coords, current_y, verbose)
        while currentPos <= len(self.rna):
            currentPos +=1
            if verbose:
                print "currentPos", currentPos
            for helix in self.helices:
                if currentPos == helix['location'][0][0]:
                    current_y = 200
                    self.__walk(helix, x_coords, current_y, verbose)
                    currentPos = helix['location'][-1][-1]
                    break

        single_strands_not_in_junctions = self.single_strands[:]
        for junction in self.junctions:
            for single_strand in junction['single_strands']:
                single_strands_not_in_junctions.remove(single_strand)
                            
        single_strands_not_in_junctions = sorted(single_strands_not_in_junctions)
        
        for single_strand in single_strands_not_in_junctions:
            if verbose:
                print "single strand not in a junction", single_strand['location']
            if single_strand['location'][0] == 1:
                for helix in self.helices:
                    if helix['location'][0][0] == single_strand['location'][-1]+1: 
                        single_strand['coords'] = [[helix['coords'][0][0]-helix['location'][0][0]*self.__residue_occupancy, helix['coords'][0][1]], [helix['coords'][0][0], helix['coords'][0][1]]]
                        break
            elif single_strand['location'][-1] == len(self.rna):
                for helix in self.helices:
                    if helix['location'][-1][-1] == single_strand['location'][0]-1: 
                        single_strand['coords'] = [[helix['coords'][0][0], helix['coords'][0][1]], [helix['coords'][0][0]+(len(self.rna)-helix['location'][-1][-1]+1)*self.__residue_occupancy, helix['coords'][0][1]]]
                        break
            else:
                first_helix = None
                second_helix = None
                for helix in self.helices:
                    if helix['location'][-1][-1] == single_strand['location'][0]-1: 
                        first_helix =  helix
                    if helix['location'][0][0] == single_strand['location'][-1]+1:
                        second_helix = helix
                    if first_helix and second_helix:
                        single_strand['coords'] = [[first_helix['coords'][0][0], first_helix['coords'][0][1]], [second_helix['coords'][0][0], second_helix['coords'][0][1]]]
                        break

    def draw_as_d3(self, stroke_width = 2, verbose = False):
        from pyrna import utils
        from numpy import mean
        all_x =[]
        all_y = []
        quantitative_values = []

        single_strands_not_in_junctions = []

        for single_strand in self.single_strands: #only the single-strands with a coords key are not in junctions
            if single_strand.has_key('coords'):
                single_strands_not_in_junctions.append(single_strand)    
        
        for single_strand in single_strands_not_in_junctions:
            if single_strand.has_key('quantitative_value'):
                quantitative_values.append(single_strand['quantitative_value'])   
            all_x.append(single_strand['coords'][0][0])
            all_y.append(single_strand['coords'][0][1])
            all_x.append(single_strand['coords'][1][0])
            all_y.append(single_strand['coords'][1][1])
        for helix in self.helices:
            if helix.has_key('quantitative_value'):
                quantitative_values.append(helix['quantitative_value'])
            all_x.append(helix['coords'][0][0])
            all_y.append(helix['coords'][0][1])
            all_x.append(helix['coords'][1][0])
            all_y.append(helix['coords'][1][1])
        for junction in self.junctions:
            if junction.has_key('quantitative_value'):
                quantitative_values.append(junction['quantitative_value'])
            all_x.append(junction['coords'][0][0])
            all_y.append(junction['coords'][0][1])

        min_x = min(all_x)
        min_y = min(all_y)
        
        for single_strand in single_strands_not_in_junctions:
            single_strand['coords'] = [[coord[0]-min_x+self.__junction_diameter, coord[1]-min_y+self.__junction_diameter] for coord in single_strand['coords']]
        for helix in self.helices:
            helix['coords'] = [[coord[0]-min_x+self.__junction_diameter, coord[1]-min_y+self.__junction_diameter] for coord in helix['coords']]
        for junction in self.junctions:
            junction['coords'] = [[coord[0]-min_x+self.__junction_diameter, coord[1]-min_y+self.__junction_diameter] for coord in junction['coords']]
        
        all_x = [x-min_x+self.__junction_diameter for x in all_x]
        all_y = [y-min_y+self.__junction_diameter for y in all_y]

        colors_d3 = """"""
        if quantitative_values:
            colors_d3 = """var colors = d3.scale.linear().domain(["""+str(min(quantitative_values))+""","""+str(mean(quantitative_values))+""","""+str(max(quantitative_values))+"""]).range(["#4daf4a",  "#377eb8", "#e41a1c"]);"""
        
        helices_d3 = """"""
        helix_color = '"steelblue"'
        if helix.has_key('quantitative_value'):
            helix_color = "colors("+str(helix['quantitative_value'])+")"
        for helix in self.helices:
            helices_d3 += """svg.append("line")
                            .style("stroke", """+helix_color+""") 
                            .style("stroke-width", """+str(stroke_width)+""")
                            .attr("x1", """+str(helix['coords'][0][0])+""")
                            .attr("y1", """+str(helix['coords'][0][1])+""")
                            .attr("x2", """+str(helix['coords'][1][0])+""")
                            .attr("y2", """+str(helix['coords'][1][1])+""");
                    """

        junctions_d3 = """"""
        for junction in self.junctions:
            if len(junction['location']) >= 3:
                junction_location = sorted(junction['location'])
                for i in range(len(junction_location)-1):
                    for h in self.helices: #next helices in junction
                        if h['location'][0][0] == junction_location[i][-1]:
                            if h['coords'][0][1] != junction['coords'][0][1]: #to avoid to redraw a vertical line
                                new_points = utils.get_points(h['coords'][0][0], h['coords'][0][1], junction['coords'][0][0], junction['coords'][0][1], distance = (self.__junction_diameter+10)/2)
                                if len(new_points) == 2:
                                    helix_color = '"steelblue"'
                                    if h.has_key('quantitative_value'):
                                        helix_color = "colors("+str(h['quantitative_value'])+")"
                                    junctions_d3 += """svg.append("line")
                                        .style("stroke-linecap", "round")
                                        .style("stroke", """+helix_color+""") 
                                        .style("stroke-width", """+str(stroke_width)+""")
                                        .attr("x1", """+str(h['coords'][0][0])+""")
                                        .attr("y1", """+str(h['coords'][0][1])+""")
                                        .attr("x2", """+str(new_points[1][0])+""")
                                        .attr("y2", """+str(new_points[1][1])+""");
                                    """

            junction_color = '"steelblue"'
            if junction.has_key('quantitative_value'):
                junction_color = "colors("+str(junction['quantitative_value'])+")"

            junctions_d3 += """svg.append("circle")
                            .style("fill", """+junction_color+""") 
                            .attr("cx", """+str(junction['coords'][0][0])+""")
                            .attr("cy", """+str(junction['coords'][0][1])+""")
                            .attr("r", """+str(self.__junction_diameter/2)+""");
                    """
            
            junction_color = '"steelblue"'
            if junction.has_key('quantitative_value'):
                junction_color = "colors("+str(junction['quantitative_value'])+")"

            junctions_d3 += """svg.append("circle")
                    .style("fill", "none")
                    .style("stroke", """+junction_color+""")
                    .style("stroke-width", """+str(stroke_width)+""")
                    .attr("cx", """+str(junction['coords'][0][0])+""")
                    .attr("cy", """+str(junction['coords'][0][1])+""")
                    .attr("r", """+str((1.5*self.__junction_diameter)/2)+""");
            """            

        single_strands_d3 = """"""
        for single_strand in single_strands_not_in_junctions:

            single_strand_color = '"steelblue"'
            if single_strand.has_key('quantitative_value'):
                single_strand_color = "colors("+str(single_strand['quantitative_value'])+")"

            single_strands_d3 += """svg.append("line")
                                    .style("stroke-linecap", "round")
                                    .style("stroke", """+single_strand_color+""") 
                                    .style("stroke-width", """+str(stroke_width)+""")
                                    .attr("x1", """+str(single_strand['coords'][0][0])+""")
                                    .attr("y1", """+str(single_strand['coords'][0][1])+""")
                                    .attr("x2", """+str(single_strand['coords'][1][0])+""")
                                    .attr("y2", """+str(single_strand['coords'][1][1])+""");
                                """

        #we end with the helices directly linked at the basis of the drawing
        directly_linked_helices_d3 = """"""
        previous_helix = self.helices[0]
        currentPos = previous_helix['location'][-1][-1]
        while currentPos <= len(self.rna):
            currentPos +=1
            if verbose:
                print "currentPos", currentPos
            for helix in self.helices:
                if currentPos == helix['location'][0][0]:
                    if previous_helix['location'][-1][-1] +1 == helix['location'][0][0]:
                        if verbose:
                            print "directly linked helices", previous_helix['location'] , helix['location']
                        directly_linked_helices_d3 += """svg.append("line")
                                        .style("stroke-linecap", "round")
                                        .style("stroke", "grey") 
                                        .style("stroke-width", """+str(stroke_width)+""")
                                        .attr("x1", """+str(previous_helix['coords'][0][0])+""")
                                        .attr("y1", """+str(previous_helix['coords'][0][1])+""")
                                        .attr("x2", """+str(helix['coords'][0][0])+""")
                                        .attr("y2", """+str(helix['coords'][0][1])+""");
                                    """
                    currentPos = helix['location'][-1][-1]
                    previous_helix = helix
                    break

        d3_description = """
        
            <div id="viz"></div>
            <script type="text/javascript">

            var svg = d3.select("#viz")
                .append("svg")
                .attr("width", """+str(max(all_x)+self.__junction_diameter)+""")
                .attr("height", """+str(max(all_y)+self.__junction_diameter)+""");
                """+colors_d3+"""
                """+junctions_d3+"""
                """+helices_d3+"""
                """+single_strands_d3+"""
                """+directly_linked_helices_d3+"""
                </script>"""

        return d3_description   

    def get_junctions(self):
        return DataFrame(self.junctions)

    def get_paired_residue(self, pos):
        for helix in self.helices:
            if pos >= helix['location'][0][0] and pos <= helix['location'][0][0] + helix['length']-1:
                return helix['location'][-1][-1] - (pos-helix['location'][0][0])
            elif pos <= helix['location'][-1][-1] and pos >= helix['location'][-1][-1] - helix['length']+1:
                return helix['location'][0][0]+ helix['location'][-1][-1] - pos
        return -1

    def find_single_strands(self):
        full_location = Location(start = 1, end = len(self.rna))
        for helix in self.helices:
            full_location =  full_location.remove_location(Location(start = helix['location'][0][0], end = helix['location'][0][-1]))
            full_location =  full_location.remove_location(Location(start = helix['location'][-1][0], end = helix['location'][-1][-1]))
        single_positions = full_location.get_single_positions()
        single_positions.sort()

        start = None
        length = 0
        single_strand_count = 1
        for index, current_pos in enumerate(single_positions):
            if index == 0 or current_pos == single_positions[index-1]+1:
                length += 1
                if index == 0:
                    start = current_pos    
            else:
                self.add_single_strand("SS_%i"%single_strand_count, start, length)
                single_strand_count +=1
                length = 1
                start = current_pos
        #the last
        self.add_single_strand("SS_%i"%single_strand_count, start, length)

    def find_junctions(self):
        self.junctions = []
        for single_strand in self.single_strands:
            if single_strand['location'][0] == 1 or single_strand['location'][-1] == len(self.rna) or len(filter(lambda junction: single_strand in junction['single_strands'], self.junctions)):
                continue
            strands = [single_strand]
            descr = self.rna[single_strand['location'][0]-1:single_strand['location'][-1]]+" "
            current_pos =  self.get_paired_residue(single_strand['location'][-1]+1)+1
            location = [[single_strand['location'][0]-1, single_strand['location'][-1]+1]] 
            next_single_strand = None           

            while current_pos >= 1 and current_pos <= len(self.rna):
                next_single_strand = filter(lambda single_strand : single_strand['location'][0] == current_pos, self.single_strands)
                if next_single_strand and next_single_strand[0] == single_strand:
                    break
                elif next_single_strand:
                    strands.append(next_single_strand[0])
                    location.append([next_single_strand[0]['location'][0]-1, next_single_strand[0]['location'][-1]+1])
                    descr += self.rna[next_single_strand[0]['location'][0]-1:next_single_strand[0]['location'][-1]]+" "
                    current_pos = self.get_paired_residue(next_single_strand[0]['location'][-1]+1)+1
                    continue
                next_helix = filter(lambda helix: current_pos == helix['location'][0][0] or current_pos == helix['location'][-1][-1]-helix['length']+1, self.helices)
                if next_helix:
                    descr += '- '
                    location.append([current_pos-1, current_pos])
                    current_pos = self.get_paired_residue(current_pos)+1

            if next_single_strand and next_single_strand[0] == single_strand:
                self.junctions.append({
                    'single_strands': strands,
                    'description': descr.strip(),
                    'location': location                    
                })

        #now we search for junctions with only directly linked helices
        for helix in self.helices:
            if helix['location'][0][0] == 1 or helix['location'][-1][-1] == len(self.rna) or len(filter(lambda junction: helix['location'][0][0] in sum(junction['location'],[]) or helix['location'][-1][-1] in sum(junction['location'],[]), self.junctions)):
                continue
            descr = ""
            location = []
            next_helix = None
            current_pos = helix['location'][-1][-1]+1

            while current_pos >= 1 and current_pos <= len(self.rna):
                next_helix = filter(lambda helix: current_pos == helix['location'][0][0] or current_pos == helix['location'][1][0], self.helices)
                if next_helix and next_helix[0] == helix:
                    descr += '- '
                    location.append([current_pos-1, current_pos])
                    break
                elif next_helix:
                    descr += '- '
                    location.append([current_pos-1, current_pos])
                    current_pos = self.get_paired_residue(current_pos)+1
                else:
                    break

            if next_helix and next_helix[0] == helix:
                self.junctions.append({
                    'single_strands': [],
                    'description': descr.strip(),
                    'location': location                    
                })  

            #the other side
            descr = ""
            location = []
            next_helix = None
            current_pos = helix['location'][0][1]+1

            while current_pos >= 1 and current_pos <= len(self.rna):
                next_helix = filter(lambda helix: current_pos == helix['location'][0][0] or current_pos == helix['location'][1][0], self.helices)
                if next_helix and next_helix[0] == helix:
                    descr += '- '
                    location.append([current_pos-1, current_pos])
                    break
                elif next_helix:
                    descr += '- '
                    location.append([current_pos-1, current_pos])
                    current_pos = self.get_paired_residue(current_pos)+1
                else:
                    break

            if next_helix and next_helix[0] == helix:
                self.junctions.append({
                    'single_strands': [],
                    'description': descr.strip(),
                    'location': location                    
                })  

        self.junctions = sorted(self.junctions, key=lambda x: x['location'][0][0])      

    def find_stem_loops(self):
        if not self.junctions:
            self.find_junctions()
        #we search for all the stem-loops. A stem loop is a set of contigous helices linked with inner loops and with an apical loop at one end.
        self.stem_loops = []
        ranges = []
        for helix in self.helices:
            #print "helix",helix['location'] 
            start = helix['location'][0][0]
            end = helix['location'][-1][-1]
            #if the helix ends are linked to a junction of degree >= 3 or not linked to any junction, this is a range to keep.
            linked_to_a_junction = False
            for junction in self.junctions:
                for i in range(0, len(junction['location'])-1):
                    if start == junction['location'][i][-1] and end == junction['location'][i+1][0]:
                        if len(junction['location']) >= 3:
                            ranges.append([start, end])
                        linked_to_a_junction = True
                if start == junction['location'][-1][-1] and end == junction['location'][0][0]: #we test the last two ends of the location (first and last values of the matrix)
                    if len(junction['location']) >= 3:
                        ranges.append([start, end])
                    linked_to_a_junction = True
            if not linked_to_a_junction:
                ranges.append([start, end])
        #print ranges

        for _range in ranges:
            start = _range[0]
            end = _range[1]
            #print '\n\n', "search between: ", start,"-", end, '\n\n' 
            enclosed_apical_loops = []
            enclosed_junctions = []
            enclosed_inner_loops = []
            enclosed_helices = []
            for _junction in self.junctions:
                _start = min(_junction['location'])[0] #the lowest end
                _end = max(_junction['location'])[-1] #the highest end
                if _start > start and _end < end:
                    if len(_junction['location']) == 1:
                        enclosed_apical_loops.append(_junction)
                        #print "found apical loop at ", _start, _end
                    elif len(_junction['location']) == 2:
                        enclosed_inner_loops.append(_junction)
                        #print "found inner loop at ", _start, _end
                        #print _junction['location']
                    elif len(_junction['location']) >= 3:
                        enclosed_junctions.append(_junction)
                        #print "found enclosed junction at ", _start, _end
                        #print _junction['location']
            for helix in self.helices:
                _start = helix['location'][0][0]
                _end = helix['location'][-1][-1]
                if _start >= start and _end <= end:
                    enclosed_helices.append(helix)
            #print "enclosed apical loops", len(enclosed_apical_loops)
            #print "enclosed junctions", len(enclosed_junctions)   
            if len(enclosed_apical_loops) == 1 and not enclosed_junctions:
                stem_loop = {'location': [[start, end]]}
                stem_loop['apical_loop'] = enclosed_apical_loops[0]
                stem_loop['inner_loops'] = enclosed_inner_loops
                stem_loop['helices'] = enclosed_helices
                self.stem_loops.append(stem_loop)

        self.stem_loops = sorted(self.stem_loops, key=lambda x: x['apical_loop']['location'][0])

    def find_connected_modules(self):
        self.connected_modules = []
        if not self.junctions:
            self.find_junctions()
        if not self.stem_loops:
            self.find_stem_loops()
               
        tertiary_interactions = self.tertiary_interactions

        for tertiary_interaction in self.tertiary_interactions:
            start = tertiary_interaction['location'][0][0]
            end = tertiary_interaction['location'][-1][-1]
            #print "Tertiary Interaction",start, end
            for stem_loop_1 in self.stem_loops:
                location_1 = Location(nested_lists = stem_loop_1['location'])
                if location_1.has_position(start):
                    for junction in self.junctions:
                        if len(junction['location']) >=3 :
                            location_2 = Location(nested_lists = junction['location'])
                            if location_2.has_position(end):
                                if location_2.end() < location_1.start() or location_2.start() > location_1.end():
                                    self.connected_modules.append((stem_loop_1, junction))
                                    #print location_1.start(),location_1.end() 
                                    #print location_2.start(),location_2.end()
                    for stem_loop_2 in self.stem_loops:
                        location_2 = Location(nested_lists = stem_loop_2['location'])
                        if location_2.has_position(end) and stem_loop_2 != stem_loop_1:
                            if location_2.end() < location_1.start() or location_2.start() > location_1.end():
                                self.connected_modules.append((stem_loop_1, stem_loop_2))
                                #print location_1.start(),location_1.end() 
                                #print location_2.start(),location_2.end()
                if location_1.has_position(end):
                    for junction in self.junctions:
                        if len(junction['location']) >=3 :
                            location_2 = Location(nested_lists = junction['location'])
                            if location_2.has_position(start):
                                if location_2.end() < location_1.start() or location_2.start() > location_1.end():
                                    self.connected_modules.append((stem_loop_1, junction))
                                    #print location_1.start(),location_1.end() 
                                    #print location_2.start(),location_2.end() 

    def add_helix(self, name, start, end, length):
        _ends = [start, start+length-1, end-length+1, end]
        #no pseudoknot allowed
        for helix in self.helices:
            ends = [helix['location'][0][0], helix['location'][0][1], helix['location'][-1][0], helix['location'][-1][-1]]
            if _ends[0] >= ends[1] and _ends[0] <= ends[2] and _ends[3] >= ends[3] or _ends[0] <= ends[0] and _ends[3] >= ends[1] and _ends[3] <= ends[2]: #pseudoknot
                for i in range(0, length):
                    self.add_tertiary_interaction('C', '(', ')', start+i, end-i)
                return None
        helix = {
            'name': name,
            'location': [[start,start+length-1],[end-length+1,end]],
            'length': length,
            'interactions': []
            }
        self.helices.append(helix)
        self.helices = sorted(self.helices, key=lambda helix: helix['location'][0][0]) #the helices are sorted according to the start position
        return helix

    def add_single_strand(self, name, start, length):
        single_strand = {
            'name': name,
            'location': [start,start+length-1]
        };
        self.single_strands.append(single_strand)
        return single_strand

    def add_tertiary_interaction(self, orientation, edge1, edge2, pos1, pos2):
        location = [[pos1, pos1], [pos2, pos2]]
        for tertiary_interaction in self.tertiary_interactions:
            if tertiary_interaction['location'] == location:
                self.tertiary_interactions.remove(tertiary_interaction)
                break    
        self.tertiary_interactions.append({
                            'orientation': orientation, 
                            'edge1': edge1, 
                            'edge2': edge2, 
                            'location': [[pos1, pos1], [pos2, pos2]]
                        })    

    def add_base_pair(self, orientation, edge1, edge2, pos1, pos2):
        is_secondary_interaction = False
        location = [[pos1, pos1], [pos2, pos2]]
        for helix in self.helices:
            start = helix['location'][0][0]
            end = helix['location'][-1][-1]
            length =  helix['length']

            if pos1 >= start and pos1 <= start+length-1:
                diff = pos1 -start
                if end - pos2 == diff:
                    #if not canonical (not AU, GC or GU, neither cWWW, we add it to the helix as a non-canonical secondary interaction
                    if not (self.rna.sequence[pos1-1] == 'A' and self.rna.sequence[pos2-1] == 'U' or \
                            self.rna.sequence[pos1-1] == 'U' and self.rna.sequence[pos2-1] == 'A' or \
                            self.rna.sequence[pos1-1] == 'G' and self.rna.sequence[pos2-1] == 'C' or \
                            self.rna.sequence[pos1-1] == 'C' and self.rna.sequence[pos2-1] == 'G' or \
                            self.rna.sequence[pos1-1] == 'G' and self.rna.sequence[pos2-1] == 'U' or \
                            self.rna.sequence[pos1-1] == 'U' and self.rna.sequence[pos2-1] == 'G') or \
                          orientation != 'C' or edge1 != '(' or edge2 != ')': #we have a non-canonical secondary-interaction
                        
                        for secondary_interaction in helix['interactions']:
                            if secondary_interaction['location'] == location:
                                helix['interactions'].remove(secondary_interaction)
                                break

                        helix['interactions'].append({
                            'orientation': orientation, 
                            'edge1': edge1, 
                            'edge2': edge2, 
                            'location': location
                        })
                    is_secondary_interaction = True
                    break

        if not is_secondary_interaction:
            #if we reach this point, its a tertiary interaction
            self.add_tertiary_interaction(orientation, edge1, edge2, pos1, pos2)

class StructuralAlignment:

    def __init__(self, json_data):
        self.json_data = json_data

    def get_source(self):
        return self.json_data['source']

    def get_aligned_sequences(self):
        rnas = []
        for rna in self.json_data['sequences']:
            rnas.append({'name':rna['name'], 'sequence':rna['sequence']})
        return DataFrame(rnas)

    def get_consensus_2d(self):
        for interaction in self.json_data['consensus2D']:
            interaction['pos1'] = int(interaction['location']['ends'][0][0]);
            interaction['pos2'] = int(interaction['location']['ends'][1][0]);
            del(interaction['location'])    
        return DataFrame(self.json_data['consensus2D']) 


class TertiaryStructure:

    def __init__(self, rna):
        self.source = 'N.A.:N.A.:N.A.'
        self.rna = rna
        self.name = "N.A."
        self.residues = {} #the keys are the absolute position of residues
        self.numbering_system = {}
        self._id = str(ObjectId())

    def get_atoms(self):
        """
        Returns:
        ------
        the description of atoms in a panda dataframe. Columns are:
        - atom name
        - residue absolute position
        - residue position label (according to the numbering system)
        - residue name
        - chain name
        - x (float)
        - y (float)
        - z (float)
        """
        _atoms = []
        keys =[]
        for k in self.residues:
            keys.append(k)

        keys.sort() #the absolute position are sorted

        for key in keys:
            atoms = self.residues[key]['atoms']
            for atom in atoms:
                _atoms.append({
                    'name': atom['name'],
                    'absolute position': key,
                    'position label': self.get_residue_label(key),
                    'residue name': self.rna.sequence[key-1],
                    'chain name': self.rna.name,
                    'x': atom['coords'][0],
                    'y': atom['coords'][1],
                    'z': atom['coords'][2]
                })
        
        return DataFrame(_atoms)

    def add_atom(self, atom_name, absolute_position, coords):
        atom_name = re.sub("\*", "'", atom_name)
        if atom_name == 'OP1':
            atom_name = 'O1P'
        elif atom_name == 'OP2':
            atom_name = 'O2P'
        elif atom_name == 'OP3':
            atom_name = 'O3P'
        if self.residues.has_key(absolute_position):
            self.residues[absolute_position]['atoms'].append({
                    'name': atom_name,
                    'coords': coords
                })
        else:
             self.residues[absolute_position] = {
                'atoms': [{
                    'name': atom_name,
                    'coords': coords
                }]
             }

    def get_residue_label(self, absolute_position):
        if self.numbering_system.has_key(str(absolute_position)):
            return self.numbering_system[str(absolute_position)]
        else:
            return str(absolute_position)          

modified_ribonucleotides = {
    "T": "U",
    "PSU": "U",
    "I": "A",
    "N": "U",
    "S": "U",
    "+A": "A",
    "+C": "C",
    "+G": "G",
    "+I": "I",
    "+T": "U",
    "+U": "U",
    "PU": "A",
    "YG": "G",
    "1AP": "G",
    "1MA": "A",
    "1MG": "G",
    "2DA": "A",
    "2DT": "U",
    "2MA": "A",
    "2MG": "G",
    "4SC": "C",
    "4SU": "U",
    "5IU": "U",
    "5MC": "C",
    "5MU": "U",
    "5NC": "C",
    "6MP": "A",
    "7MG": "G",
    "A23": "A",
    "AD2": "A",
    "AET": "A",
    "AMD": "A",
    "AMP": "A",
    "APN": "A",
    "ATP": "A",
    "AZT": "U",
    "CCC": "C",
    "CMP": "A",
    "CPN": "C",
    "DAD": "A",
    "DCT": "C",
    "DDG": "G",
    "DG3": "G",
    "DHU": "U",
    "DOC": "C",
    "EDA": "A",
    "G7M": "G",
    "GDP": "G",
    "GNP": "G",
    "GPN": "G",
    "GTP": "G",
    "GUN": "G",
    "H2U": "U",
    "HPA": "A",
    "IPN": "U",
    "M2G": "G",
    "MGT": "G",
    "MIA": "A",
    "OMC": "C",
    "OMG": "G",
    "OMU": "U",
    "ONE": "U",
    "P2U": "P",
    "PGP": "G",
    "PPU": "A",
    "PRN": "A",
    "PST": "U",
    "QSI": "A",
    "QUO": "G",
    "RIA": "A",
    "SAH": "A",
    "SAM": "A",
    "T23": "U",
    "T6A": "A",
    "TAF": "U",
    "TLC": "U",
    "TPN": "U",
    "TSP": "U",
    "TTP": "U",
    "UCP": "U",
    "VAA": "A",
    "YYG": "G",
    "70U": "U",
    "12A": "A",
    "2MU": "U",
    "127": "U",
    "125": "U",
    "126": "U",
    "MEP": "U",
    "TLN": "U",
    "ADP": "A",
    "TTE": "U",
    "PYO": "U",
    "SUR": "U",
    "PSD": "A",
    "S4U": "U",
    "CP1": "C",
    "TP1": "U",
    "NEA": "A",
    "GCK": "C",
    "CH": "C",
    "EDC": "G",
    "DFC": "C",
    "DFG": "G",
    "DRT": "U",
    "2AR": "A",
    "8OG": "G",
    "IG": "G",
    "IC": "C",
    "IGU": "G",
    "IMC": "C",
    "GAO": "G",
    "UAR": "U",
    "CAR": "C",
    "PPZ": "A",
    "M1G": "G",
    "ABR": "A",
    "ABS": "A",
    "S6G": "G",
    "HEU": "U",
    "P": "G",
    "DNR": "C",
    "MCY": "C",
    "TCP": "U",
    "LGP": "G",
    "GSR": "G",
    "X": "G",
    "R": "A",
    "Y": "A",
    "E": "A",
    "GSS": "G",
    "THX": "U",
    "6CT": "U",
    "TEP": "G",
    "GN7": "G",
    "FAG": "G",
    "PDU": "U",
    "MA6": "A",
    "UMP": "U",
    "SC": "C",
    "GS": "G",
    "TS": "U",
    "AS": "A",
    "ATD": "U",
    "T3P": "U",
    "5AT": "U",
    "MMT": "U",
    "SRA": "A",
    "6HG": "G",
    "6HC": "C",
    "6HT": "U",
    "6HA": "A",
    "55C": "C",
    "U8U": "U",
    "BRO": "U",
    "BRU": "U",
    "5IT": "U",
    "ADI": "A",
    "5CM": "C",
    "IMP": "G",
    "THM": "U",
    "URI": "U",
    "AMO": "A",
    "FHU": "P",
    "TSB": "A",
    "CMR": "C",
    "RMP": "A",
    "SMP": "A",
    "5HT": "U",
    "RT": "U",
    "MAD": "A",
    "OXG": "G",
    "UDP": "U",
    "6MA": "A",
    "5IC": "C",
    "SPT": "U",
    "TGP": "G",
    "BLS": "A",
    "64T": "U",
    "CB2": "C",
    "DCP": "C",
    "ANG": "G",
    "BRG": "G",
    "Z": "A",
    "AVC": "A",
    "5CG": "G",
    "UDP": "U",
    "UMS": "U",
    "BGM": "G",
    "SMT": "U",
    "DU": "U",
    "CH1": "C",
    "GH3": "G",
    "GNG": "G",
    "TFT": "U",
    "U3H": "U",
    "MRG": "G",
    "ATM": "U",
    "GOM": "A",
    "UBB": "U",
    "A66": "A",
    "T66": "U",
    "C66": "C",
    "3ME": "A",
    "A3P": "A",
    "ANP": "A",
    "FA2": "A",
    "9DG": "G",
    "GMU": "U",
    "UTP": "U",
    "5BU": "U",
    "APC": "A",
    "DI": "I",
    "UR3": "U",
    "3DA": "A",
    "DDY": "C",
    "TTD": "U",
    "TFO": "U",
    "TNV": "U",
    "MTU": "U",
    "6OG": "G",
    "E1X": "A",
    "FOX": "A",
    "CTP": "C",
    "D3T": "U",
    "TPC": "C",
    "7DA": "A",
    "7GU": "U",
    "2PR": "A",
    "CBR": "C",
    "I5C": "C",
    "5FC": "C",
    "GMS": "G",
    "2BT": "U",
    "8FG": "G",
    "MNU": "U",
    "AGS": "A",
    "NMT": "U",
    "NMS": "U",
    "UPG": "U",
    "G2P": "G",
    "2NT": "U",
    "EIT": "U",
    "TFE": "U",
    "P2T": "U",
    "2AT": "U",
    "2GT": "U",
    "2OT": "U",
    "BOE": "U",
    "SFG": "G",
    "CSL": "I",
    "PPW": "G",
    "IU": "U",
    "D5M": "A",
    "ZDU": "U",
    "DGT": "U",
    "UD5": "U",
    "S4C": "C",
    "DTP": "A",
    "5AA": "A",
    "2OP": "A",
    "PO2": "A",
    "DC": "C",
    "DA": "A",
    "LOF": "A",
    "ACA": "A",
    "BTN": "A",
    "PAE": "A",
    "SPS": "A",
    "TSE": "A",
    "A2M": "A",
    "NCO": "A",
    "A5M": "C",
    "M5M": "C",
    "S2M": "U",
    "MSP": "A",
    "P1P": "A",
    "N6G": "G",
    "MA7": "A",
    "FE2": "G",
    "AKG": "G",
    "SIN": "G",
    "PR5": "G",
    "GOL": "G",
    "XCY": "G",
    "5HU": "U",
    "CME": "C",
    "EGL": "G",
    "LC": "C",
    "LHU": "U",
    "LG": "G",
    "PUY": "U",
    "PO4": "U",
    "PQ1": "U",
    "ROB": "U",
    "O2C": "C",
    "C30": "C",
    "C31": "C",
    "C32": "C",
    "C33": "C",
    "C34": "C",
    "C35": "C",
    "C36": "C",
    "C37": "C",
    "C38": "C",
    "C39": "C",
    "C40": "C",
    "C41": "C",
    "C42": "C",
    "C43": "C",
    "C44": "C",
    "C45": "C",
    "C46": "C",
    "C47": "C",
    "C48": "C",
    "C49": "C",
    "C50": "C",
    "A30": "A",
    "A31": "A",
    "A32": "A",
    "A33": "A",
    "A34": "A",
    "A35": "A",
    "A36": "A",
    "A37": "A",
    "A38": "A",
    "A39": "A",
    "A40": "A",
    "A41": "A",
    "A42": "A",
    "A43": "A",
    "A44": "A",
    "A45": "A",
    "A46": "A",
    "A47": "A",
    "A48": "A",
    "A49": "A",
    "A50": "A",
    "G30": "G",
    "G31": "G",
    "G32": "G",
    "G33": "G",
    "G34": "G",
    "G35": "G",
    "G36": "G",
    "G37": "G",
    "G38": "G",
    "G39": "G",
    "G40": "G",
    "G41": "G",
    "G42": "G",
    "G43": "G",
    "G44": "G",
    "G45": "G",
    "G46": "G",
    "G47": "G",
    "G48": "G",
    "G49": "G",
    "G50": "G",
    "T30": "U",
    "T31": "U",
    "T32": "U",
    "T33": "U",
    "T34": "U",
    "T35": "U",
    "T36": "U",
    "T37": "U",
    "T38": "U",
    "T39": "U",
    "T40": "U",
    "T41": "U",
    "T42": "U",
    "T43": "U",
    "T44": "U",
    "T45": "U",
    "T46": "U",
    "T47": "U",
    "T48": "U",
    "T49": "U",
    "T50": "U",
    "U30": "U",
    "U31": "U",
    "U32": "U",
    "U33": "U",
    "U34": "U",
    "U35": "U",
    "U36": "U",
    "U37": "U",
    "U38": "U",
    "U39": "U",
    "U40": "U",
    "U41": "U",
    "U42": "U",
    "U43": "U",
    "U44": "U",
    "U45": "U",
    "U46": "U",
    "U47": "U",
    "U48": "U",
    "U49": "U",
    "U50": "U",
    "UFP": "U",
    "UFR": "U",
    "UCL": "U",
    "3DR": "U",
    "CBV": "C",
    "HFA": "A",
    "MMA": "A",
    "DCZ": "C",
    "GNE": "C",
    "A1P": "A",
    "6IA": "A",
    "CTG": "G",
    "5FU": "U",
    "2AD": "A",
    "T2T": "U",
    "XUG": "G",
    "2ST": "U",
    "5PY": "U",
    "4PC": "C",
    "US1": "U",
    "M5C": "C",
    "DG": "G",
    "DA": "A",
    "DT": "U",
    "DC": "C",
    "P5P": "A",
    "FMU": "U"
}