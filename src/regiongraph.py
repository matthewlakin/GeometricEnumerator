########################################################################
#
# regiongraph.py
#
########################################################################
# 
# GeometricEnumerator
# Copyright (C) 2023 Sarika Kumar & Matthew Lakin
# 
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>. 
# 
########################################################################

try:
    import graphviz
    GraphvizAvailable = True
except:
    GraphvizAvailable = False

from structures import CartesianCoords
import math
import probio_lib as lib


#
# Classes and functions to translate strand graph to region graph.
#

# Class to represent end position of any strand as a pair (s,w) where s is the site and w is a tag w ∈ {5', 3'}
class Position:
    
    def __init__(self, site, tag):
        self.s = site
        self.w = tag

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self.s) +" "+ str(self.w) 

    def __eq__(self, other):
        return (self.s, self.w) == (other.s, other.w)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return (self.s, self.w) < (other.s, other.w)

    def __gt__(self, other):
        return (self.s, self.w) > (other.s, other.w)

############################################################################################################
# Class to represent regions implemented for representing condensing the graph and generating constraints.
#

class Region:

    def __init__(self, sites, comp_sites, totalNucleotideLength, label):
        self.sites = sites
        self.comp_sites = comp_sites
        self.totalNucleotideLength = totalNucleotideLength
        self.label = label

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self.sites) + " " + str(self.comp_sites) + " " + str(self.totalNucleotideLength)

    def __eq__(self, other):
        return (self.sites, self.comp_sites, self.totalNucleotideLength) == (other.sites, other.comp_sites, self.totalNucleotideLength)

    def __ne__(self, other):
        return not self.__eq__(other)

    def isBoundRegion(self):
        return False if self.comp_sites is None else True

    def getRegionLabel(self):
        return self.label

    def regionsites(self):
        if self.isBoundRegion():
            return (self.sites + self.comp_sites)
        else:
            return self.sites
 
############################################################################################################

class RegionEdge:

    def __init__(self, v1, v2, doubleStranded, label, totalNucleotideLength):
        self.v1 = v1
        self.v2 = v2
        self.doubleStranded = doubleStranded
        self.label = label
        self.totalNucleotideLength = totalNucleotideLength

    def getLength(self):
        return self.totalNucleotideLength

############################################################################################################
############################################################################################################
class RegionGraph:

    def __init__(self, vertices_list, edge_list, region_list):
        self.vertices_list = vertices_list
        self.edge_list = edge_list
        self.region_list = region_list # list of edges



    # Make a graphical representation of this region graph.
    def makeGraphicalRepresentation(self):
        if not GraphvizAvailable:
            lib.error('Graphviz module not available!')

        d = graphviz.Graph('d') 
        d.attr('node')

        for s in self.vertices_list:
            d.node(str(s), label= str(' '), width='0.001', height='0.001',
                    fontname='Arial',shape = 'circle',  xlabel = str(s))

        for e in self.edge_list:
            thisColor = 'black:black' if e.doubleStranded else 'black'
            d.edge(str(e.v1), str(e.v2), color = thisColor, label = str(e.label))
        return d

    # Display graphical representation of this region graph.
    def displayRepresentation(self):
        try:
            testJupyter = display
            useGraphviz = True
        except NameError:
            useGraphviz = False
        if useGraphviz:
            display(self.makeGraphicalRepresentation())

    # def findMaxDegreeVertices(self):
    #     max_deg_vertices = []
    #     max_degree = 0
    #     for vertex in self.vertices_list:
    #         num = 0
    #         for edge in self.edge_list:
    #             flag_v1 = True
    #             for l1 in edge.v1:
    #                 if(l1 not in vertex):
    #                     flag_v1 = False
    #             flag_v2 = True
    #             for l2 in edge.v2:
    #                 if(l2 not in vertex):
    #                     flag_v2 = False
    #             if(flag_v1):
    #                 num += 1
    #             if(flag_v2):
    #                 num += 1
    #         if (max_degree < num):
    #             max_degree = num
    #             max_deg_vertices = [vertex]
    #         elif(max_degree == num):
    #             max_deg_vertices.append(vertex)
    #     return max_deg_vertices

    #Finding the vertexes with maximum degree
    def findMaxDegreeVertices(self):
        vertex_degree = {}
        # Go through the list of edges once and record the degree of vertexes in a dict.
        for edge in self.edge_list:
            for v in [edge.v1, edge.v2]:
                if(str(v) in vertex_degree.keys()):
                    item = vertex_degree[str(v)]
                    vertex_degree[str(v)] = (v, (item[1] + 1))
                else:
                    vertex_degree[str(v)] = (v, 1)
        # Find the vertexes with maximum degree
        max_deg_vertices = []
        max_value = 0
        for k, val in vertex_degree.values():
            if (val > max_value):
                max_deg_vertices = [k]
                max_value = val
            elif(val == max_value):
                max_deg_vertices.append(k)
        return max_deg_vertices

    # Compute the angle between the double bonded regions.        
    def computeNickedAngles(self, sampled_strucutres): # max_allowed_Angle  
        nicked_angles = {}     
        for e1 in self.edge_list:
            for e2 in self.edge_list:
                if (e1 != e2):
                    if((not (e1.v1 == e1.v2 or e2.v1 == e2.v2)) and (e1.doubleStranded and e2.doubleStranded)):#self.isNickedRegion(e1, e2)): 
                        coord1 =  sampled_strucutres[str(e1.v1)][0]
                        coord2 =  sampled_strucutres[str(e1.v2)][0]
                        coord3 =  sampled_strucutres[str(e2.v1)][0]
                        coord4 =  sampled_strucutres[str(e2.v2)][0]
                        theta = None
 
                        if(coord1 == coord3):
                            theta = computeAngleBetweenRegions(coord1, coord2, coord4)
                        elif(coord1 == coord4):
                            theta = computeAngleBetweenRegions(coord1, coord2, coord3)
                        elif(coord2 == coord3):
                            theta = computeAngleBetweenRegions(coord2, coord1, coord4)
                        elif(coord2 == coord4):
                            theta = computeAngleBetweenRegions(coord2, coord1, coord3)
                        if (theta is not None): 
                            nicked_angles[str(e1.label) + str(e2.label)] = theta
        return nicked_angles

###########################################################

# Following functions were formerly in the RegionMapping class

# Returns the set of position located at 5' end of given region r 
# In single stranded case, there is only one 5' 
# In double stranded case, the function returns 5' of one end and 5' of complementary strand. 
def five_prPosns(r):
    if r.isBoundRegion():
        return [Position(r.sites[0], "5'"), Position(r.comp_sites[-1], "5'")]
    else:
        return [Position(r.sites[0], "5'")]

# Returns the set of position located at 3' end of given region r  
def three_prPosns(r):
    if r.isBoundRegion():
        return [Position(r.sites[-1], "3'"), Position(r.comp_sites[0], "3'")]
    else:
        return [Position(r.sites[-1], "3'")]

# A and B are the two ends of the regions. 
# This function returns the set of positions at one end of the region
def APosns(r):
    if r.isBoundRegion():
        if (r.sites[0] < r.comp_sites[-1]):
            return [Position(r.sites[0], "5'"), Position(r.comp_sites[0], "3'")]
        elif (r.comp_sites[-1] < r.sites[0]):
            return [Position(r.comp_sites[-1], "5'"), Position(r.sites[-1], "3'")]
        else:
            assert False
    else:
        return [Position(r.sites[0], "5'")]

# This function returns the set of positions at other end of the region. 
def BPosns(r):
    if r.isBoundRegion():
        if (r.sites[-1] < r.comp_sites[0]):
            return [Position(r.sites[-1], "3'"), Position(r.comp_sites[-1], "5'")]
        elif (r.comp_sites[0] < r.sites[-1]):
            return [Position(r.comp_sites[0], "3'"), Position(r.sites[0], "5'")]
        else:
            assert False
    else:
        return [Position(r.sites[-1], "3'")]

# This function makes a check of rules while creating the region graph. 
# These rules are defined in the document. 
def check_rules(v1, v2, region_list, debug=False):
    def debugPrint(x):
        if debug:
            print(x)
    for p1 in v1:
        for p2 in v2:                          
            for i in range(len(region_list)):
                if(p1 in APosns(region_list[i]) and
                   p1 in five_prPosns(region_list[i]) and
                   p2 in APosns(region_list[i]) and
                   p2 in three_prPosns(region_list[i])):
                    return True
                if(p1 in BPosns(region_list[i]) and
                   p1 in five_prPosns(region_list[i]) and
                   p2 in BPosns(region_list[i]) and
                   p2 in three_prPosns(region_list[i])):
                    return True
                if((p1.s.v == p2.s.v) and ((p1.s.n - 1  == p2.s.n) and p1.w == "5'" and p2.w == "3'")):
                    return True
                if((p1.s.v == p2.s.v) and (p1.s.n + 1 == p2.s.n) and (p1.w == "3'" and p2.w == "5'")):
                    return True
    return False

####################################################################################################

# Check to find which vertex from vertices_list contains ALL the positions from posns, if any.
# In a well-constructed region graph, all positions from posns should be in the same vertex.
def findVertexContainingPosns(vertices_list, posns):
    for v1 in vertices_list:
        flag = True
        for p1 in posns:
            if (v1 != []) and (p1 not in v1):
                flag = False
        if (flag):
            return v1
    return None # Flags an error (should not happen in a well-constructed region graph). Dealt with at call site.

# Main function to create region graph from a given strand graph.
def regionGraphFromStrandGraph(sg, debug=False):
    assert sg.isConnected()
    """Creating region graph from strand graph"""
    vertices_list = [] 
    edge_list = []
    region_list = condenseStrandGraphToRegionList(sg) ## TO DO --- FOLD THIS FUNCTION IN HERE.

    ######## Find all positions at end of some region in the list ##########
    all_posns = []
    for i in range(len(region_list)):
        all_posns += five_prPosns(region_list[i])
        all_posns += three_prPosns(region_list[i])
    all_posns = sorted(all_posns)
  
    ######## Create Vertices ##########
    for j in range(len(all_posns)):
        vertices_list.append([all_posns[j]])

    ####### Merge Vertices ########
    while (True):
        merged = False
        for v1 in vertices_list:
            for v2 in vertices_list:
                    if (v1 != v2):
                        if(check_rules(v1, v2, region_list)):
                            v1 += v2
                            vertices_list.remove(v2)
                            merged = True

        if(not merged):
            break

    ######### create edges #############
    for i in range(len(region_list)):
        v1 = findVertexContainingPosns(vertices_list, APosns(region_list[i]))
        v2 = findVertexContainingPosns(vertices_list, BPosns(region_list[i]))
        assert (v1 is not None and v2 is not None)
        edge_list.append(RegionEdge(v1, v2, region_list[i].isBoundRegion(), i, region_list[i].totalNucleotideLength))       

    #rg.displayRepresentation()
    return RegionGraph(vertices_list, edge_list, region_list)

# This function condenses the given sites into a list of regions.
# NB: this function is only called from the regionGraphFromStrandGraph function above,
# and could arguably just be folded into that function as it has no meaning beyond that context.
# TO DO: FIX THIS!!!
def condenseStrandGraphToRegionList(sg, debug=False):
    def debugPrint(x):
        if debug:
            print(x)
    totalNucleotideLength = 0
    regions = []
    visited_sites= []
    sites = sg.getSites()
    i = 0
    label = 0
    while(i < len(sites)):
        if(not sg.siteIsBound(sites[i])):
            reg_ssDNA = [sites[i]]
            debugPrint(i)
            debugPrint(sg.getDomain(sites[i]).name)
            #totalNucleotideLength = sg.getDomain(sites[i]).domainNucleotideLength
            totalNucleotideLength = sg.domainLength[str(sg.getDomain(sites[i]).name)][1]
            while(((i + 1) < len(sites)) and (sites[i].v == sites[i + 1].v) and (not sg.siteIsBound(sites[i+1]))):
                #totalNucleotideLength +=  sg.getDomain(sites[i+1]).domainNucleotideLength
                totalNucleotideLength += sg.domainLength[str(sg.getDomain(sites[i+1]).name)][1]
                reg_ssDNA.append(sites[i+1])
                debugPrint(i+1)
                debugPrint(sg.getDomain(sites[i+1]).name)
                i = i + 1

            regions.append(Region(reg_ssDNA, None, totalNucleotideLength, label))
            debugPrint("label: "+str(label))
            label += 1
            visited_sites += reg_ssDNA
        else:
                d1_comp = sg.getBindingPartner(sites[i])
                reg_dsDNA_s1 = [sites[i]]
                debugPrint(i)
                debugPrint(sg.getDomain(sites[i]).name)
                reg_dsDNA_s2 = [d1_comp]
                debugPrint(sg.getDomain(d1_comp).name)
                
                if(sites[i] not in visited_sites and  d1_comp not in visited_sites):
                    #totalNucleotideLength = sg.getDomain(sites[i]).domainNucleotideLength
                    totalNucleotideLength = sg.domainLength[str(sg.getDomain(sites[i]).name)][1]
                    while(((i+1) < len(sites)) and (sites[i].v == sites[i + 1].v)):
                        d2_comp = sg.getBindingPartner(sites[i + 1])
                        if((d2_comp is not None) and (d1_comp.v == d2_comp.v) and (d1_comp.n == d2_comp.n +1) and (sites[i+1] not in visited_sites and  d2_comp not in visited_sites)):
                            assert (sg.getDomain(sites[i+1]).domainNucleotideLength ==  sg.getDomain(d2_comp).domainNucleotideLength)
                            #totalNucleotideLength +=  sg.getDomain(sites[i+1]).domainNucleotideLength
                            totalNucleotideLength += sg.domainLength[str(sg.getDomain(sites[i+1]).name)][1]
                            reg_dsDNA_s1.append(sites[i + 1])
                            debugPrint(i+1)
                            debugPrint(sg.getDomain(sites[i+1]).name)
                            reg_dsDNA_s2.append(d2_comp)
                            debugPrint(sg.getDomain(d2_comp).name)
                            i = i + 1
                            d1_comp = d2_comp
                        else:
                            break
                    regions.append(Region(reg_dsDNA_s1, reg_dsDNA_s2, totalNucleotideLength, label))
                    debugPrint("label: "+str(label))
                    label += 1
                    visited_sites += reg_dsDNA_s1
                    visited_sites += reg_dsDNA_s2
        i = (i + 1)
    return regions

# This function computes the angle between the regions given their (x, y, z) coordinates.
def computeAngleBetweenRegions(coord1, coord2, coord3):
    vect1 = CartesianCoords(coord2.x - coord1.x, coord2.y - coord1.y, coord2.z - coord1.z)
    vect2 =  CartesianCoords(coord1.x - coord3.x, coord1.y - coord3.y, coord1.z - coord3.z)
    vect1_mag = math.sqrt(vect1.x ** 2 + vect1.y ** 2 + vect1.z ** 2)
    vect2_mag = math.sqrt(vect2.x ** 2 + vect2.y ** 2 + vect2.z ** 2)

    assert (vect1_mag != 0 and vect2_mag != 0)
    val = (vect1.x * vect2.x + vect1.y * vect2.y + vect1.z * vect2.z)/ (vect1_mag * vect2_mag)
    if (math.isclose(val, 1)):
        val = 1
    elif (math.isclose(val, -1)):
        val = -1 
    theta = math.degrees(math.acos(val))
    
    return theta
