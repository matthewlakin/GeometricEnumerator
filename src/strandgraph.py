########################################################################
#
# strandgraph.py
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

from process import *
import probio_lib as lib
import os
try:
    import graphviz
    GraphvizAvailable = True
except:
    GraphvizAvailable = False


############################################################################################################


#
# Class for representing the "sites" used in the strand graph implementation
#
class Site(object):

    def __init__(self, v, n, nmax):
        self.v = v
        self.n = n
        self.nmax = nmax
        # assert self.isValid()

    # def isValid(self, debug=True):
    #     def debugPrint(x):
    #         if debug:
    #             print(x)
    #     if not isinstance(self.v, int):
    #         debugPrint('v is not an integer: '+str(self.v))
    #         return False
    #     if not isinstance(self.n, int):
    #         debugPrint('n is not an integer: '+str(self.n))
    #         return False
    #     if not isinstance(self.nmax, int):
    #         debugPrint('nmax is not an integer: '+str(self.nmax))
    #         return False
    #     if self.v < 0:
    #         debugPrint('v is <0: '+str(self.v))
    #         return False
    #     if self.n < 0:
    #         debugPrint('n is <0: '+str(self.n))
    #         return False
    #     if self.nmax < 1:
    #         debugPrint('nmax is <1: '+str(self.nmax))
    #         return False
    #     if self.n > (self.nmax-1):
    #         debugPrint('n is >(nmax-1): '+str(self.n)+' whereas (nmax-1)='+str(self.nmax-1))
    #         return False
    #     return True

    def __str__(self):
        return self.__repr__()

    #def __repr__(self):
    #    return str((self.v,self.n))+'{'+str(self.nmax)+'}'
    def __repr__(self):
        return '('+str(self.v)+','+str(self.n)+')'

    def __eq__(self, other):
        if self.v == other.v:
            return self.n == other.n
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return (self.v, self.n) < (other.v, other.n)

    def __gt__(self, other):
        return (self.v, self.n) > (other.v, other.n)

    def allFivePrimeSites(self):
        return [Site(self.v, n, self.nmax) for n in range(0, self.n)]
    
    def allThreePrimeSites(self):
        return [Site(self.v, n, self.nmax) for n in range(self.n+1, self.nmax)]

    def fivePrimeAdjacentSite(self):
        return None if self.n <= 0 else Site(self.v, self.n-1, self.nmax)

    def threePrimeAdjacentSite(self):
        return None if self.n >= (self.nmax-1) else Site(self.v, self.n+1, self.nmax)

    # Return all sites between two sites, which should be on the same vertex
    def interveningSitesOnSameVertex(self, other):
        assert self.v == other.v
        assert self.nmax == other.nmax
        if self.n == other.n:
            return []
        (fivePrimeSite,threePrimeSite) = (self,other) if self < other else (other,self)
        return [Site(fivePrimeSite.v, n, fivePrimeSite.nmax) for n in range(fivePrimeSite.n+1, threePrimeSite.n)]

    # Relabel ths site according to the supplied mapping, "vmap".
    # Vmap is a list of indexes. The LIST INDEX of each value in the list
    # is the value that it should be replaced with.
    # Values can be "None", which means that they are there for the purpose
    # of padding the list but should not be looked-up.
    # The vmap can just be the "Visited" list returned from the "enumerateEdges" method below.
    def __relabeled__(self, vmap):
        assert lib.distinct([z for z in vmap if z is not None])
        assert self.v in vmap
        assert self.v is not None
        return Site(vmap.index(self.v), self.n, self.nmax)

############################################################################################################

#
# Class for representing the "edges" used in the strand graph implementation
#
class Edge(object):

    
    def __init__(self, s1, s2): ## Use the ordering on sites to store edges canonically
        if s1 < s2:
            self.s1 = s1
            self.s2 = s2
        else:
            self.s1 = s2
            self.s2 = s1
        # assert self.isValid()
        
    # def isValid(self, debug=True):
    #     def debugPrint(x):
    #         if debug:
    #             print(x)
    #     if not self.s1.isValid():
    #         debugPrint('Site s1 is not valid: '+str(s1))
    #         return False
    #     if not self.s2.isValid():
    #         debugPrint('Site s2 is not valid: '+str(s2))
    #         return False
    #     if self.s1 == self.s2:
    #         debugPrint('Sites s1 and s2 are identical: '+str(s1)+' vs '+str(s2))
    #         return False
    #     return True
        
    def __str__(self):
        return self.__repr__()

    #def __repr__(self):
    #    return str(self.s1)+'<-->'+str(self.s2)
    def __repr__(self):
        return str(self.s1)+'->'+str(self.s2)

    def __eq__(self, other):
        return (self.s1, self.s2) == (other.s1, other.s2)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return (self.s1, self.s2) < (other.s1, other.s2)
    
    def __gt__(self, other):
        return (self.s1, self.s2) > (other.s1, other.s2)

    # Return a new version of this edge, relabeled according to the supplied mapping, "vmap".
    # Vmap is a list of indexes. The LIST INDEX of each value in the list
    # is the value that it should be replaced with.
    # Values can be "None", which means that they are there for the purpose
    # of padding the list but should not be looked-up.
    # The vmap can just be the "Visited" list returned from the "enumerateEdges" method below.
    def __relabeled__(self, vmap):
        assert lib.distinct([z for z in vmap if z is not None])
        return Edge(self.s1.__relabeled__(vmap), self.s2.__relabeled__(vmap))

    def getOutgoingSite(self):
        return self.s1

    def getIncomingSite(self):
        return self.s2

    def getSites(self):
        return (self.s1, self.s2)

    def bothWaysRound(self):
        return [(self.s1, self.s2), (self.s2, self.s1)]

    def withinOneStrand(self):
        return self.s1.v == self.s2.v
###########################################################################################################

############################################################################################################

#
# Class for representing strand graphs
#
class StrandGraph(object):

    def __init__(self, colors_info, vertex_colors, admissible_edges, toehold_edges, current_edges, domainLength):
        self.colors_info = colors_info
        self.vertex_colors = vertex_colors
        self.admissible_edges = admissible_edges
        self.toehold_edges = toehold_edges
        self.current_edges = current_edges
        self.domainLength = domainLength
        # assert self.isValid()

    ####################################################################################################
    #
    # # 
    # # TURNING VALIDITY CHECKING OFF BECAUSE DOING SO SPEEDS UP REACTION ENUMERATION __MASSIVELY__!!!
    # # 
    #
    # def isValid(self, debug=True):
    #     def debugPrint(x):
    #         if debug:
    #             print(x)
    #     # colors_info
    #     if not isinstance(self.colors_info, list):
    #         debugPrint('colors_info is not a list: '+str(self.colors_info))
    #         return False
    #     if len(self.colors_info) == 0:
    #         debugPrint('There are no colors!')
    #         return False
    #     for (c,info) in enumerate(self.colors_info):
    #         if not isinstance(info, dict):
    #             debugPrint('colors_info entry '+str(c)+' is not a dict: '+str(info))
    #             return False
    #         if sorted(info.keys()) != sorted(['length', 'strand_type']):
    #             debugPrint('colors_info info keys are not correct: '+str(sorted(info.keys())))
    #             return False
    #         if c not in self.getColorNumbers():
    #             debugPrint('colors_info color value is not a valid color value: '+str(c))
    #             return False
    #         if not isinstance(info['length'], int):
    #             debugPrint('colors_info length value is not an int: '+str(info['length']))
    #             return False
    #         if info['length'] != len(info['strand_type'].domains):
    #             debugPrint('colors_info length value is incorrect: '+str(info['length'])+' vs '+str(len(info['strand_type'].domains)))
    #             return False
    #         if not isinstance(info['strand_type'], Strand):
    #             debugPrint('colors_info strand_type is not a Strand: '+str(info['strand_type']))
    #             return False
    #     # vertex_colors
    #     if not isinstance(self.vertex_colors, list):
    #         debugPrint('self.vertex_colors is not a list: '+str(self.vertex_colors))
    #         return False
    #     if len(self.vertex_colors) == 0:
    #         debugPrint('There are no vertexes!')
    #         return False
    #     for c in self.vertex_colors:
    #         if not isinstance(c, int):
    #             debugPrint('color is not an int: '+str(c))
    #             return False
    #         if c < 0:
    #             debugPrint('color is <0: '+str(c))
    #             return False
    #         if c >= self.numColors():
    #             debugPrint('color is >num_colors: '+str(c))
    #             return False
    #     # admissible_edges
    #     for e in self.admissible_edges:
    #         if not isinstance(e, Edge):
    #             debugPrint('Admissible edge is not of type Edge: '+str(e))
    #             return False
    #         if e.s1.v not in self.getVertexNumbers():
    #             debugPrint('Edge vertex not in self.getVertexNumbers(): '+str(e.s1.v))
    #             return False
    #         if not (0 <= e.s1.n < self.colors_info[self.vertex_colors[e.s1.v]]['length']):
    #             debugPrint('Edge node index outside of range: '+str(e.s1.n)+' --- should be between 0 and '+str(self.colors_info[self.vertex_colors[e.s1.v]]['length']))
    #             return False
    #         if e.s2.v not in self.getVertexNumbers():
    #             debugPrint('Edge vertex not in self.getVertexNumbers(): '+str(e.s2.v))
    #             return False
    #         if not (0 <= e.s2.n < self.colors_info[self.vertex_colors[e.s2.v]]['length']):
    #             debugPrint('Edge node index outside of range: '+str(e.s2.n)+' --- should be between 0 and '+str(self.colors_info[self.vertex_colors[e.s2.v]]['length']))
    #             return False
    #         if not lib.distinct(self.admissible_edges):
    #             debugPrint('self.admissible_edges contains duplicate elements: '+str(self.admissible_edges))
    #             return False
    #     # toehold_edges
    #     for e in self.toehold_edges:
    #         if not isinstance(e, Edge):
    #             debugPrint('Toehold edge is not of type Edge: '+str(e))
    #             return False
    #         if e not in self.admissible_edges:
    #             debugPrint('Toehold edge is not an admissible edge: '+str(e))
    #             return False
    #         if not self.colors_info[self.vertex_colors[e.s1.v]]['strand_type'].domains[e.s1.n].istoehold:
    #             debugPrint('Toehold edge vertex is not marked as such in colors_info: '+str(e.s1.v))
    #             return False
    #         if not self.colors_info[self.vertex_colors[e.s2.v]]['strand_type'].domains[e.s2.n].istoehold:
    #             debugPrint('Toehold edge vertex is not marked as such in colors_info: '+str(e.s2.v))
    #             return False
    #         if not lib.distinct(self.toehold_edges):
    #             debugPrint('self.toehold_edges contains duplicate elements: '+str(self.toehold_edges))
    #             return False
    #     # current_edges
    #     for e in self.current_edges:
    #         if not isinstance(e, Edge):
    #             debugPrint('Current edge is not of type Edge: '+str(e))
    #             return False
    #         if e not in self.admissible_edges:
    #             debugPrint('Current edge is not an admissible edge: '+str(e))
    #             return False
    #         if not lib.distinct(self.current_edges):
    #             debugPrint('self.current_edges contains duplicate elements: '+str(self.current_edges))
    #             return False
    #     return True
    #
    ####################################################################################################

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        output =  ('{'+os.linesep)
        output += ('  colors_info: '+str(self.colors_info)+os.linesep)
        output += ('  vertex_colors: '+str(self.vertex_colors)+os.linesep)
        output += ('  admissible_edges: '+str(self.admissible_edges)+os.linesep)
        output += ('  toehold_edges: '+str(self.toehold_edges)+os.linesep)
        output += ('  current_edges: '+str(self.current_edges)+os.linesep)
        output += ('  domain length: '+str(self.domainLength)+os.linesep)
        output += ('}')
        return output

    def toProcess(self):
        # assert self.isValid()
        strands = [self.colors_info[self.vertex_colors[vdx]]['strand_type'].copyStrand() for vdx in self.getVertexNumbers()]
        nextBondIdx = 1
        for e in self.current_edges:
            thisBondName = 'i'+str(nextBondIdx)
            nextBondIdx += 1
            new_d1 = strands[e.s1.v].domains[e.s1.n].updateBond(thisBondName)
            strands[e.s1.v] = strands[e.s1.v].modifyDomain(e.s1.n, new_d1)
            new_d2 = strands[e.s2.v].domains[e.s2.n].updateBond(thisBondName)
            strands[e.s2.v] = strands[e.s2.v].modifyDomain(e.s2.n, new_d2)
        p = Process(strands)
        assert p.wellFormed()
        return p

    def printAsProcess(self, useNewlines=False):
        return self.toProcess().compactString(useNewlines=useNewlines)

    def __metric__(self):
        # NB: ordering and equality __CURRENTLY__ only defined for connected strand graphs!
        assert self.isConnected()
        return (self.vertex_colors, self.admissible_edges, self.toehold_edges, self.current_edges)
    
    # def __eq__(self, other):
    #     # NB: equality only defined between strand graphs with compatible colors!
    #     # NB: equality __CURRENTLY__ only defined for connected strand graphs!
    #     assert self.compatibleColors(other)
    #     assert self.isConnected()
    #     assert other.isConnected()
    #     return ((self.vertex_colors == other.vertex_colors) and
    #             (self.admissible_edges == other.admissible_edges) and
    #             (self.toehold_edges == other.toehold_edges) and
    #             (self.current_edges == other.current_edges))

    def __eq__(self, other):
        # NB: equality only defined between strand graphs with compatible colors!
        # NB: equality __CURRENTLY__ only defined for connected strand graphs!
        assert self.compatibleColors(other) 
        assert self.isConnected() 
        assert other.isConnected()
        return self.__metric__() == other.__metric__()

    def __ne__(self, other):
        # NB: inequality only defined between strand graphs with compatible colors!
        # NB: inequality __CURRENTLY__ only defined for connected strand graphs!
        return not(self.__eq__(other))
    
    def __lt__(self, other):
        # NB: ordering only defined between strand graphs with compatible colors!
        # NB: ordering __CURRENTLY__ only defined for connected strand graphs!
        assert self.isConnected()
        assert other.isConnected()
        return self.__metric__() < other.__metric__()

    def __gt__(self, other):
        # NB: ordering only defined between strand graphs with compatible colors!
        # NB: ordering __CURRENTLY__ only defined for connected strand graphs!
        assert self.compatibleColors(other)
        assert self.isConnected()
        assert other.isConnected()
        return self.__metric__() > other.__metric__()

    # Relabel this strand graph according to the supplied mapping, "vmap".
    # Vmap is a list of indexes. The LIST INDEX of each value in the list
    # is the value that it should be replaced with.
    # Values can be "None", which means that they are there for the purpose
    # of padding the list but should not be looked-up.
    # The vmap can just be the "Visited" list returned from the "enumerateEdges" method below.
    def __relabel__(self, vmap):
        assert lib.distinct([z for z in vmap if z is not None])
        self.vertex_colors = [self.vertex_colors[vmap[idx]] for idx in range(len(self.vertex_colors))]
        self.admissible_edges = [e.__relabeled__(vmap) for e in self.admissible_edges]
        self.toehold_edges = [e.__relabeled__(vmap) for e in self.toehold_edges]
        self.current_edges = [e.__relabeled__(vmap) for e in self.current_edges]
        # assert self.isValid()
    
    # RETURN A NEW VERSION of this strand graph that is relabeled according to the supplied mapping, "vmap".
    # Vmap is a list of indexes. The LIST INDEX of each value in the list
    # is the value that it should be replaced with.
    # Values can be "None", which means that they are there for the purpose
    # of padding the list but should not be looked-up.
    # The vmap can just be the "Visited" list returned from the "enumerateEdges" method below.
    def __relabeled__(self, vmap):
        assert lib.distinct([z for z in vmap if z is not None])
        new_colors_info = list(self.colors_info)
        new_vertex_colors = [self.vertex_colors[vmap[idx]] for idx in range(len(self.vertex_colors))]
        new_admissible_edges = [e.__relabeled__(vmap) for e in self.admissible_edges]
        new_toehold_edges = [e.__relabeled__(vmap) for e in self.toehold_edges]
        new_current_edges = [e.__relabeled__(vmap) for e in self.current_edges]
        return StrandGraph(new_colors_info, new_vertex_colors, new_admissible_edges, new_toehold_edges, new_current_edges, self.domainLength)

    # Relabel and sort the graph into a canonical form to simplify equality checking.
    # See Oury 2013 for more details.
    def __getCanonicalRelabeling__(self):
        assert self.isConnected()
        # First, need to figure out the starting vertex(es)
        # Figure out how many strands there are associated with each color in the strand graph,
        # and find the smallest color with the minimal (but non-zero) number of strands.
        min_count = float('inf')
        min_colors = []
        for c in self.getColorNumbers():
            this_count = self.vertex_colors.count(c)
            if this_count > 0: # Don't allow colors with zero occurrences - no vertexes to start from!
                if this_count < min_count:
                    min_count = this_count
                    min_colors = [c]
                elif this_count == min_count:
                    min_colors += [c]
        starting_color = min(min_colors)
        # Then, compute the edge and vertex enumerations starting from each vertex with the chosen color.
        starting_vertexes = [i for (i,vc) in enumerate(self.vertex_colors) if vc == starting_color]
        enumerations = [self.enumerateEdges(sv) for sv in starting_vertexes]
        # Finally, compare the resulting enumerations and return a minimal one.
        enum_min = None
        alpha_min = None
        relabeled_enum_min = None
        # How to compare enumerations: the enumeration is just a list of edges, which have a suitable ordering defined already.
        # They just need to be relabeled under the corresponding relabeling before testing them for equality?!
        def relabelEnum(this_enum, this_alpha):
            return [e.__relabeled__(this_alpha) for e in this_enum]
        for (enum,alpha) in enumerations:
            relabeled_enum = relabelEnum(enum, alpha)
            if enum_min is None:
                assert alpha_min is None
                assert relabeled_enum_min is None
                enum_min = enum
                alpha_min = alpha
                relabeled_enum_min = relabeled_enum
            else:
                if relabeled_enum < relabeled_enum_min:
                    enum_min = enum
                    alpha_min = alpha
                    relabeled_enum_min = relabeled_enum
        assert enum_min is not None
        assert alpha_min is not None
        assert relabeled_enum_min is not None
        return alpha_min

    def __convertToCanonicalForm__(self):
        alpha_min = self.__getCanonicalRelabeling__()
        #print('&&&&&&&&&& canonical alpha = '+str(alpha_min))
        self.__relabel__(alpha_min)
        self.admissible_edges.sort()
        self.toehold_edges.sort()
        self.current_edges.sort()
        
    def numVertexes(self):
        return len(self.vertex_colors)

    def getVertexNumbers(self):
        return list(range(self.numVertexes()))

    def numColors(self):
        return len(self.colors_info)

    def getColorNumbers(self):
        return list(range(self.numColors()))

    def numCurrentEdges(self):
        return len(self.current_edges)

    # def siteIsValid(self, s, debug=False):
    #     def debugPrint(x):
    #         if debug:
    #             print(x)
    #     if not isinstance(s, Site):
    #         debugPring('s is not a Site: '+str(s))
    #         return False
    #     if not (0 <= s.v < self.numVertexes()):
    #         print('s.v is outside allowable range of vertexes: '+str(s))
    #         return False
    #     if not (0 <= s.n < self.colors_info[self.vertex_colors[s.v]]['length']):
    #         print('s.n is outside allowable range of sites for specified vertex: '+str(s))
    #         return False
    #     if not (s.nmax == self.colors_info[self.vertex_colors[s.v]]['length']):
    #         print('s.nmax does not match: '+str(s)+' vs '+str(self.colors_info[self.vertex_colors[s.v]]['length']))
    #         return False
    #     return True
    
    def getSites(self):
        res = []
        for vidx in self.getVertexNumbers():
            this_num_sites = self.colors_info[self.vertex_colors[vidx]]['length']
            for nidx in range(this_num_sites):
                this_site = Site(vidx,nidx,this_num_sites)
                # assert self.siteIsValid(this_site)
                res += [this_site]
        return res

    def getDomain(self, s):
        # assert self.siteIsValid(s)
        return self.colors_info[self.vertex_colors[s.v]]['strand_type'].domains[s.n]
    
    def compatibleColors(self, other, debug=False):
        def debugPrint(x):
            if debug:
                print(x)
        assert isinstance(other, StrandGraph)
        # assert self.isValid()
        # assert other.isValid()
        if len(self.colors_info) != len(other.colors_info):
            debugPrint('Lengths of colors_info lists do not match: '+str(self.colors_info)+' vs '+str(other.colors_info))
            return False
        for (c1,c2) in zip(self.colors_info, other.colors_info):
            if c1['length'] != c2['length']:
                debugPrint('Lengths assigned to following colors do not match: '+str(c1)+' vs '+str(c2))
                return False
            if c1['strand_type'] != c2['strand_type']:
                debugPrint('Strand types assigned to following colors do not match: '+str(c1)+' vs '+str(c2))
                return False
        return True

    def getExternalEdgesFromVertexes(self, vdxs):
        res = []
        for e in self.current_edges:
            if ((e.s1.v in vdxs and e.s2.v not in vdxs) or
                (e.s2.v in vdxs and e.s1.v not in vdxs)):
                res += [e]
        return res

    def __makeVertexPartitions__(self):
        def tryToMergeVertexPartitions(vertex_partitions):
            for (idx,p1) in enumerate(vertex_partitions):
                for (jdx,p2) in enumerate(vertex_partitions):
                    if idx <= jdx:
                        pass
                    else:
                        for e in self.current_edges:
                            if ((e.s1.v in vertex_partitions[idx] and e.s2.v in vertex_partitions[jdx]) or
                                (e.s2.v in vertex_partitions[idx] and e.s1.v in vertex_partitions[jdx])):
                                return lib.mergeInnerLists(vertex_partitions, [idx, jdx])
            return None
        vertex_partitions = lib.explode(self.getVertexNumbers())
        while True:
            new_vertex_partitions = tryToMergeVertexPartitions(vertex_partitions)
            if new_vertex_partitions is None:
                break
            else:
                vertex_partitions = new_vertex_partitions
        return vertex_partitions

    def isConnected(self):
        return len(self.__makeVertexPartitions__()) == 1
    
    def connectedComponents(self):
        def filterConvertAndMaybeCheckEdges(edges, vs, doCheck):
            res = []
            for e in edges:
                if ((e.s1.v in vs) and (e.s2.v in vs)):
                    res += [e.__relabeled__(vs)]
                else:
                    if doCheck: # Make sure that edge is completely inside or completely outside the component, if doCheck is True...
                        assert ((e.s1.v not in vs) and (e.s2.v not in vs))
            return res
        def makeStrandGraphFromVertexPartition(vs):
            new_vertex_colors = [self.vertex_colors[v] for v in vs]
            new_admissible_edges = filterConvertAndMaybeCheckEdges(self.admissible_edges, vs, False)
            new_toehold_edges = filterConvertAndMaybeCheckEdges(self.toehold_edges, vs, False)
            new_current_edges = filterConvertAndMaybeCheckEdges(self.current_edges, vs, True)

            new_sg = StrandGraph(self.colors_info, new_vertex_colors, new_admissible_edges, new_toehold_edges, new_current_edges, self.domainLength)
            assert new_sg.isConnected()
            new_sg.__convertToCanonicalForm__()
            return new_sg
        return [makeStrandGraphFromVertexPartition(vs) for vs in self.__makeVertexPartitions__()]

    # def connectedComponents_another(self):
    #     def filterConvertAndMaybeCheckEdges(edges, vs, doCheck):
    #         res = []
    #         for e in edges:
    #             if ((e.s1.v in vs) and (e.s2.v in vs)):
    #                 res += [e.__relabeled__(vs)]
    #             else:
    #                 if doCheck: # Make sure that edge is completely inside or completely outside the component, if doCheck is True...
    #                     assert ((e.s1.v not in vs) and (e.s2.v not in vs))
    #         return res
    #     def makeStrandGraphFromVertexPartition(vs):
    #         new_vertex_colors = [self.vertex_colors[v] for v in vs]
    #         new_admissible_edges = filterConvertAndMaybeCheckEdges(self.admissible_edges, vs, False)
    #         new_toehold_edges = filterConvertAndMaybeCheckEdges(self.toehold_edges, vs, False)
    #         new_current_edges = filterConvertAndMaybeCheckEdges(self.current_edges, vs, True)

    #         new_sg = StrandGraph(self.colors_info, new_vertex_colors, new_admissible_edges, new_toehold_edges, new_current_edges, self.domainLength)
    #         assert new_sg.isConnected()
    #         return new_sg
    #     return [makeStrandGraphFromVertexPartition(vs) for vs in self.__makeVertexPartitions__()]


    def compose(self, other):
        assert self.compatibleColors(other)
        vmap = ([None] * self.numVertexes()) + list(other.getVertexNumbers())
        new_vertex_colors = list(self.vertex_colors) + list(other.vertex_colors)
        extra_admissible_edges = []
        extra_toehold_edges = []
        for s1 in self.getSites():
            d1 = self.getDomain(s1)
            for s2 in other.getSites():
                d2 = other.getDomain(s2)
                if d1.isComplementaryTo(d2):
                    new_edge = Edge(s1, s2.__relabeled__(vmap))
                    extra_admissible_edges += [new_edge]
                    if d1.istoehold and d2.istoehold:
                        extra_toehold_edges += [new_edge]
        new_admissible_edges = list(self.admissible_edges) + [e.__relabeled__(vmap) for e in other.admissible_edges] + extra_admissible_edges
        new_toehold_edges = list(self.toehold_edges) + [e.__relabeled__(vmap) for e in other.toehold_edges] + extra_toehold_edges
        new_current_edges = list(self.current_edges) + [e.__relabeled__(vmap) for e in other.current_edges]
        new_sg = StrandGraph(self.colors_info, new_vertex_colors, new_admissible_edges, new_toehold_edges, new_current_edges, self.domainLength)
        return new_sg

    # def getOutEdges(self, v):
    #     assert v in self.getVertexNumbers()
    #     edges = []
    #     for e in self.current_edges:
    #         if e.s1.v == v:
    #             edges += [e]
    #     return edges

    # def getInEdges(self, v):
    #     assert v in self.getVertexNumbers()
    #     edges = []
    #     for e in self.current_edges:
    #         if e.s2.v == v:
    #             edges += [e]
    #     return edges

    # This returns a list of the current edges that interact with the vertex v in this strand graph.
    # They are ordered based on the site number that they join to (into / out of, doesn't matter)
    # on the specified vertex v.
    def getLocallySortedCurrentEdges(self, v):
        assert v in self.getVertexNumbers()
        edges = []
        for e in self.current_edges:
            if e.s1.v == v or e.s2.v == v:
                edges += [e]
        def sortingKey(e):
            if e.s1.v == v:
                return e.s1.n
            elif e.s2.v == v:
                return e.s2.n
            else:
                assert False
        return sorted(edges, key=sortingKey)

    # # Edges are assigned colors by lifting the colors of vertexes.
    # # By convention, we use the color of the "out" vertex (i.e., the smaller site).
    # def getEdgeColor(self, e):
    #     assert e in self.admissible_edges
    #     return self.vertex_colors[e.s1.v]

    # Enumerate edges (and vertexes) from the given starting vertex,
    # via a deterministic graph traversal algorithm.
    # See Oury 2013 for more details.
    def enumerateEdges(self, startVertex):
        assert self.isConnected()
        assert startVertex in self.getVertexNumbers()
        Enum = []
        Q = [startVertex]
        Visited = [startVertex]
        while Q != []:
            v = Q.pop(0)
            ############################################################################################################
            # outEdgesSortedByOutgoingSite = sorted(self.getOutEdges(v), key=lambda e: e.getOutgoingSite())
            # inEdgesSortedByOutgoingSite = sorted(self.getInEdges(v), key=lambda e: e.getOutgoingSite())
            # edges = outEdgesSortedByOutgoingSite + inEdgesSortedByOutgoingSite
            ############################################################################################################
            # # IT SEEMS THAT A PROBLEM WITH THE ABOVE IS THAT WE SORT THE SITES IN AN EDGE SO THAT THE "OUTGOING SITE"
            # # HAS THE SMALLER VERTEX NUMBER WHEREAS THE "INCOMING SITE" HAS THE LARGER VERTEX NUMBER.
            # # SINCE WE WERE DOING OUT-EDGES BEFORE IN-EDGES, THIS MEANT THAT THE ORDER IN WHICH WE EXPLORE THE EDGES
            # # CAN CHANGE DEPENDING ON THE VERTEX NUMBER OF THE OTHER VERTEXES IN THE SYSTEM!!!
            # # I THINK THAT THE UNDERLYING GRAPH IMPLEMENTATION IN THE OURY PAPER MUST HAVE BEEN SLIGHTLY DIFFERENT.
            # # THEREFORE, IN THE VERSION BELOW I AM TRYING A DIFFERENT ORDERING ON THESE EDGES, WHICH IS DERIVED
            # # ONLY FROM THEIR RELATIONSHIP WITH THE CURRENT VERTEX IN QUESTION...
            ############################################################################################################
            edges = self.getLocallySortedCurrentEdges(v)
            ############################################################################################################
            for e in edges:
                if e not in Enum:
                    Enum += [e]
                    if v == e.s2.v:
                        vnew = e.s1.v
                    elif v == e.s1.v:
                        vnew = e.s2.v
                    else:
                        assert False
                    if vnew not in Visited:
                        Visited += [vnew]
                        Q += [vnew]
        assert len(Enum) == self.numCurrentEdges()
        assert lib.distinct(Enum)
        #print('&&&&& Enum = '+str(Enum))
        #print('&&&&& Visited = '+str(Visited))
        return (Enum, Visited) # Enum is ordering on edges, Visited is ordering on vertexes -> "vertex alpha-renaming" from the Oury paper.

    def siteIsBound(self, s):
        assert s in self.getSites()
        for e in self.current_edges:
            if s in [e.s1, e.s2]:
                return True
        return False

    def currentlyUnboundSites(self):
        unboundSites = self.getSites()
        for s in self.currentlyBoundSites():
            unboundSites.remove(s)
        return unboundSites
    
    def currentlyBoundSites(self):
        currently_bound_sites = []
        for e in self.current_edges:
            currently_bound_sites += [e.s1, e.s2] ## Don't need to check for duplicates here...
        return currently_bound_sites

    def possibleNewEdges(self):
        possible_new_edges = []
        for e in self.admissible_edges:
            if e not in self.current_edges:
                possible_new_edges.append(e)
        return possible_new_edges

    def addEdgeToCurrentEdges(self, e):
        assert e in self.admissible_edges
        new_current_edges = list(self.current_edges) + [e]
        return StrandGraph(self.colors_info, self.vertex_colors, self.admissible_edges, self.toehold_edges, new_current_edges,self.domainLength)

    def removeEdgeFromCurrentEdges(self, e):
        assert e in self.admissible_edges
        assert e in self.current_edges
        new_current_edges = []
        for ce in self.current_edges:
            if ce != e:
                new_current_edges.append(ce)
        assert len(new_current_edges) == (len(self.current_edges) - 1)
        return StrandGraph(self.colors_info, self.vertex_colors, self.admissible_edges, self.toehold_edges, new_current_edges, self.domainLength)

    def possibleAdjacentEdges(self, e):
        s1_5pr = e.s1.fivePrimeAdjacentSite()
        s1_3pr = e.s1.threePrimeAdjacentSite()
        s2_5pr = e.s2.fivePrimeAdjacentSite()
        s2_3pr = e.s2.threePrimeAdjacentSite()
        edges_to_check = []
        if s1_5pr is not None and s2_3pr is not None:
            edges_to_check.append(Edge(s1_5pr, s2_3pr))
        if s1_3pr is not None and s2_5pr is not None:
            edges_to_check.append(Edge(s1_3pr, s2_5pr))
        res = []
        for new_edge in edges_to_check:
            if new_edge in self.admissible_edges:
                res.append(new_edge)
        return res

    def has_adjacent(self, e):
        assert e in self.admissible_edges
        adjacent_possibilities = self.possibleAdjacentEdges(e)
        for poss_edge in adjacent_possibilities:
            if poss_edge in self.current_edges:
                return True
        return False

    def getBindingPartner(self, s):
        for e in self.current_edges:
            if e.s1 == s:
                return e.s2
            elif e.s2 == s:
                return e.s1
        return None

    # Find bound sites on same vertex as a given site.
    # For convenience later on, we also split these out depending on
    # whether they are located toward the 5' or 3' end from the specified site.
    def boundSitesFivePrimeFrom(self, s):
        return [r for r in s.allFivePrimeSites() if self.siteIsBound(r)]
    def boundSitesThreePrimeFrom(self, s):
        return [r for r in s.allThreePrimeSites() if self.siteIsBound(r)]
    def boundSitesOnSameVertexAs(self, s):
        return self.boundSitesFivePrimeFrom(s) + self.boundSitesThreePrimeFrom(s)

    ################################################################################################################################################
    #
    # NB: Giving up using subgraph clusters to draw reactions, so don't need the commented code below, at least for now.
    # 
    ################################################################################################################################################

    # # Inner function to make a graphical representation of this strand graph!
    # # This function takes a graphviz subgraph in which to draw the representation.
    # # This allows it to be used to also draw reaction subgraphs, for example
    # # Kwargs specify cluster style (by default, transparent).
    # def __fillInGraphicalRepresentation__(self, cl, prefix='', style='filled', outlinecolor='transparent', fillcolor='transparent'):
    #     strandColors = ['red', 'blue', 'darkgreen', 'goldenrod', 'darkviolet',
    #                     'orange', 'maroon', 'orchid', 'limegreen', 'salmon',
    #                     'brown', 'aquamarine3', 'slateblue', 'hotpink2', 'darkgoldenrod1',
    #                     'darkorchid4', 'darkolivegreen3', 'firebrick2', 'deeppink2', 'mediumturquoise']
    #     def getColor(s):
    #         thisVertexColor = self.vertex_colors[s.v]
    #         return strandColors[thisVertexColor % len(strandColors)]
    #     cl.attr(color=outlinecolor)
    #     cl.attr(style=style)
    #     cl.attr(fillcolor=fillcolor)
    #     cl.attr(nodesep='0.1')
    #     nodeNames = []
    #     for s in self.getSites():
    #         thisColor = getColor(s)
    #         siteIsToehold = self.colors_info[self.vertex_colors[s.v]]['strand_type'].domains[s.n].istoehold
    #         thisShape = 'ellipse' if siteIsToehold else 'rectangle'
    #         thisNodeName = prefix+str(s)
    #         cl.node(thisNodeName, str(self.getDomain(s)), width='0.5', height='0.3',
    #                 color=thisColor, fontcolor=thisColor, fontname='Arial', shape=thisShape)
    #         nodeNames += [thisNodeName]
    #     for s in self.getSites():
    #         t = s.threePrimeAdjacentSite()
    #         if t is not None:
    #             assert s.v == t.v
    #             c = getColor(s)
    #             cl.edge(prefix+str(s), prefix+str(t), color=c, len='1.0')
    #     for e in self.current_edges:
    #         thisStyle = 'dashed'
    #         thisColor = 'lightgrey' if e in self.toehold_edges else 'grey'
    #         cl.edge(prefix+str(e.s1), prefix+str(e.s2), dir='none', color=thisColor, style=thisStyle, len='0.75')
    #     return nodeNames
    
    # # Make a graphical representation of this strand graph!
    # # Returns the graphviz object. This allows them to be rendered in a Jupyter notebook via display(), for example.
    # # See here for documentation: https://graphviz.readthedocs.io/en/stable/manual.html
    # def makeGraphicalRepresentation(self):
    #     if not GraphvizAvailable:
    #         lib.error('Graphviz module not available!')
    #     d = graphviz.Digraph(engine='neato') # 'fdp' -- maybe better for implementing clusters etc?
    #     d.attr(colorscheme='X11')
    #     #d.attr(nodesep='1.0')
    #     #d.attr(rankdir='LR')
    #     with d.subgraph(name='cluster0') as cl:
    #         self.__fillInGraphicalRepresentation__(cl)
    #     return d

    ############################################################################################################

    # Make a graphical representation of this strand graph!
    # Returns the graphviz object. This allows them to be rendered in a Jupyter notebook via display(), for example.
    # See here for documentation: https://graphviz.readthedocs.io/en/stable/manual.html
    def makeGraphicalRepresentation(self):
        if not GraphvizAvailable:
            lib.error('Graphviz module not available!')
        strandColors = ['red', 'blue', 'darkgreen', 'goldenrod', 'darkviolet', 'orange', 'maroon', 'orchid', 'limegreen', 'salmon',
                        'brown', 'aquamarine3', 'slateblue', 'hotpink2', 'darkgoldenrod1', 'darkorchid4', 'darkolivegreen3', 'firebrick2', 'deeppink2', 'mediumturquoise']
        def getColor(s):
            thisVertexColor = self.vertex_colors[s.v]
            return strandColors[thisVertexColor % len(strandColors)]
        d = graphviz.Digraph(engine='neato') # 'fdp' -- maybe better for implementing clusters etc?
        d.attr(colorscheme='X11')
        #d.attr(rankdir='LR')
        d.attr(nodesep='0.1')
        for s in self.getSites():
            thisColor = getColor(s)
            siteIsToehold = self.colors_info[self.vertex_colors[s.v]]['strand_type'].domains[s.n].istoehold
            thisShape = 'ellipse' if siteIsToehold else 'rectangle'
            d.node(str(s), str(self.getDomain(s)), width='0.5', height='0.3',
                   color=thisColor, fontcolor=thisColor, fontname='Arial', shape=thisShape)
        for s in self.getSites():
            t = s.threePrimeAdjacentSite()
            if t is not None:
                assert s.v == t.v
                c = getColor(s)
                d.edge(str(s), str(t), color=c, len='1.0')
        for e in self.current_edges:
            thisStyle = 'dashed'
            #thisColor = 'black:black' if e in self.toehold_edges else 'black:black'
            thisColor = 'lightgrey' if e in self.toehold_edges else 'grey'
            d.edge(str(e.s1), str(e.s2), dir='none', color=thisColor, style=thisStyle, len='0.75')
        return d
    
    # Display graphical representation of this strand graph.
    # Tries to check if display() is bound, so if we are in the terminal, so if we are in a Jupyter notebook, this should draw the visualization.
    # If we are running from the terminal, this should not do anything.
    def displayRepresentation(self):
        try:
            testJupyter = display
            useGraphviz = True
        except NameError: # Should jump back to here if not in Jupyter, as display won't be defined by default. NB: any other NameError will be picked up too, however...
            useGraphviz = False
        if useGraphviz:
            display(self.makeGraphicalRepresentation())
        else:
            print(self)

    def sameSpecies(self, s1, s2):
        vertex_list = self.__makeVertexPartitions__()
        for vertex in vertex_list:
            if(s1.v in vertex and s2.v in vertex):
                return True
        return False


############################################################################################################

def makeEmptyStrandGraph(colors_info):
    return StrandGraph(colors_info, [], [], [], [], {})

############################################################################################################

def strandGraphComponentsFromProcess(p):
    assert isinstance(p, Process)
    assert p.wellFormed()
    colors_info = []
    strand_types_seen = []
    for s in p.strands:
        this_strand_type = s.strandType()
        if this_strand_type not in strand_types_seen:
            colors_info.append({'strand_type':this_strand_type, 'length':len(s.domains)})
            strand_types_seen.append(this_strand_type)
    colors_info.sort(key=lambda d: d['strand_type']) # Ensure a canonical color representation for a given set of strand types
    def findColor(strand_type, colors_info):
        color_idx = None
        for (idx,d) in enumerate(colors_info):
            if d['strand_type'] == strand_type:
                color_idx = idx
                break
        assert color_idx is not None
        return color_idx
    vertex_colors = [findColor(s.strandType(), colors_info) for s in p.strands]
    admissible_edges = []
    for (vdx1,c1) in enumerate(vertex_colors):
        for (ddx1,d1) in enumerate(p.strands[vdx1].domains):
            for (vdx2,c2) in enumerate(vertex_colors):
                for (ddx2,d2) in enumerate(p.strands[vdx2].domains):
                    if d1.stripBond().isComplementaryTo(d2.stripBond()):
                        this_edge = Edge(Site(vdx1, ddx1, colors_info[c1]['length']),
                                         Site(vdx2, ddx2, colors_info[c2]['length']))
                        if this_edge not in admissible_edges:
                            admissible_edges.append(this_edge)
    toehold_edges = []
    for e in admissible_edges:
        d1 = p.strands[e.s1.v].domains[e.s1.n]
        d2 = p.strands[e.s2.v].domains[e.s2.n]
        if d1.istoehold and d2.istoehold:
            toehold_edges.append(e)
    current_edges = []
    for e in admissible_edges:
        d1 = p.strands[e.s1.v].domains[e.s1.n]
        d2 = p.strands[e.s2.v].domains[e.s2.n]
        if d1.wellFormedBondTo(d2):
            current_edges.append(e)
    domainLength = {}
    return (colors_info, vertex_colors, admissible_edges, toehold_edges, current_edges, domainLength)

def strandGraphFromProcess(p):
    assert isinstance(p, Process)
    assert p.wellFormed()
    (colors_info, vertex_colors, admissible_edges, toehold_edges, current_edges, domainLength) = strandGraphComponentsFromProcess(p)
    return StrandGraph(colors_info, vertex_colors, admissible_edges, toehold_edges, current_edges, domainLength)

def connectedStrandGraphsFromProcess(p):
    assert isinstance(p, Process)
    assert p.wellFormed()
    sg = strandGraphFromProcess(p)
    return [sg] if sg.isConnected() else sg.connectedComponents()

############################################################################################################

def composeMultipleStrandGraphs(sgs):
    assert sgs != []
    new_sg = None
    for this_sg in sgs:
        if new_sg is None:
            new_sg = this_sg
        else:
            new_sg = new_sg.compose(this_sg)
    assert new_sg is not None
    return new_sg

def parseDomainLength(domainLengthString):
    domainLength = {}
    # split the string based on spaces
    split_string = domainLengthString.split()
    for count, item in enumerate(split_string):
        if item == 'toeholdDomain':
            if (split_string[count + 2] == 'length'):
                domainLength[split_string[count + 1]] = (True, int(split_string[count + 3]))
        elif item == 'longDomain':
            if (split_string[count + 2] == 'length'):
                domainLength[split_string[count + 1]] = (False, int(split_string[count + 3]))
    return domainLength




