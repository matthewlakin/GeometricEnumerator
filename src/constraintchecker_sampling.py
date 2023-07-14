########################################################################
#
# constraintchecker_sampling.py
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

import math
import random
import matplotlib.pyplot as plt
from constraintchecker_abstract import *
from constants import *
from distributions import *
from domain import *
from angle_distributions import *
from length_distributions import *
from structures import *
from regiongraph import *


class ConstraintChecker_Sampling(ConstraintChecker_Abstract):

    def __init__(self, seed=None):
        super().__init__()
        self.reseed(seed=seed)
        self.ssDomainLengthDist = WormLikeChainLengthDistribution() #UniformLengthDistribution()
        self.dsDomainLengthDist = MaxLengthDistribution()
        self.tetherAngleDist = UniformSphereAngleDistribution() # UniformHemisphereAngleDistribution() # No tethering
        self.ssDomainAngleDist = UniformSphereAngleDistribution()
        self.dsdsDomainAngleDist = UniformSphereAngleDistribution() #NickedAngleDistribution() #UniformSphereAngleDistribution()

    def reseed(self, seed=None):
        if seed is None:
            self.prng = random.Random()
        else:
            self.prng = random.Random(seed)
            
    # To check whether a strand graph is physically possible or not.
    def isPlausible(self, sg, debug=False):
        def debugPrint(x):
            if debug:
                print(x)

        unsuccessful_trials = 0
        #sg.displayRepresentation()
        # Plausability is only checked for connected strand graph
        if (sg.isConnected()):
            # Convert strand graph to region graph
            rg = regionGraphFromStrandGraph(sg)
            for i in range(SAMPLING_TRIALS):             
                # Find the physical coordinates of the vertices in the region graph
                sampled_structures = self.sampleCoordinates(rg)
                # check if all of the constraints are satisfied simultaneously 
                flag = self.checkConstraints(rg, sampled_structures)
                if (flag):
                    debugPrint("Satisfiable!!!!---Sampling")
                    debugPrint("numer of unsuccessful trials  " + str(i))  
                    sampling_info = {'sampling_unsuccessful_trials': unsuccessful_trials}
                    return (True, sampling_info)
                unsuccessful_trials += 1
        debugPrint("UnSatisfiable!!!!---Sampling")
        debugPrint("number of unsuccessful trials  " + str(unsuccessful_trials))
        sampling_info = {'sampling_unsuccessful_trials': unsuccessful_trials}
        # Structure is not plausible
        return (False, sampling_info)
  
    # Find the coordiantes of vertices in a given region graph.
    def sampleCoordinates(self, rg):
 
        # create four set of regions
        dsDNA_regions = []
        ssDNA_regions = []
        unprocessed_regions = []
        skipped_regions = []

        # find maximum degree vertex
        max_deg_vertices = rg.findMaxDegreeVertices()
        max_deg_vertex = self.prng.choice(max_deg_vertices)
        
        # Split regions among ssDNA_region, dsDNA_region and other_region
        for edge in rg.edge_list:
            if (edge.v1 == max_deg_vertex or edge.v2 == max_deg_vertex):
                if (edge.doubleStranded):
                    dsDNA_regions.append(edge)
                else:
                    ssDNA_regions.append(edge)
            else:
                unprocessed_regions.append(edge)

        # sampled_structures = {"vertex_label" : (Coordinates(x, y, z), previousDomainInfo=None)}
        # Initial values
        sampled_structures = {str(max_deg_vertex) : (CartesianCoords(0, 0, 0), None)}

        dist = Distributions(self.ssDomainLengthDist, self.dsDomainLengthDist, self.tetherAngleDist, self.ssDomainAngleDist, self.dsdsDomainAngleDist)

        while ((len(dsDNA_regions) > 0) or (len(ssDNA_regions) > 0)  or (len(unprocessed_regions) > 0)): 
            if(len(dsDNA_regions) > 0):
                region_edge = dsDNA_regions.pop(self.prng.randrange(0,len(dsDNA_regions)))
                if (str(region_edge.v1) in sampled_structures.keys() and str(region_edge.v2) in sampled_structures.keys()):
                    skipped_regions.append(region_edge)
                else:
                    sampled_structures, ssDNA_regions, dsDNA_regions, unprocessed_regions = self.sampleJunctionBetweenRegions(sampled_structures, dist, region_edge, ssDNA_regions, dsDNA_regions, unprocessed_regions)
            elif(len(ssDNA_regions) > 0 ):
                region_edge = ssDNA_regions.pop(self.prng.randrange(0,len(ssDNA_regions)))
                if (str(region_edge.v1) in sampled_structures.keys() and str(region_edge.v2) in sampled_structures.keys()):
                    skipped_regions.append(region_edge)
                else:
                    sampled_structures, ssDNA_regions, dsDNA_regions, unprocessed_regions = self.sampleJunctionBetweenRegions(sampled_structures, dist, region_edge, ssDNA_regions, dsDNA_regions, unprocessed_regions)
            else:
                assert False
        return sampled_structures


    def sampleJunctionBetweenRegions(self, sampled_structures, dist, e, ssDNA_regions, dsDNA_regions, unprocessed_regions):
        flag1 = str(e.v1) in sampled_structures.keys() 
        flag2 = str(e.v2) in sampled_structures.keys()  
        assert flag1 or flag2
              
        currentDomain = RegionDomain(e.doubleStranded, e.totalNucleotideLength)
        if(flag1):
            previousCoord = sampled_structures[str(e.v1)][0]
            previousDomainInfo = sampled_structures[str(e.v1)][1]
            label = e.v1
        else:
            previousCoord = sampled_structures[str(e.v2)][0]
            previousDomainInfo = sampled_structures[str(e.v2)][1]
            label = e.v2

        domainUnitVec, domainLengthNM, sampledAngle = samplePoint(previousDomainInfo, currentDomain, dist, self.prng)
        coord = CartesianCoords(previousCoord.x + domainUnitVec.x * domainLengthNM, previousCoord.y + domainUnitVec.y * domainLengthNM, previousCoord.z + domainUnitVec.z * domainLengthNM)
        previousDomainInfo = {}
        previousDomainInfo['unitVec'] = domainUnitVec
        previousDomainInfo['domain'] = currentDomain
        previousDomainInfo['sampledAngle'] = sampledAngle
        previousDomainInfo['prev_label'] = label
        if (flag1):
            sampled_structures[str(e.v2)] = (coord, previousDomainInfo)
        else:
            sampled_structures[str(e.v1)] = (coord, previousDomainInfo)

        # If any regions from other regions to be processed set involve that point,
        # add them to the set of ssDNA/dsDNA regions to be processed next as appropriate.
        del_reg = []
        for rg_edge in unprocessed_regions:
            if(str(rg_edge.v1) in sampled_structures.keys() or str(rg_edge.v2) in sampled_structures.keys()):
                if (rg_edge.doubleStranded):
                    dsDNA_regions.append(rg_edge)
                else:
                    ssDNA_regions.append(rg_edge)
                del_reg.append(rg_edge)
        for item in del_reg:
            unprocessed_regions.remove(item)

        return sampled_structures, ssDNA_regions, dsDNA_regions, unprocessed_regions
    
    # Check whether the constraints are all satisfied simultaneously
    def checkConstraints(self, rg, sampled_structures, debug = False):

        if(self.checkDistanceConstraints(rg.edge_list, sampled_structures) and self.checkAngleConstraints(rg, sampled_structures)):
            if(debug):
                print("Plotting the coordinates from the sampled structures...")
                self.plot_sampled_regiongraph(rg, sampled_structures)
            return True
        else : 
            return False        

    # Check for distance constraints between vertices in the region graph
    def checkDistanceConstraints(self, edge_list, sampled_structures):

        for edges in edge_list:
            c1 = sampled_structures[str(edges.v1)][0]
            c2 = sampled_structures[str(edges.v2)][0]
            d = (((c1.x - c2.x) ** 2) + ((c1.y - c2.y) ** 2) + ((c1.z - c2.z) ** 2)) ** 0.5

            # if single stranded, then inequality equation otherwise equality equation
            if (edges.doubleStranded):
                # Form an equality equation
                l = (edges.totalNucleotideLength * DS_LENGTH)
                if (not math.isclose(d, l)):
                    return False
            else:
                l_ss = edges.totalNucleotideLength * SS_LENGTH
                if (not ((d <= l_ss) or math.isclose(d, l_ss))):
                    return False
        return True

    # Check for angle constraints between vertices in the region graph
    def checkAngleConstraints(self, rg, sampled_strucutres):
        if not NICKED_FLAG: 
            return True
        nicked_angles = rg.computeNickedAngles(sampled_strucutres)
        for key, angles in nicked_angles.items():
            if (angles > NICKEDANGLE_UPPER_BOUND):
                return False   
        return True 

    # For debugging purposes, plot region graph from given sampled coordinates
    def plot_sampled_regiongraph(self, rg, sampled_structures):
        
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(111, projection='3d')
        for edge in rg.edge_list:
            c1 = sampled_structures[str(edge.v1)][0]
            c2 = sampled_structures[str(edge.v2)][0]
            ax.plot([c1.x, c2.x],[c1.y, c2.y ],[c1.z, c2.z], label= str(edge.label))
            length = round(edge.totalNucleotideLength * DS_LENGTH) if edge.doubleStranded else round(edge.totalNucleotideLength * SS_LENGTH)
            ax.text((c1.x +c2.x)/2, (c1.y + c2.y)/2, (c1.z + c2.z)/2, (str(edge.label) + " (nt: " + str(edge.totalNucleotideLength))+", "+ str(length) +")" )       
        plt.legend()
        plt.show()
