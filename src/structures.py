########################################################################
#
# structures.py
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

#
# NOTES ON PICKING NEXT UNIT VECTOR GIVEN SAMPLED ANGLE
# =====================================================
#
# These notes pertain to the implementation of makeNextUnitVec below.
#
# URL to helpful page: https://math.stackexchange.com/questions/73237/parametric-equation-of-a-circle-in-3d-space
#
# Given v, our previous unit vector, we need to use one of the methods defined below
# to randomly choose the offset angle of the next unit vector relative to that one.
#
# This is just a single angle, and we assume that it could be in any direction
# therefore forming a cone centered on the previous unit vector whose "steepness"
# will be determined by the sampled angle.
#
# To turn the sampled angle into 3D we must use it to form the
# parameteric equations of the circle representing the rim of the cone,
# and then sample the parametric angle to give us a point on that circle.
# From this point we can then back out the new unit vector.
#
# Given previous unit vector v = (v1,v2,v3), assume v starts at the origin (0,0,0).
# Thus the "far end" of v is at (v1,v2,v3).
# Let the sampled angle be \phi.
# Then the centre of the new circle is c = (c1,c2,c3) where:
#  c1 = v1 + v1*cos(\phi)
#  c2 = v2 + v2*cos(\phi)
#  c3 = v3 + v3*cos(\phi)
# Also, the radius of the circle is r = sin(\phi).
#
# Now we set up the parametric equations of the circle by finding two
# unit vectors a and b which are perpendicular to each other and to v.
# These will serve as basis vectors for our parametric equations.
# To do this, we find a suitable vector for a, and then use
# a cross product of that vector with v to get the other one, b.
#
# The unit vector a must be such that a \cdot v = 0.
# In other words, (a1 \cdot v1) + (a2 \cdot v2) + (a3 \cdot v3) = 0.
# We note that at least one element of v must be non-zero,
# as the unit vector can't have length zero.
# We let vi be the non-zero element and let vj,vk be the others.
# Then we fix aj,ak (e.g., aj = ak = 1) and solve for ai using:
#  ai = (-aj*vj - ak*vk)/vi
# This will always work since we picked vi to be non-zero.
# Finally we normalize a = (a1,a2,a3) to be a unit vector.
#
# Then we can compute b = a \cross{} v.
# There is a numpy.cross function that could be used for this.
# See here for an explicit formula: https://mathinsight.org/cross_product_formula
# It should be easy enough as we have a and b wrt the orthonormal
# basis of the x,y,z,-axes.
#
# Then, the parametric equation of the circle that represents
# the rim of our cone is given by:
#  x(\theta) = c1 + r*cos(\theta)*a1 + r*sin(\theta)*b1
#  y(\theta) = c2 + r*cos(\theta)*a2 + r*sin(\theta)*b2
#  z(\theta) = c3 + r*cos(\theta)*a3 + r*sin(\theta)*b3
#
# Now we just need to sample a uniformly distributed value for
# \theta between 0 and 2*\pi radians, and plug that into the
# parametric equations above to produce a uniformly-sampled
# point p = (p1,p2,p3) on the circle.
# 
# To get the next unit vector u = (u1,u2,u3),
# we simply subtract v from p, so that:
#  u1 = p1 - v1
#  u2 = p2 - v2
#  u3 = p3 - v3
# Finally, we normalize u to be a unit vector, as required.
#

########################################################################

#from constants import *
#from distributions import *
from region_domain import *
#from angle_distributions import *
#from length_distributions import *
import numpy
import math
import os

########################################################################

# Simple class for representing cartesian coordinates in 3D space
class CartesianCoords:
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __repr__(self):
        return('('+str(self.x) + "," + str(self.y) + "," + str(self.z)+')')

    def __str__(self):
        return self.__repr__()

    def __hash__(self):
        return hash((self.x, self.y, self.z))

    def __eq__(self, other):
        return ((self.x == other.x) and (self.y == other.y) and (self.z == other.z))

    def __ne__(self, other):
        return not self.__eq__(other)

    def print_xyz(self):
        return(str(self.x) + "," + str(self.y) + "," + str(self.z))

########################################################################

# # Simple class for representing a UNIT vector in 3D space
class UnitVector:
    def __init__(self, x, y, z):
        magnitude = math.sqrt((x**2) + (y**2) + (z**2))
        self.x = float(x) / float(magnitude)
        self.y = float(y) / float(magnitude)
        self.z = float(z) / float(magnitude)
 
    def __repr__(self):
        return('['+str(self.x) + "," + str(self.y) + "," + str(self.z)+']')
 
    def __str__(self):
        return self.__repr__()
 
    def __hash__(self):
        return hash((self.x, self.y, self.z))
 
    def __eq__(self, other):
        return ((self.x == other.x) and (self.y == other.y) and (self.z == other.z))

    def __ne__(self, other):
        return not self.__eq__(other)

    def print_xyz(self):
        return(str(self.x) + "," + str(self.y) + "," + str(self.z))

########################################################################

# Abstracted representation of a "linear" structure with a single terminal tether
# We don't worry about issues like 5' to 3' orientation etc (yet!)
class AbstractLinearStructure:
    def __init__(self, tetherCoords, domainsList):
        self.tetherCoords = tetherCoords
        self.domainsList = domainsList

    def __repr__(self):
        res = '{tether@'+str(self.tetherCoords)+'}'
        for domain in self.domainsList:
            res += '{'+str(domain)+'}'
        return res
        
    def __str__(self):
        return self.__repr__()

########################################################################

# Function to use sampled deviation angle and previous unit vector
# to create a new unit vector that deviates by that angle,
# in a randomly chosen direction.
def makeNextUnitVec(previousUnitVec, sampledAngle, prng):

    #Center of the circle
    center_of_circle = [previousUnitVec.x + previousUnitVec.x * math.cos(sampledAngle), 
                          previousUnitVec.y + previousUnitVec.y * math.cos(sampledAngle),
                          previousUnitVec.z + previousUnitVec.z * math.cos(sampledAngle)]

    #Finding two basis axis a and b
    if previousUnitVec.x != 0:
        a2 = 1
        a3 = 1
        a1 = ((-1) * (a2 * previousUnitVec.y + a3 * previousUnitVec.z)/ previousUnitVec.x)
    elif previousUnitVec.y != 0:
        a1 = 1
        a3 = 1
        a2 = ((-1) * (a1 * previousUnitVec.x + a3 * previousUnitVec.z)/ previousUnitVec.y)
    elif previousUnitVec.z != 0:
        a1 = 1
        a2 = 1
        a3 = ((-1) * (a1 * previousUnitVec.x + a2 * previousUnitVec.y)/ previousUnitVec.z)
    else:
        assert False

    basis_vector_a = UnitVector(a1, a2, a3)
    basis_vector_b_numpy = numpy.cross([basis_vector_a.x, basis_vector_a.y, basis_vector_a.z], [previousUnitVec.x, previousUnitVec.y, previousUnitVec.z]) 
    basis_vector_b = UnitVector(basis_vector_b_numpy[0], basis_vector_b_numpy[1], basis_vector_b_numpy[2])

    theta = prng.uniform(0, 2 * math.pi)
    radius = math.sin(sampledAngle)
 
    point_on_circle = [center_of_circle[0] + radius * math.cos(theta)* basis_vector_a.x + radius * math.sin(theta) * basis_vector_b.x,
                       center_of_circle[1] + radius * math.cos(theta)* basis_vector_a.y + radius * math.sin(theta) * basis_vector_b.y,
                       center_of_circle[2] + radius * math.cos(theta)* basis_vector_a.z + radius * math.sin(theta) * basis_vector_b.z]

    u1 = point_on_circle[0] - previousUnitVec.x
    u2 = point_on_circle[1] - previousUnitVec.y
    u3 = point_on_circle[2] - previousUnitVec.z
    
    nextUnitVec = UnitVector(u1, u2, u3)
    return nextUnitVec, sampledAngle


# Function to sample the next unit vector given the previous one.
# The angle distribution to use is distributions.ssDomainAngleDist
# if either of the current and previous domains is single-stranded.
# If both of the current and previous domains are double-stranded,
# then distributions.dsdsDomainAngleDist should be used instead.
def sampleNextUnitVec(previousDomainInfo, currentDomain, distributions, prng):
    angleDistToUse = (distributions.dsdsDomainAngleDist
                      if previousDomainInfo['domain'].isDS and currentDomain.isDS
                      else distributions.ssDomainAngleDist)
    sampledAngle = angleDistToUse.sampleAngle(prng)
    return makeNextUnitVec(previousDomainInfo['unitVec'], sampledAngle, prng)

# Function to sample the initial unit vector from a tether.
# The angle distribution to use should be distributions.tetherAngleDist
def sampleInitialUnitVec(distributions, prng):
    sampledAngle = distributions.tetherAngleDist.sampleAngle(prng)
    dummyPreviousUnitVec = UnitVector(0,0,1) ## X=0, Y=0, Z=1
    return makeNextUnitVec(dummyPreviousUnitVec, sampledAngle, prng)

# Function to sample the length of a domain.
# The length distribution will be either distributions.ssDomainLengthDist
# or distributions.dsDomainLengthDist, depending on whether the domain
# is single-stranded or double-stranded.
# Both of these will take as an argument the length of the domain in nt.
def sampleDomainLength(domain, distributions, prng):
    lengthDistToUse = (distributions.dsDomainLengthDist
                          if domain.isDS
                       else distributions.ssDomainLengthDist)
    return lengthDistToUse.sampleLengthNm(domain, prng)

# Main function to generate sampled domain info for a structure,
# given a distributions object (see above).
def samplePoint(previousDomainInfo, currentDomain, distributions, prng):

    if previousDomainInfo is None:
        thisDomainUnitVec, sampledAngle = sampleInitialUnitVec(distributions, prng)
    else:
        thisDomainUnitVec, sampledAngle = sampleNextUnitVec(previousDomainInfo, currentDomain, distributions, prng)
    
    thisDomainLengthNm = sampleDomainLength(currentDomain, distributions, prng)
     
    return (thisDomainUnitVec, thisDomainLengthNm, sampledAngle)

def sampleStructure(absLinStruct, distributions):
    domainUnitVecs = []
    domainLengthsNm = []
    firstCoord = absLinStruct.tetherCoords
    previousDomainInfo = None
    #print(absLinStruct.domainsList)
    for currentDomain in absLinStruct.domainsList:
        if previousDomainInfo is None:
            thisDomainUnitVec = sampleInitialUnitVec(distributions)
            previousDomainInfo = {}
        else:
            thisDomainUnitVec = sampleNextUnitVec(previousDomainInfo, currentDomain, distributions)
        
        #print(currentDomain)
        thisDomainLengthNm = sampleDomainLength(currentDomain, distributions)
                    
        domainUnitVecs += [thisDomainUnitVec]
                    
        domainLengthsNm += [thisDomainLengthNm]
        #print("domainLengthsNm")
        #print(domainLengthsNm)
                    
                    
        previousDomainInfo['unitVec'] = thisDomainUnitVec
        previousDomainInfo['domain'] = currentDomain

    """print("domainUnitVec")
    print(domainUnitVecs)
    print("domainLengthsNm")
    print(domainLengthsNm)
    """     
            
    return (domainUnitVecs, domainLengthsNm)


# Class to store the result of sampling a physical conformation for an abstracted structure.
class SampledLinearStructure:
    def __init__(self, absLinStruct, distributions):
        (domainUnitVecsNm,domainLengthsNm) = sampleStructure(absLinStruct, distributions)
        self.tetherCoords = absLinStruct.tetherCoords
        self.domainUnitVecs = domainUnitVecsNm
        self.domainLengthsNm = domainLengthsNm
        # Figure out the joint coords and midpoints and store them
        previousCoords = CartesianCoords(self.tetherCoords.x, self.tetherCoords.y, self.tetherCoords.z)
        self.jointCoords = []
        self.domainMidpointCoords = []
        for (v,l) in zip(self.domainUnitVecs, self.domainLengthsNm):
            newX = previousCoords.x + v.x * l
            newY = previousCoords.y + v.y * l
            newZ = previousCoords.z + v.z * l
            newCoords = CartesianCoords(newX, newY, newZ)
            self.jointCoords += [newCoords]
            midpointX = previousCoords.x + v.x * 0.5 * l
            midpointY = previousCoords.y + v.y * 0.5 * l
            midpointZ = previousCoords.z + v.z * 0.5 * l
            midpointCoords = CartesianCoords(midpointX, midpointY, midpointZ)
            self.domainMidpointCoords += [midpointCoords]
            previousCoords = newCoords

    def __repr__(self):
        res = '<tether@'+str(self.tetherCoords)+'>'
        for (domainUnitVec,domainLengthNm) in zip(self.domainUnitVecs, self.domainLengthsNm):
            res += '<'+str(domainLengthNm)+'toward'+str(domainUnitVec)+'>'
        return res

    def __str__(self):
        return self.__repr__()

    def getJointCoords(self, index):
        return self.jointCoords[index]

    def getDomainMidpointCoords(self, index):
        return self.domainMidpointCoords[index]

    def getLastDomainMidpointCoords(self):
        return self.getDomainMidpointCoords(-1)

    def writeToFile(self, file_name):
        res = ''
        #for (domainUnitVec,domainLengthNm) in zip(self.domainUnitVecs, self.domainLengthsNm):
        #    res += domainUnitVec.print_xyz() + os.linesep
        for vec in self.jointCoords:
            res += vec.print_xyz() + os.linesep
        with open(file_name, 'a') as out_file:
            out_file.write(res+os.linesep)

    def writeToFileLastMidpointCoords(self, file_name, m):
        res = ''
        res += self.getLastDomainMidpointCoords().print_xyz()
        with open(file_name, m) as out_file:
            out_file.write(res+os.linesep)

########################################################################
