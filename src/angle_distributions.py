########################################################################
#
# angle_distributions.py
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
# NOTES ON SAMPLING THE ANGLE
# ===========================
#
# See this link on sphere point picking: http://mathworld.wolfram.com/SpherePointPicking.html
# This link validates our approach in structures.py of selecting \theta uniformly between 0 and 2*\pi.
# The only question then is how to sample \phi, the angle which sets the "steepness" of the cone.
#
# For the nicked distribution, we just draw from the histogram since that is non-uniform anyway,
# and convert to radians before we return.
#
# For the distributions that are intended to be uniform, we must compensate as per the link above.
#
# For the uniformly distributed full-sphere:
#  Choose v uniformly distributed between 0 and 1
#  Then, \phi = \arccos(2*v - 1)
#
# For the uniformly distriuted hemi-sphere:
#  Choose v uniformly distributed between 0 and 0.5
#  Then, \phi = \arccos(2*v - 1)
#

########################################################################

#from constants import *
import math
import numpy as np
from scipy.stats import gaussian_kde

########################################################################

class UniformSphereAngleDistribution:

    def sampleAngle(self, prng):
        v = prng.random() # Uniformly distributed between 0 and 1
        phi = math.acos(2*v - 1) # Circle compensation
        return phi

########################################################################

class UniformHemisphereAngleDistribution:

    def sampleAngle(self, prng):
        v = prng.uniform(0.0, 0.5) # Uniformly distributed between 0 and 0.5
        phi = math.acos(2*v - 1) # Circle compensation
        return phi

########################################################################

class NickedAngleDistribution:

    def __init__(self):
        # nicked_angles_data.tsv is raw data from Chatterjee et al SI.
        # First column is true nick data.
        # Second column is with one or more based deleted from the nick.
        # So we only really care about the first column.
        angles = []
        with open('nicked_angles_data.tsv') as f:
            next(f)
            for line in f:
                data = line.split("\t")
                angles.append(float(data[1]))
        self.angles = np.array(angles)
        self.x_grid = np.linspace(min(self.angles), max(self.angles), len(self.angles))
        self.kdepdf = self.kde(self.angles, self.x_grid, bandwith=0.1)
        self.cdf = self.initalize_cdf()

    def kde(self, x,x_grid,bandwith=0.2, **kwargs):
        # Kernel Density Estimation with scipy
        kde = gaussian_kde(x,bw_method=bandwith/x.std(ddof=1),**kwargs)
        return kde.evaluate(x_grid)

    def initalize_cdf(self):
        cdf = np.cumsum(self.kdepdf)
        cdf = cdf / cdf[-1]
        return cdf

    def sampleAngle(self, prng):
        value = prng.random()
        value_bins = np.searchsorted(self.cdf, value)
        degrees_random_from_cdf = self.x_grid[value_bins]
        return math.radians(degrees_random_from_cdf)

########################################################################
