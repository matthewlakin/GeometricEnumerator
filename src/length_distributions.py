########################################################################
#
# length_distributions.py
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

from constants import *
import math
import numpy as np
import scipy.integrate

########################################################################

class UniformLengthDistribution:

    def sampleLengthNm(self, domain, prng):
        res = prng.uniform(0, domain.maxLength())
        assert 0.0 <= res <= domain.maxLength()
        return res

########################################################################

class MaxLengthDistribution:

    def sampleLengthNm(self, domain, prng):
        res = domain.maxLength()
        #print(res)        
        assert 0.0 <= res <= domain.maxLength()
        return res

########################################################################

class WormLikeChainLengthDistribution:

    def __init__(self):
        # This is a cache for precomputed WLC cumulative probability distributions
        self.__wlc_cache__ = {}
        
    # Compute probability density of worm-like chain (WLC) of length L with
    # persistence length s being extended to length R.
    # R, s, and L should all have the same length unit.
    def __wlc_prob_density__(self, R, s, L):
        if R > L:
            raise ValueError('R > L')
        elif R == L:
            raise ValueError('R == L')
        R = float(R)
        s = float(s)
        L = float(L)
        t = L / s
        r = R / L
        A = ((4.0 * ((0.75*t)**1.5) * math.exp(0.75*t)) /
             ((math.pi**1.5) * (4.0 + 12/(0.75*t) + 15/((0.75*t)**2))))
        res = (1.0 / L) * ((4.0*math.pi*A*(r**2)*math.exp(-0.75*t/(1.0-(r**2)))) /
                           ((1.0-r**2)**4.5))
        return res

    # This method serves as an interface to get WLC cumulative probability distributions
    # If nothing is stored in the cache for the supplied argument combination,
    # the cumulative probabilities are computed and stored in the cache.
    # If there is an entry in the cache for the supplied argument combination,
    # it is simply returned.
    def __get_wlc_cumul_probs__(self, s, L, num_slices):
        if (s,L,num_slices) not in self.__wlc_cache__:
            these_xs = []
            these_probs = []
            delta = L / num_slices
            x = delta
            while True:
                cumul_prob = scipy.integrate.quad(lambda R: self.__wlc_prob_density__(R,s,L), 0.0, x)[0]
                these_xs += [x]
                these_probs += [cumul_prob]
                if x >= L:
                    break
                x = L if x+delta > L else x+delta
            self.__wlc_cache__[(s,L,num_slices)] = (list(these_xs), list(these_probs))
        return self.__wlc_cache__[(s,L,num_slices)]

    # Sample a single value from the WLC distribution.
    # The persistence length s and maximum length L are derived from the supplied domain information.
    # s and L should both have the same length unit.
    # Assume that the distribution has been approximated via discretization
    # into num_slices slices, which is 1000 by default.
    def sampleLengthNm(self, domain, prng, num_slices=1000):
        s = domain.persistenceLength()
        L = domain.maxLength()
        (xs, cumul_probs) = self.__get_wlc_cumul_probs__(s, L, num_slices=num_slices)
        pGenerated = float(prng.random())
        res = np.interp(pGenerated, cumul_probs, xs)
        assert 0.0 <= res <= domain.maxLength()
        return res

########################################################################
