########################################################################
#
# distributions.py
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

# The Distributions class encapsulates the different distributions that we
# may need to sample structures:
#   ssDomainLengthDist - distribution of ss domain lengths in nm (given length in nt)
#   dsDomainLengthDist - distribution of ds domain lengths in nm (given length in nt)
#   tetherAngleDist - distribution of angles for domains immediately attached to a tether
#   ssDomainAngleDist - distribution of angles between an ss domain and any other domain
#   dsdsDomainAngleDist - distribution of angles between two ds domains (i.e., a nick)
# Each of these distributions will be represented by a Python object with an appropriate method to call
#   sampleAngle() for angle distributions
#   sampleLengthNm(n) for length distributions, where n is length of domain in nt

class Distributions():

    def __init__(self, ssDomainLengthDist, dsDomainLengthDist, tetherAngleDist,
                 ssDomainAngleDist, dsdsDomainAngleDist):
        self.ssDomainLengthDist = ssDomainLengthDist
        self.dsDomainLengthDist = dsDomainLengthDist
        self.tetherAngleDist = tetherAngleDist
        self.ssDomainAngleDist = ssDomainAngleDist
        self.dsdsDomainAngleDist = dsdsDomainAngleDist

########################################################################
