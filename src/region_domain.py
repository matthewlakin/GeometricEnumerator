########################################################################
#
# region_domain.py
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

########################################################################

# Basic representation of individual domains: just ssDNA vs dsDNA and length in nucleotides
class RegionDomain:
    def __init__(self, isDS, lengthNT):
        self.isDS = isDS
        self.lengthNT = lengthNT
    
    def __repr__(self):
        return ('DS' if self.isDS else 'SS') + str(self.lengthNT)

    def __str__(self):
        return self.__repr__()

    def maxLength(self):
        return self.lengthNT * (DS_LENGTH if self.isDS else SS_LENGTH)

    def persistenceLength(self):
        return (DSDNA_PERSISTENCE_LENGTH
                if self.isDS
                else SSDNA_PERSISTENCE_LENGTH)

########################################################################
