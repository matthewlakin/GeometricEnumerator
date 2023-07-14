########################################################################
#
# strand.py
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

from domain import *

class Strand(object):

    def __init__(self, domains):
        assert isinstance(domains, list)
        for d in domains:
            assert isinstance(d, Domain)
        self.domains = domains

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        output = '<'
        for (idx, d) in enumerate(self.domains):
            if idx != 0:
                output += ' '
            output += str(d)
        output += '>'
        return output

    def __eq__(self, other):
        return self.domains == other.domains

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return self.domains < other.domains

    def __gt__(self, other):
        return self.domains > other.domains

    def strandType(self):
        return Strand([d.stripBond() for d in self.domains])

    def copyStrand(self):
        return Strand([d.copyDomain() for d in self.domains])

    def modifyDomain(self, idx, newDomain):
        assert 0 <= idx <= len(self.domains) - 1
        new_domains = [d.copyDomain() for d in self.domains]
        new_domains[idx] = newDomain
        return Strand(new_domains)
