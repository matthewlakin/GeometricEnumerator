########################################################################
#
# process.py
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

from strand import *
import os

class Process(object):

    def __init__(self, strands):
        assert isinstance(strands, list)
        for s in strands:
            assert isinstance(s, Strand)
        self.strands = strands

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return self.compactString()

    # (Maybe) compact string representation
    def compactString(self, useNewlines=True):
        output = ''
        for (idx,s) in enumerate(self.strands):
            if idx == 0:
                output += '( '
            else:
                output += '| '
            output += str(s)
            if idx == len(self.strands) - 1:
                output += ' )'
            else:
                if useNewlines:
                    output += os.linesep
                else:
                    output += ' '
        return output

    # We assume that all processes P are well-formed in that each
    # bond i in P appears exactly twice and is shared between complementary domains.
    def wellFormed(self):
        bonds_dict = {}
        for s in self.strands:
            for d in s.domains:
                if d.bond is not None:
                    bondname = d.bond
                    if bondname not in bonds_dict:
                        bonds_dict[bondname] = [d]
                    else:
                        bonds_dict[bondname] += [d]
        for (bondname, domains) in bonds_dict.items():
            if (len(domains) == 2) and (domains[0].wellFormedBondTo(domains[1])):
                pass
            else:
                return False
        return True
