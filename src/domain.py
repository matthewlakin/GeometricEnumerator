########################################################################
#
# domain.py
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
# Class for strand graph domains (with labeled bonds)
#

defaultToeholdDomainNucleotideLength = 5
defaultLongDomainNucleotideLength = 20

class Domain(object):

    def __init__(self, name, istoehold, complemented, bond, domainNucleotideLength=None):
        assert isinstance(name, str)
        assert isinstance(istoehold, bool)
        assert isinstance(complemented, bool)
        assert isinstance(bond, str) or isinstance(bond, type(None))
        self.name = name
        self.istoehold = istoehold
        self.complemented = complemented
        self.bond = bond
        if(domainNucleotideLength is None):
            self.domainNucleotideLength = defaultToeholdDomainNucleotideLength if (istoehold) else defaultLongDomainNucleotideLength
        else:
            self.domainNucleotideLength = domainNucleotideLength


    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        output = self.name
        if self.istoehold:
            output += '^'
        if self.complemented:
            output += '*'
        if self.bond is not None:
            output += ('!' + self.bond)
        return output

    def __eq__(self, other):
        return self.__str__() == other.__str__()
    
    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return self.__str__() < other.__str__()

    def __gt__(self, other):
        return self.__str__() > other.__str__()
    
    def getComplement(self):
        return Domain(self.name, self.istoehold, not self.complemented, self.bond, domainNucleotideLength=self.domainNucleotideLength)

    def updateBond(self, newBond):
        return Domain(self.name, self.istoehold, self.complemented, newBond, domainNucleotideLength=self.domainNucleotideLength)

    def stripBond(self):
        return self.updateBond(None)

    def copyDomain(self):
        return Domain(self.name, self.istoehold, self.complemented, self.bond, domainNucleotideLength=self.domainNucleotideLength)

    def isComplementaryTo(self, other):
        return ((self.name == other.name) and
                (self.istoehold == other.istoehold) and
                (self.complemented == (not other.complemented)))

    def wellFormedBondTo(self, other):
        return (self.isComplementaryTo(other) and
                (self.bond is not None) and 
                (self.bond == other.bond))
