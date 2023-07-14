########################################################################
#
# species.py
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
from strandgraph import *
import probio_lib as lib
import sgparser
import os
import math
   
###############################################################################################

#
# Species class (subclass of StrandGraph)
#
# - Represents a specific subset of strand graphs that are connected, and therefore represent a single species.
# - Currently implemented as a subclass of StrandGraph.
# - We could do away with this by generalizing graph enumeration to non-connected strand graphs...

class Species(StrandGraph):

    # Initializer essentially just checks that the species is connected
    def __init__(self, colors_info, vertex_colors, admissible_edges, toehold_edges, current_edges, domainLength):
        super().__init__(colors_info, vertex_colors, admissible_edges, toehold_edges, current_edges, domainLength)
        if self.isConnected():
            self.__convertToCanonicalForm__()
        else:
            errMsg = 'Tried to create a Species object from the following non-connected strand graph:'+os.linesep+str(self)
            lib.error(errMsg)
            #print('ERROR: '+str(errMsg)) # Commented this out for testing purposes. Ultimately want to crash if this happens!

#
# Additional helper functions for species
#

# Check whether something is a list of species
def isListOfSpecies(xs):
    if not isinstance(xs, list):
        return False
    for x in xs:
        if not isinstance(x, Species):
            return False
    return True

# Given a connected strand graph, convert it into a Species object
def speciesFromStrandGraph(sg):
    assert sg.isConnected()
    return Species(sg.colors_info, sg.vertex_colors, sg.admissible_edges, sg.toehold_edges, sg.current_edges, sg.domainLength)

# Given a process, convert it into a list of species (NB: there may be some duplicates?!)
def speciesListFromProcess(p):
    assert isinstance(p, Process)
    assert p.wellFormed()
    return [speciesFromStrandGraph(sg) for sg in connectedStrandGraphsFromProcess(p)]

###############################################################################################

