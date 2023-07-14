########################################################################
#
# reaction.py
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

from species import *
from strandgraph import *
import sgparser
import probio_lib as lib

###############################################################################################

class Reaction(object):

    def __init__(self, reactants, fwdrate, products, bwdrate=None, metadata={}):
        self.reactants = sorted(reactants)
        self.fwdrate = fwdrate
        self.bwdrate = bwdrate
        self.products = sorted(products)
        self.metadata = metadata
        assert self.isValid()

    def isValid(self):
        if not isListOfSpecies(self.reactants):
            print('In Reaction.isValid: specified reactants is not a list of Species: found '+str(self.reactants))
            return False
        if not isinstance(self.products, list):
            print('In Reaction.isValid: specified products is not a list of Species: found '+str(self.products))
            return False
        if not isinstance(self.fwdrate, float):
            print('In Reaction.isValid: specified forward rate is not a float: found '+str(self.fwdrate))
            return False
        if self.fwdrate <= 0.0:
            print('In Reaction.isValid: specified forward rate is less than or equal to zero: found '+str(self.fwdrate))
            return False
        if self.bwdrate is not None:
            if not isinstance(self.bwdrate, float):
                print('In Reaction.isValid: specified backward rate is not a float: found '+str(self.bwdrate))
                return False
            if self.bwdrate <= 0.0:
                print('In Reaction.isValid: specified backward rate is less than or equal to zero: found '+str(self.bwdrate))
                return False
        if self.reactants == self.products:
            print('In Reaction.isValid: specified reactants and products are identical: found '+str(self.reactants))
            return False
        if len(self.reactants) > 2:
            print('In Reaction.isValid: more than two reactants: found '+str(self.reactants))
            return False
        if type(self.metadata) is not dict:
            print('In Reaction.isValid: metadata should be a dict: found '+str(self.metadata))
            return False
        return True
        
    def __str__(self):
        return self.__repr__()

    def __repr__(self): # NB: ignoring metadata here
        stringOfReactants = str(self.reactants)
        stringOfRate = '  ->{'+str(self.fwdrate)+'}' if self.bwdrate is None else '  {'+str(self.bwdrate)+'}<->{'+str(self.fwdrate)+'}'
        stringOfProducts = str(self.products)
        return stringOfReactants + os.linesep + stringOfRate + os.linesep + stringOfProducts

    # def prettyPrint(self, speciesNames): # NB: ignoring metadata here
    #     stringOfReactants = prettyPrintSpeciesList(self.reactants, speciesNames)
    #     stringOfRate = ' ->{'+str(self.fwdrate)+'} ' if self.bwdrate is None else ' {'+str(self.bwdrate)+'}<->{'+str(self.fwdrate)+'} '
    #     stringOfProducts = prettyPrintSpeciesList(self.products, speciesNames)
    #     return stringOfReactants + stringOfRate + stringOfProducts
    
    def __metric__(self):
        return (self.reactants, self.fwdrate, -1.0 if self.bwdrate is None else self.bwdrate, self.products)
    
    def __eq__(self, other):
        return self.__metric__() == other.__metric__()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return self.__metric__() < other.__metric__()

    def __gt__(self, other):
        return self.__metric__() > other.__metric__()

    def tryToCombineWith(self, other):
        if (self.reactants == other.reactants) and (self.products == other.products):
            # Found two identically oriented reactions: combine them!
            newFwdRate = self.fwdrate + other.fwdrate
            if self.bwdrate is None and other.bwdRate is None:
                newBwdRate = None
            else:
                newBwdRate = (0.0 if self.fwdrate is None else self.fwdrate) + (0.0 if self.bwdrate is Noe else self.bwdrate)
            newMetadata = dict(self.metadata)  # NB: in future, might want a more nuanced approach to combining metadata when combining reactions
            newMetadata.update(other.metadata) # NB: in future, might want a more nuanced approach to combining metadata when combining reactions
            return Reaction(self.reactants, newFwdRate, self.products, bwdrate=newBwdRate, metadata=newMetadata)
        elif (self.reactants == other.products) and (self.products == other.reactants):
            # Found two oppositely oriented reactions: combine them!
            newFwdRate = self.fwdrate + (0.0 if other.bwdrate is None else other.bwdrate)
            newBwdRate = (0.0 if self.bwdrate is None else self.bwdrate) + other.fwdrate
            newMetadata = dict(self.metadata)  # NB: in future, might want a more nuanced approach to combining metadata when combining reactions
            newMetadata.update(other.metadata) # NB: in future, might want a more nuanced approach to combining metadata when combining reactions
            return Reaction(self.reactants, newFwdRate, self.products, bwdrate=newBwdRate, metadata=newMetadata)
        else:
            return None

    def listOfSpeciesInvolved(self):
        res = []
        for x in self.reactants + self.products:
            if x not in res:
                res += [x]
        return res

    def checkAndGetColorsInfo(self):
        xs = self.listOfSpeciesInvolved()
        assert xs != []
        the_colors_info = xs[0].colors_info
        for x in xs[1:]:
            assert x.colors_info == the_colors_info
        return the_colors_info

    # Make a graphical representation of this reaction!
    # Returns the graphviz object. This allows them to be rendered in a Jupyter notebook via display(), for example.
    # See here for documentation: https://graphviz.readthedocs.io/en/stable/manual.html
    def makeGraphicalRepresentation(self):
        if not GraphvizAvailable:
            lib.error('Graphviz module not available!')
        lib.error('VISUALIZING INDIVIDUAL REACTIONS AS SINGLE GRAPHVIZ OBJECTS IS NOT IMPLEMENTED!!! (YET...)')
        ################################################################################################################################################
        #
        # NB: Giving up on the commented code below, at least for now. Subgraph placement seems rather random and not very attractive.
        #     Can probably get by with text-based reactions plus a key of visualized species, at least for the time being...
        # 
        ################################################################################################################################################
        # def composeSpeciesList(xs):
        #     sgcomp = makeEmptyStrandGraph(self.checkAndGetColorsInfo())
        #     for x in xs:
        #         sgcomp = sgcomp.compose(x)
        #     return sgcomp
        # d = graphviz.Digraph(engine='neato') # 'fdp' -- maybe better for implementing clusters etc?
        # d.attr(colorscheme='X11')
        # #d.attr(nodesep=1.0)
        # d.attr(rankdir='TB') # LR
        # d.attr(ranksep='1.0')
        # d.attr(compound='true')
        # with d.subgraph(name='cluster_reactants') as cl:
        #     reactantNodeNames = composeSpeciesList(self.reactants).__fillInGraphicalRepresentation__(cl, prefix='reactants_', outlinecolor='blue')
        # with d.subgraph(name='cluster_products') as cl:
        #     productNodeNames = composeSpeciesList(self.products).__fillInGraphicalRepresentation__(cl, prefix='products_', outlinecolor='green')
        # assert reactantNodeNames != []
        # assert productNodeNames != []
        # d.edge(reactantNodeNames[0], productNodeNames[1], ltail='cluster_reactants', lhead='cluster_products', label=str(self.fwdrate))
        # return d
        ################################################################################################################################################

    # Display graphical representation of this reaction!
    # Tries to check if display() is bound, so if we are in the terminal, so if we are in a Jupyter notebook, this should draw the visualization.
    # If we are running from the terminal, this should not do anything.
    # NB: for now, am hacking this up from representations of individual strand graphs.
    # def displayRepresentation(self):
    #     try:
    #         display(self.makeGraphicalRepresentation())
    #     except NameError: # Should jump back to here if not in Jupyter, as display won't be defined (by default)
    #         pass
    def displayRepresentation(self): # NB: ignoring metadata here
        try:
            testJupyter = display
            useGraphviz = True
        except NameError: # Should jump back to here if not in Jupyter, as display won't be defined by default. NB: any other NameError will be picked up too, however...
            useGraphviz = False
        if useGraphviz:
            def displayRepresentationOfSpeciesList(xs):
                sgcomp = makeEmptyStrandGraph(self.checkAndGetColorsInfo())
                for x in xs:
                    sgcomp = sgcomp.compose(x)
                sgcomp.displayRepresentation()
            spacer = '==================================================================='
            print(spacer)
            print('REACTANTS:')
            displayRepresentationOfSpeciesList(self.reactants)
            print('RATE:')
            if self.bwdrate is None:
                print(str(self.fwdrate))
            else:
                print('Forward: '+str(self.fwdrate))
                print('Backward: '+str(self.bwdrate))
            print('PRODUCTS:')
            displayRepresentationOfSpeciesList(self.products)
            print(spacer)
        else:
            print(self)
    
###############################################################################################
