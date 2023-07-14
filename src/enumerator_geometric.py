########################################################################
#
# enumerator_geometric.py
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

# Reaction enumeration settings and code

import probio_lib as lib
from species import * #speciesFromStrandGraph
from reaction import *
from crn import *
from enumerator_abstract import *

#
############################################################################
# A class that encapsulates the geometric constraints rules
#

class ReactionEnumerator_Geometric(ReactionEnumerator_Abstract):

    ########################################################################
    
    def __init__(self, settings):
        super().__init__()
        self.settings = settings
        self.plausible_species = []
        self.implausible_species = []
        assert self.validSettings()

    ########################################################################
    
    # Check that the settings we provided are valid (wrt types etc)
    def validSettings(self):
        VALID_threeWayModeOptions = ['adjacent']
        VALID_unbindingModeOptions = ['adjacent']
        VALID_enumerationModeOptions = ['detailed']
        VALID_rateOptions = ['bind', 'unbind', 'migrate','displace']
        if sorted(self.settings.keys()) != sorted(['name', 'debug', 'maxComplexSize', 'threeWayMode',
                                                   'unbindingMode', 'enumerationMode', 'rate', 'constraintChecker']):
            print('Settings error: wrong keys: found '+str(self.settings.keys()))
            return False
        if type(self.settings['name']) != str:
            print('Settings error: wrong name option type: found '+str(self.settings['name']))
            return False
        if type(self.settings['debug']) != bool:
            print('Settings error: wrong debug option type: found '+str(self.settings['debug']))
            return False
        if self.settings['enumerationMode'] not in VALID_enumerationModeOptions:
            print('Settings error: wrong enumerationMode option: found '+str(self.settings['enumerationMode'])+' with type '+str(type(self.settings['enumerationMode'])))
            return False
        if type(self.settings['maxComplexSize']) not in [float, int]:
            print('Settings error: wrong maxComplexSize option type: found '+str(self.settings['maxComplexSize'])+' with type '+str(type(self.settings['maxComplexSize'])))
            return False
        if self.settings['threeWayMode'] not in VALID_threeWayModeOptions:
            print('Settings error: illegal option for threeWayMode: found '+str(self.settings['threeWayMode'])+' with type '+str(type(self.settings['threeWayMode'])))
            return False
        if self.settings['unbindingMode'] not in VALID_unbindingModeOptions:
            print('Settings error: illegal option for unbindingMode: found '+str(self.settings['unbindingMode'])+' with type '+str(type(self.settings['unbindingMode'])))
            return False
        if self.settings['constraintChecker'] == None:
            print("Settings error: Constraint Checker object is None")
            return False
        if sorted(self.settings['rate'].keys()) != sorted(VALID_rateOptions):
            print('Settings error: illegal option for rate: found '+str(self.settings['rate'])+' with type '+str(type(self.settings['rate'])))
            return False            
        return True
        
    def debugPrint(self, x, debug=False):
        if(debug):
            print(x)  

    ########################################################################
    
    # 
    # Helper methods that define this enumeration algorithm
    #

    # Compute the set of bound reachable sites via bonds from a given starting site.
    #  * The start site itself is __not__ included in the returned list.
    #  * There are two kinds of traversal step possible: going along a bond or along a strand.
    #  * IN THIS VERSION: the first and last steps MUST be bond traversal steps.
    #  * This means that this traversal will not include any sites on the same strand as the
    #    start site, __unless__ they are connected to the rest of the graph via a bond.
    def boundSitesReachableVersionOne(self, this, startSite):
        assert startSite in this.getSites()
        ReachedSites = [startSite] # Put it in here to prevent looping
        Q = [startSite]
        while Q != []:
            s = Q.pop(0)
            sPrime = this.getBindingPartner(s)
            if sPrime is not None:
                if sPrime not in ReachedSites:
                    ReachedSites += [sPrime]
                    Q += [sPrime]
                if sPrime.v != startSite.v: # Only traverse along the vertex to find other bond start sites if it is NOT the starting vertex!
                    for pns in this.boundSitesOnSameVertexAs(sPrime):
                        if pns not in ReachedSites:
                            ReachedSites += [pns]
                            Q += [pns]
        ReachedSites.remove(startSite) # Don't want to return the starting site in this list
        return ReachedSites

    # Method to check if the structure is plausible 
    def checkPlausibility(self, this):
        cc = self.settings['constraintChecker']
        comps = this.connectedComponents()
        for item in comps:
            for species, sampling_info in self.plausible_species:
                if (species == item):
                   return True 
            for species, sampling_info in self.implausible_species:
                if (species == item):     
                    return False        
            flag, sampling_info = cc.isPlausible(item)
            if (flag):
                self.plausible_species.append((item, sampling_info))
                return True
            else:
                self.implausible_species.append((item, sampling_info))
                return False
        return False

    def allBindingTransitions(self, this):        
        possible_new_edges = this.possibleNewEdges()
        currently_bound_sites = this.currentlyBoundSites()
        all_binding_transitions = []
        for a in possible_new_edges:
            if a.s1 not in currently_bound_sites and a.s2 not in currently_bound_sites :
                edges_added_in_transition = [a]
                edges_removed_in_transition = []
                all_edges_involved_in_transition = sorted(edges_added_in_transition + edges_removed_in_transition)
                new_strand_graph = this.addEdgeToCurrentEdges(a)
                new_strand_graph.domainLength = this.domainLength        
                flag = self.checkPlausibility(new_strand_graph) 
                if(flag):
                    self.debugPrint("checked plausible in binding!!!")
                    this_binding_transition = {'type':'BINDING',
                                            'edges_added':edges_added_in_transition,
                                            'edges_removed':edges_removed_in_transition,
                                            'all_edges_involved':all_edges_involved_in_transition,
                                            'new_strand_graph':new_strand_graph,
                                            'rate': self.settings['rate']['bind']}
                    all_binding_transitions.append(this_binding_transition)
        return all_binding_transitions

    def allUnbindingTransitions(self, this, debug = False):       
        assert this.isConnected()
        all_unbinding_transitions = []
        for e in this.current_edges:
            if e in this.toehold_edges:
                edges_added_in_transition = []
                edges_removed_in_transition = [e]
                all_edges_involved_in_transition = sorted(edges_added_in_transition + edges_removed_in_transition)
                new_strand_graph = this.removeEdgeFromCurrentEdges(e)
                new_strand_graph.domainLength = this.domainLength
                ## NB: Can't use "sameSpecies" method here as that method creates a new graph, and the indexes may change.
                #      Also, new_strand_graph may have different indexes relative to the original.
                #      (Only alternative would be to do a graph traversal on the new strand graph.)
                if(not new_strand_graph.isConnected()):
                    ##################################################################
                    #
                    # ## NB: Commenting this out since just removing an edge should preserve plausibility.
                    # ##     If we assume that the initial strand graph is plausible then the constraints can be satisfied.
                    # ##     Removing an edge just removes some constraints, so any structure that satisfied the constraints
                    # ##     beforehand will also satisfy this reduced set of constraints afterwards.
                    # ##     Hence plausibility is preserved and we can therefore skip rechecking here.
                    #
                    #flag, plausible_species, implausible_species = self.checkPlausibility(new_strand_graph, plausible_species, implausible_species)
                    #if(flag):
                    #    debugPrint("checked plausible in unbinding!!!")
                    #    this_unbinding_transition = {'type':'UNBINDING',
                    #                                 'edges_added':edges_added_in_transition,
                    #                                 'edges_removed':edges_removed_in_transition,
                    #                                 'all_edges_involved':all_edges_involved_in_transition,
                    #                                 'new_strand_graph':new_strand_graph}
                    #    all_unbinding_transitions.append(this_unbinding_transition)
                    ##################################################################
                    this_unbinding_transition = {'type':'UNBINDING',
                                                 'edges_added':edges_added_in_transition,
                                                 'edges_removed':edges_removed_in_transition,
                                                 'all_edges_involved':all_edges_involved_in_transition,
                                                 'new_strand_graph':new_strand_graph,
                                                 'rate': self.settings['rate']['unbind']}
                    all_unbinding_transitions.append(this_unbinding_transition)
        return all_unbinding_transitions
 
    def allThreeWayMigrationTransitions(self, this):
        possible_new_edges = this.possibleNewEdges()
        currently_unbound_sites = this.currentlyUnboundSites()
        all_threeway_migration_transitions = []
        for edge_to_remove in this.current_edges:
            for (s1, s2) in edge_to_remove.bothWaysRound():
                for s in currently_unbound_sites:                
                    edge_to_add = Edge(s, s2)
                    if edge_to_add in possible_new_edges and this.sameSpecies(s, s2):
                        new_strand_graph = this.removeEdgeFromCurrentEdges(edge_to_remove).addEdgeToCurrentEdges(edge_to_add)
                        new_strand_graph.domainLength = this.domainLength
                        flag = self.checkPlausibility(new_strand_graph)      
                        if(flag):                            
                            self.debugPrint("checked plausible in three way migration!!!")
                            edges_added_in_transition = [edge_to_add]
                            edges_removed_in_transition = [edge_to_remove]
                            all_edges_involved_in_transition = sorted(edges_added_in_transition + edges_removed_in_transition)
                            this_threeway_migration_transition = {'type':'THREE_WAY_MIGRATION',
                                                                'edges_added':edges_added_in_transition,
                                                                'edges_removed':edges_removed_in_transition,
                                                                'all_edges_involved':all_edges_involved_in_transition,
                                                                'new_strand_graph':new_strand_graph,
                                                                'rate': self.settings['rate']['displace']}                          
                            all_threeway_migration_transitions.append(this_threeway_migration_transition)
        return all_threeway_migration_transitions

    # Look for 4-way branch migration transitions from "this" species.
    # This function only supports 4-way branch migration, and not n-way migration for n>4.
    # TO DO: maybe generalize to n-way migration?!
    def allFourWayMigrationTransitions(self, this):
        possible_new_edges = this.possibleNewEdges()
        currently_bound_sites = this.currentlyBoundSites()
        all_fourway_migration_transitions = []
        for edge in this.current_edges:
            for (s1,s2) in edge.bothWaysRound():
                self.debugPrint('Testing edge where s1='+str(s1)+' and s2='+str(s2))
                s1pr = s1.threePrimeAdjacentSite()
                s2pr = s2.fivePrimeAdjacentSite()
                self.debugPrint('s1pr = '+str(s1pr))
                self.debugPrint('s2pr = '+str(s2pr))
                if s1pr is not None and s2pr is not None:
                    s3 = this.getBindingPartner(s1pr)
                    self.debugPrint('s3 = '+str(s3))
                    if s3 is not None:
                        s3pr = s3.threePrimeAdjacentSite()
                        self.debugPrint('s3pr = '+str(s3pr))
                        if s3pr is not None:
                            s4 = this.getBindingPartner(s3pr)
                            self.debugPrint('s4 = '+str(s4))
                            if s4 is not None:
                                s4pr = s4.threePrimeAdjacentSite()
                                self.debugPrint('s4pr = '+str(s4pr))
                                if s4pr is not None:
                                    s4pr_binding_partner = this.getBindingPartner(s4pr)
                                    self.debugPrint('s4pr_binding_partner = '+str(s4pr_binding_partner))
                                    if s4pr_binding_partner is not None and s4pr_binding_partner == s2pr:
                                        first_edge_to_add = Edge(s1pr, s2pr)
                                        second_edge_to_add = Edge(s3, s4pr)
                                        self.debugPrint('first_edge_to_add = '+str(first_edge_to_add)+'    ...Possible? '+str(first_edge_to_add in possible_new_edges))
                                        self.debugPrint('second_edge_to_add = '+str(second_edge_to_add)+'    ...Possible? '+str(second_edge_to_add in possible_new_edges))
                                        if first_edge_to_add in possible_new_edges and second_edge_to_add in possible_new_edges:
                                            first_edge_to_remove = Edge(s1pr, s3)
                                            second_edge_to_remove = Edge(s2pr, s4pr)
                                            self.debugPrint('first_edge_to_remove = '+str(first_edge_to_remove))
                                            self.debugPrint('second_edge_to_remove = '+str(second_edge_to_remove))
                                            edges_added_in_transition = sorted([first_edge_to_add, second_edge_to_add])
                                            edges_removed_in_transition = sorted([first_edge_to_remove, second_edge_to_remove])
                                            all_edges_involved_in_transition = sorted(edges_added_in_transition + edges_removed_in_transition)
                                            new_strand_graph = this.removeEdgeFromCurrentEdges(first_edge_to_remove) \
                                                                   .removeEdgeFromCurrentEdges(second_edge_to_remove) \
                                                                   .addEdgeToCurrentEdges(first_edge_to_add) \
                                                                   .addEdgeToCurrentEdges(second_edge_to_add)
                                            if all_edges_involved_in_transition not in [d['all_edges_involved'] for d in all_fourway_migration_transitions]:
                                                new_strand_graph.domainLength = this.domainLength
                                                flag = self.checkPlausibility(new_strand_graph)      
                                                if(flag):
                                                        self.debugPrint("checked plausible in four way migration!!!")
                                                        this_fourway_migration_transition = {'type':'FOUR_WAY_MIGRATION',
                                                                                            'edges_added':edges_added_in_transition,
                                                                                            'edges_removed':edges_removed_in_transition,
                                                                                            'all_edges_involved':all_edges_involved_in_transition,
                                                                                            'new_strand_graph':new_strand_graph,
                                                                                            'rate': self.settings['rate']['displace']}
                                                        self.debugPrint('TRANSITION INFO: '+str(this_fourway_migration_transition))
                                                        all_fourway_migration_transitions.append(this_fourway_migration_transition)
     
        self.debugPrint(all_fourway_migration_transitions)
        return all_fourway_migration_transitions

    # Get all unimolecular transitions possible from "this" species                                                                                                                                         
    def allUnimolecularTransitions(self, this):
        bindingTransition = self.allBindingTransitions(this)
        unbindingTransitions = self.allUnbindingTransitions(this)
        threeWayMigrationTransitions = self.allThreeWayMigrationTransitions(this)
        fourWayMigrationTransitions = self.allFourWayMigrationTransitions(this)
        allTransitions = bindingTransition + unbindingTransitions + threeWayMigrationTransitions + fourWayMigrationTransitions
        return allTransitions

    ########################################################################
    
    # Compute all unimolecular reactions possible starting from "this" species
    def unimolecularReactions(self, this):
        allTransitions = self.allUnimolecularTransitions(this)
        allReactions = []
        reactants = [this]
        for t in allTransitions:
            thisFwdRate = t['rate']
            theseProducts = [speciesFromStrandGraph(sg) for sg in t['new_strand_graph'].connectedComponents()]
            thisMetadata = {'type':t['type'], 'edges_added':t['edges_added'], 'edges_removed':t['edges_removed'], 'all_edges_involved':t['all_edges_involved']}
            thisReaction = Reaction(reactants, thisFwdRate, theseProducts, bwdrate=None,metadata=thisMetadata)
            if thisReaction not in allReactions:
                allReactions += [thisReaction]
        return allReactions

    # Compute all bimolecular reactions possible when "this" species is paired with "that" species
    def bimolecularReactions(self, this, that):
        allTransitions = self.allBindingTransitions(this.compose(that))
        allReactions = []
        reactants = [this, that]
        for t in allTransitions:
            thisFwdRate = t['rate']
            theseProducts = [speciesFromStrandGraph(sg) for sg in t['new_strand_graph'].connectedComponents()]
            thisMetadata = {'type':t['type'], 'edges_added':t['edges_added'], 'edges_removed':t['edges_removed'], 'all_edges_involved':t['all_edges_involved']}
            thisReaction = Reaction(reactants, thisFwdRate, theseProducts, bwdrate=None, metadata=thisMetadata)
            if thisReaction not in allReactions:
                allReactions += [thisReaction]
        return allReactions

    def enumerateReactions(self, species_list):
        assert self.validSettings()                                                                                                                                                                                   
        if not isListOfSpecies(species_list):
            lib.error('In ReactionEnumerator_Original.enumerateReactions: expected list of species as argument, but found: '+str(species_list))
        if not lib.distinct(species_list):
            lib.error('In ReactionEnumerator_Original.enumerateReactions: expected all species in argument list to be unique, but found: '+str(species_list))
        allReactions = []
        species_processed = []
        species_pairs_processed_SORTED = []
        species_to_process = list(species_list)
        self.plausible_species = []
        self.implausible_species = []
        iterationcount = 1                                                                                                                                                                                 
        while (species_to_process != []):
            x = species_to_process.pop(0) # Remove and return first species in the list
            flag_in_plausible_species = self.checkPlausibility(x)
            if (not flag_in_plausible_species):
                continue

            if x.numVertexes() > self.settings['maxComplexSize']:
                lib.error('In enumerateReactions: check for possible polymers! Specified max complex size ('+str(self.settings['maxComplexSize'])+') exceeded by following species: '+str(x))
            #debugPrint('SPECIES X FOR THIS ITERATION:')                                                                                                                                                    
            #debugPrint(x)   
            if self.settings['enumerationMode'] == 'detailed':
                newReactions = self.unimolecularReactions(x)
            elif self.settings['enumerationMode'] == 'infinite':
                assert self.allUnimolecularTransitions(x) == []
                newReactions = []
            else:
                assert False
            for y in species_processed:
                this_sorted_pair = tuple(sorted((x,y)))
                if this_sorted_pair not in species_pairs_processed_SORTED:
                    if self.settings['enumerationMode'] == 'detailed':
                      reacs = self.bimolecularReactions(x, y)
                      newReactions += reacs
                    else:
                        assert False
                    species_pairs_processed_SORTED += [this_sorted_pair]                                                                                                                                                                            
            possiblyNewSpecies = []

            for r in newReactions:
                assert r not in allReactions                                                                                                                                                                 
                allReactions += [r]
                possiblyNewSpecies += r.listOfSpeciesInvolved()

            species_processed += [x] # Do this before the next loop so we don't double-count species!                                                                                                       
            for pns in possiblyNewSpecies:
                if (pns not in species_processed) and (pns not in species_to_process):
                    species_to_process += [pns]
            iterationcount += 1  
        return CRN(species_processed, allReactions)
