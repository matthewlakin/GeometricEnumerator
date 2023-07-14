########################################################################
#
# enumerator_original.py
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

import math
import probio_lib as lib
from species import speciesFromStrandGraph
from reaction import *
from crn import *
from enumerator_abstract import *

############################################################################

# 
# A class that encapsulates the "original" reaction enumeration algorithms (several variants)
# 

class ReactionEnumerator_Original(ReactionEnumerator_Abstract):

    ########################################################################

    def __init__(self, settings):
        super().__init__()
        self.settings = settings
        assert self.validSettings()

    ########################################################################

    # Check that the settings we provided are valid (wrt types etc)
    def validSettings(self):
        VALID_threeWayModeOptions = ['adjacent', 'anchored_strgsd'] ## , 'anchored_traversal']
        VALID_unbindingModeOptions = ['adjacent', 'anchored_strgsd'] ## , 'anchored_traversal']
        VALID_enumerationModeOptions = ['detailed', 'infinite']
        if sorted(self.settings.keys()) != sorted(['name', 'debug', 'maxComplexSize', 'threeWayMode', 'unbindingMode', 'enumerationMode']):
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
        return True

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

    # Implement a "hidden" predicate similar to (same as?) that from strgsd.
    # This uses the boundSitesReachableVersionOne() method from above to look for graph
    # traversals that enclose either end of the specified edge.
    def hidden(self, this, e):
        assert e not in this.current_edges # Should only be checking a prospective new edge here
        (s1, s2) = e.getSites()
        # Special case: check if both sites are on the same strand and have no bonds in between - forming such edges is OK (i.e., the edge is NOT hidden and the function returns False).
        if e.withinOneStrand(): #if s1.v == s2.v:
            numInterveningSitesBound = 0
            for interveningSite in s1.interveningSitesOnSameVertex(s2):
                if this.siteIsBound(interveningSite):
                    numInterveningSitesBound += 1
            if numInterveningSitesBound == 0:
                return False
        for edgeSite in (s1, s2):
            fivePrimeSites = this.boundSitesFivePrimeFrom(edgeSite)
            threePrimeSites = this.boundSitesThreePrimeFrom(edgeSite)
            for startSite in fivePrimeSites:
                reachableSites = self.boundSitesReachableVersionOne(this, startSite)
                for possibleEndSite in threePrimeSites:
                    if possibleEndSite in reachableSites:
                        return True
        return False

    # Definition of anchored predicate from strgsd (searching in 3' direction only)
    def anchored_3prime_strgsd(self, this, startSite, targetSite):
        currentSite = startSite
        while True:
            threePrimeSite = currentSite.threePrimeAdjacentSite()
            if threePrimeSite is None:
                return False
            else:
                if threePrimeSite == targetSite:
                    return True
                else:
                    threePrimeSiteBinder = this.getBindingPartner(threePrimeSite)
                    if threePrimeSiteBinder is None:
                        return False
                    else:
                        if threePrimeSiteBinder is targetSite:
                            return True
                        else:
                            currentSite = threePrimeSiteBinder
        return False

    # # Compute the set of bound reachable sites via bonds from a given starting site.
    # #  * The start site itself is __not__ included in the returned list.
    # #  * There are two kinds of traversal step possible: going along a bond or along a strand.
    # #  * IN THIS VERSION: the first and last steps MUST be strand traversal steps.
    # #  * This means that this traversal will not include any site bound to the start site,
    # #    __unless__ there is some other route between those two sites.
    # #  * This is used for the anchored_3prime_traversal test (currently disabled).
    # def boundSitesReachableVersionTwo(self, this, startSite):
    #     assert startSite in this.getSites()
    #     ReachedSites = [startSite] # Put it in here to prevent looping
    #     Q = [startSite]
    #     while Q != []:
    #         s = Q.pop(0)
    #         for pns in this.boundSitesOnSameVertexAs(s):
    #             if pns not in ReachedSites:
    #                 ReachedSites += [pns]
    #                 Q += [pns]
    #         if s != startSite:
    #             sPrime = this.getBindingPartner(s)
    #             if sPrime is not None:
    #                 if sPrime not in ReachedSites:
    #                     ReachedSites += [sPrime]
    #                     Q += [sPrime]
    #     ReachedSites.remove(startSite) # Don't want to return the starting site in this list
    #     return ReachedSites

    # # Definition of more general anchored using graph traversal (searching in 3' direction only)
    # def anchored_3prime_traversal(self, this, startSite, targetSite):
    #     return targetSite in self.boundSitesReachableVersionTwo(this, startSite)

    # Anchored predicate: can pick from two different versions using the "unbindingMode" setting ('anchored_strsgd' or 'anchored_traversal')
    # The 'anchored_strgsd' option produces the same behavior as in our TCS paper on strand graphs.
    # The 'anchored_traversal' option uses graph traversal to implement a more permissive notion of being 'anchored'.
    def anchored(self, this, e):
        # Try the sites both way round. This means it's OK to only consider searching in the 3' direction (without loss of generality)
        for (startSite, targetSite) in e.bothWaysRound():
            assert startSite != targetSite
            if self.settings['unbindingMode'] == 'anchored_strgsd':
                if self.anchored_3prime_strgsd(this, startSite, targetSite): # This version uses the strgsd definition of anchored predicate
                    return True
            # elif self.settings['unbindingMode'] == 'anchored_traversal':
            #     if self.anchored_3prime_traversal(this, startSite, targetSite): # This version uses a more general definition of the anchored predicate, via graph traversal
            #         return True
            else:
                assert False
        return False

    ########################################################################
    
    # 
    # Key methods for defining valid transitions between strand graphs under this enumeration algorithm
    #

    def allBindingTransitions(self, this):
        possible_new_edges = this.possibleNewEdges()
        currently_bound_sites = this.currentlyBoundSites()
        all_binding_transitions = []
        for a in possible_new_edges:
            # hidden predicate checks for binding in hairpin loops etc. Is it too restrictive???
            if a.s1 not in currently_bound_sites and a.s2 not in currently_bound_sites and not self.hidden(this, a):
                edges_added_in_transition = [a]
                edges_removed_in_transition = []
                all_edges_involved_in_transition = sorted(edges_added_in_transition + edges_removed_in_transition)
                new_strand_graph = this.addEdgeToCurrentEdges(a)
                this_binding_transition = {'type':'BINDING',
                                           'edges_added':edges_added_in_transition,
                                           'edges_removed':edges_removed_in_transition,
                                           'all_edges_involved':all_edges_involved_in_transition,
                                           'new_strand_graph':new_strand_graph}
                all_binding_transitions.append(this_binding_transition)
        return all_binding_transitions

    ############################################################################################################

    def allUnbindingTransitions(self, this):
        all_unbinding_transitions = []
        for e in this.current_edges:
            if e in this.toehold_edges:
            #if e in this.toehold_edges and not self.anchored(this, e): # Using the strgsd definition of anchored predicate, at least for now...
            ##if e in this.toehold_edges and not this.has_adjacent(e): # A bit more permissive than commented test
                if self.settings['unbindingMode'] == 'adjacent':
                    canUnbind = not this.has_adjacent(e)
                elif self.settings['unbindingMode'] in ['anchored_strgsd']: ## , 'anchored_traversal']:
                    canUnbind = not self.anchored(this, e)
                else:
                    assert False
                if canUnbind:
                    edges_added_in_transition = []
                    edges_removed_in_transition = [e]
                    all_edges_involved_in_transition = sorted(edges_added_in_transition + edges_removed_in_transition)
                    new_strand_graph = this.removeEdgeFromCurrentEdges(e)
                    this_unbinding_transition = {'type':'UNBINDING',
                                                 'edges_added':edges_added_in_transition,
                                                 'edges_removed':edges_removed_in_transition,
                                                 'all_edges_involved':all_edges_involved_in_transition,
                                                 'new_strand_graph':new_strand_graph}
                    all_unbinding_transitions.append(this_unbinding_transition)
        return all_unbinding_transitions

    # This "adjacent" version only does 3-way reactions where the invader is immediately adjacent to the incumbent.
    # It does not permit 3-way initiated 4-way migration, or remote toehold-style reactions.
    def allThreeWayMigrationTransitions_Adjacent(self, this):
        possible_new_edges = this.possibleNewEdges()
        currently_bound_sites = this.currentlyBoundSites()
        all_threeway_migration_transitions = []
        for direction in ['toward5prime', 'toward3prime']:
            for edge in this.current_edges:
                for (s1,s2) in edge.bothWaysRound():
                    if direction == 'toward5prime': # Test for displacement toward 5' end
                        s1pr = s1.threePrimeAdjacentSite()
                        s2pr = s2.fivePrimeAdjacentSite()
                    elif direction == 'toward3prime': # Test for displacement toward 3' end
                        s1pr = s1.fivePrimeAdjacentSite()
                        s2pr = s2.threePrimeAdjacentSite()
                    else:
                        assert False
                    if s1pr is not None and s2pr is not None:
                        s3 = this.getBindingPartner(s2pr)
                        if s3 is not None: # I.e., s2' is already bound to something
                            if s1pr not in currently_bound_sites: ## I.e., s1' is NOT already bound to something
                                edge_to_add = Edge(s1pr, s2pr)
                                if edge_to_add in possible_new_edges:
                                    edge_to_remove = Edge(s2pr, s3)
                                    new_strand_graph = this.removeEdgeFromCurrentEdges(edge_to_remove).addEdgeToCurrentEdges(edge_to_add)
                                    edges_added_in_transition = [edge_to_add]
                                    edges_removed_in_transition = [edge_to_remove]
                                    all_edges_involved_in_transition = sorted(edges_added_in_transition + edges_removed_in_transition)
                                    this_threeway_migration_transition = {'type':'THREE_WAY_MIGRATION',
                                                                          'edges_added':edges_added_in_transition,
                                                                          'edges_removed':edges_removed_in_transition,
                                                                          'all_edges_involved':all_edges_involved_in_transition,
                                                                          'new_strand_graph':new_strand_graph}
                                    all_threeway_migration_transitions.append(this_threeway_migration_transition)
        return all_threeway_migration_transitions

    # This version uses the "anchored" predicate to decide if the invader is close enough to the incumbent for displacement to proceed.
    # The "threeWayMode" setting determines which version of the "anchored" predicate to use, which determines the behavior of this function.
    # If 'anchored_strgsd' then the version of "anchored" from strgsd is used.
    #   The 'strgsd' version permits 3-way initiated 4-way migration, but not remote toehold-style reactions.
    # If 'anchored_traversal' then a more permissive version of "anchored" based on graph traversal is used.
    #   The 'traversal' version permits both 3-way initiated 4-way migration and remote toehold-style reactions.
    def allThreeWayMigrationTransitions_Anchored(self, this):
        possible_new_edges = this.possibleNewEdges()
        currently_unbound_sites = this.currentlyUnboundSites()
        all_threeway_migration_transitions = []
        for edge_to_remove in this.current_edges:
            for (s1, s2) in edge_to_remove.bothWaysRound():
                for s in currently_unbound_sites:
                    edge_to_add = Edge(s, s2)
                    if edge_to_add in possible_new_edges:
                        new_strand_graph = this.removeEdgeFromCurrentEdges(edge_to_remove).addEdgeToCurrentEdges(edge_to_add)
                        if self.anchored(new_strand_graph, edge_to_add):
                            ####################################################################################################################################
                            #                                                                                                                                  #
                            # TO DO: if both ends of the new edge are on the same strand, check that no bonds will become hidden as a result of this reaction! #
                            #        Note that this is only really necessary for anchored_traversal mode, which is disabled until I can add this test...       #
                            #                                                                                                                                  #
                            ####################################################################################################################################
                            #if edge_to_add.withinOneStrand():
                            edges_added_in_transition = [edge_to_add]
                            edges_removed_in_transition = [edge_to_remove]
                            all_edges_involved_in_transition = sorted(edges_added_in_transition + edges_removed_in_transition)
                            this_threeway_migration_transition = {'type':'THREE_WAY_MIGRATION',
                                                                  'edges_added':edges_added_in_transition,
                                                                  'edges_removed':edges_removed_in_transition,
                                                                  'all_edges_involved':all_edges_involved_in_transition,
                                                                  'new_strand_graph':new_strand_graph}
                            all_threeway_migration_transitions.append(this_threeway_migration_transition)
        return all_threeway_migration_transitions

    def allThreeWayMigrationTransitions(self, this):
        if self.settings['threeWayMode'] == 'adjacent':
            return self.allThreeWayMigrationTransitions_Adjacent(this)
        elif self.settings['threeWayMode'] in ['anchored_strgsd']: ## , 'anchored_traversal']:
            return self.allThreeWayMigrationTransitions_Anchored(this)
        else:
            assert False

    # Look for 4-way branch migration transitions from "this" species.
    # This function only supports 4-way branch migration, and not n-way migration for n>4.
    # TO DO: maybe generalize to n-way migration?!
    def allFourWayMigrationTransitions(self, this, debug=False):
        def debugPrint(x):
            if debug:
                print(x)
        possible_new_edges = this.possibleNewEdges()
        currently_bound_sites = this.currentlyBoundSites()
        all_fourway_migration_transitions = []
        for edge in this.current_edges:
            for (s1,s2) in edge.bothWaysRound():
                debugPrint('%%%')
                debugPrint('Testing edge where s1='+str(s1)+' and s2='+str(s2))
                s1pr = s1.threePrimeAdjacentSite()
                s2pr = s2.fivePrimeAdjacentSite()
                debugPrint('s1pr = '+str(s1pr))
                debugPrint('s2pr = '+str(s2pr))
                if s1pr is not None and s2pr is not None:
                    s3 = this.getBindingPartner(s1pr)
                    debugPrint('s3 = '+str(s3))
                    if s3 is not None:
                        s3pr = s3.threePrimeAdjacentSite()
                        debugPrint('s3pr = '+str(s3pr))
                        if s3pr is not None:
                            s4 = this.getBindingPartner(s3pr)
                            debugPrint('s4 = '+str(s4))
                            if s4 is not None:
                                s4pr = s4.threePrimeAdjacentSite()
                                debugPrint('s4pr = '+str(s4pr))
                                if s4pr is not None:
                                    s4pr_binding_partner = this.getBindingPartner(s4pr)
                                    debugPrint('s4pr_binding_partner = '+str(s4pr_binding_partner))
                                    if s4pr_binding_partner is not None and s4pr_binding_partner == s2pr:
                                        first_edge_to_add = Edge(s1pr, s2pr)
                                        second_edge_to_add = Edge(s3, s4pr)
                                        debugPrint('first_edge_to_add = '+str(first_edge_to_add)+'    ...Possible? '+str(first_edge_to_add in possible_new_edges))
                                        debugPrint('second_edge_to_add = '+str(second_edge_to_add)+'    ...Possible? '+str(second_edge_to_add in possible_new_edges))
                                        if first_edge_to_add in possible_new_edges and second_edge_to_add in possible_new_edges:
                                            first_edge_to_remove = Edge(s1pr, s3)
                                            second_edge_to_remove = Edge(s2pr, s4pr)
                                            debugPrint('first_edge_to_remove = '+str(first_edge_to_remove))
                                            debugPrint('second_edge_to_remove = '+str(second_edge_to_remove))
                                            edges_added_in_transition = sorted([first_edge_to_add, second_edge_to_add])
                                            edges_removed_in_transition = sorted([first_edge_to_remove, second_edge_to_remove])
                                            all_edges_involved_in_transition = sorted(edges_added_in_transition + edges_removed_in_transition)
                                            new_strand_graph = this.removeEdgeFromCurrentEdges(first_edge_to_remove) \
                                                                   .removeEdgeFromCurrentEdges(second_edge_to_remove) \
                                                                   .addEdgeToCurrentEdges(first_edge_to_add) \
                                                                   .addEdgeToCurrentEdges(second_edge_to_add)
                                            if all_edges_involved_in_transition not in [d['all_edges_involved'] for d in all_fourway_migration_transitions]:
                                                this_fourway_migration_transition = {'type':'FOUR_WAY_MIGRATION',
                                                                                     'edges_added':edges_added_in_transition,
                                                                                     'edges_removed':edges_removed_in_transition,
                                                                                     'all_edges_involved':all_edges_involved_in_transition,
                                                                                     'new_strand_graph':new_strand_graph}
                                                debugPrint('TRANSITION INFO: '+str(this_fourway_migration_transition))
                                                all_fourway_migration_transitions.append(this_fourway_migration_transition)
        return all_fourway_migration_transitions

    # Get all unimolecular transitions possible from "this" species
    def allUnimolecularTransitions(self, this):
        unbindingTransitions = self.allUnbindingTransitions(this)
        threeWayMigrationTransitions = self.allThreeWayMigrationTransitions(this)
        fourWayMigrationTransitions = self.allFourWayMigrationTransitions(this)
        allTransitions = unbindingTransitions + threeWayMigrationTransitions + fourWayMigrationTransitions
        return allTransitions

    ########################################################################

    # Compute all unimolecular reactions possible starting from "this" species
    def unimolecularReactions(self, this):
        allTransitions = self.allUnimolecularTransitions(this)
        allReactions = []
        reactants = [this]
        for t in allTransitions:
            thisFwdRate = 1.0 # NB: will probably want more realistic rate constants when it comes to actually running simulations.
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
            thisFwdRate = 1.0 # NB: will probably want more realistic rate constants when it comes to actually running simulations.
            theseProducts = [speciesFromStrandGraph(sg) for sg in t['new_strand_graph'].connectedComponents()]
            thisMetadata = {'type':t['type'], 'edges_added':t['edges_added'], 'edges_removed':t['edges_removed'], 'all_edges_involved':t['all_edges_involved']}
            thisReaction = Reaction(reactants, thisFwdRate, theseProducts, bwdrate=None, metadata=thisMetadata)
            if thisReaction not in allReactions:
                allReactions += [thisReaction]
        return allReactions

    # Compute all merged unimolecular reactions starting from "this" species
    def unimolecularMergedReactionProducts(self, this):
        assert self.settings['enumerationMode'] == 'infinite'
        # Debug printing
        def debugPrint(x):
            if self.settings['debug']:
                print(x)
        # A helper function, to find state indexes given their contents.
        # We assume here that these_species is sorted.
        def find_state(state_space_repr, these_species):
            for (pnsidx,pns) in enumerate(state_space_repr):
                if sorted(pns['species_list']) == these_species:
                    return pnsidx
            return None
        #####
        # Step 1: do a "state space" enumeration of the possible unimolecular reactions, starting from this species
        #####
        states_to_process = [0]
        states_processed = []
        state_space_repr = [{'species_list':[this], 'out_edge_destinations':[]}]
        while states_to_process != []:
            # Get next state off the list
            sidx = states_to_process.pop(0)
            next_states = []
            these_species = state_space_repr[sidx]['species_list']
            # For each species, find all possible unimolecular transitions and form new (sorted) species lists for them
            for (xidx,x) in enumerate(state_space_repr[sidx]['species_list']):
                for t in self.allUnimolecularTransitions(x):
                    theseProducts = [speciesFromStrandGraph(sg) for sg in t['new_strand_graph'].connectedComponents()]
                    thisNextState = these_species[0:xidx] + theseProducts + these_species[(xidx+1):]
                    next_states += [sorted(thisNextState)]
            # Try to find each of the (sorted) species lists in the existing state space
            for ns in next_states:
                nsidx = find_state(state_space_repr, ns)
                if nsidx == None:
                    # nsidx == None means that ns is a new state, so we need to:
                    #  - Add it to the end of the state_state_repr list, creating a new index for it
                    #  - Add that index to the list of states that need to be processed
                    state_space_repr += [{'species_list':ns, 'out_edge_destinations':[]}]
                    nsidx = len(state_space_repr)-1
                    states_to_process += [nsidx]
                # Add the next state index (now sure to exist) to the out edge destinations for the current state
                assert nsidx != sidx
                state_space_repr[sidx]['out_edge_destinations'] += [nsidx]
                # If the next state index has not been processed and is not queued up already to be processed:
                #  - Add the next state index to the list of states waiting to be processed
                if (nsidx not in states_processed) and (nsidx not in states_to_process):
                    states_to_process += [nsidx]
            # Add this state index to the list of states already processed
            states_processed += [sidx]
        assert len(states_processed) == len(state_space_repr)
        debugPrint('STATE SPACE REPRESENTATION:')
        debugPrint(state_space_repr)
        #####
        # Step 2: find all trajectories through the state space, starting from state 0 and ending either in state 0 or another terminal state
        #####
        total_num_edges = sum([len(s['out_edge_destinations']) for s in state_space_repr])
        debugPrint('TOTAL NUMBER OF EDGES: '+str(total_num_edges))
        all_trajectories = None
        for i in range(total_num_edges):
            trajectories_modified_in_this_iteration = 0
            debugPrint('All trajectories: '+str(all_trajectories))
            if i == 0:
                ## Find the initial steps, starting from state 0
                assert all_trajectories is None
                first_new_states = state_space_repr[0]['out_edge_destinations']
                assert len(first_new_states) != 0
                all_trajectories = [[0,s] for s in first_new_states]
            else:
                ## Find steps to extend the existing trajectories
                assert all_trajectories is not None
                new_trajectories = []
                for tidx in range(len(all_trajectories)):
                    fsidx = all_trajectories[tidx][-1]
                    debugPrint('Trajectory index '+str(tidx)+' has final state '+str(fsidx))
                    new_states = state_space_repr[fsidx]['out_edge_destinations']
                    debugPrint('New states reachable from state '+str(fsidx)+': '+str(new_states))
                    if new_states == []:
                        new_trajectories += [all_trajectories[tidx]]
                    else:
                        for nsidx in new_states:
                            new_trajectories += [list(all_trajectories[tidx])+[nsidx]]
                            trajectories_modified_in_this_iteration += 1
                debugPrint('New trajectories: '+str(new_trajectories))
                all_trajectories = new_trajectories
                # Stop if we have got to a point where no more steps are possible in any trajectory
                if trajectories_modified_in_this_iteration == 0:
                    debugPrint('No trajectories modified: breaking out of loop')
                    break
        # SPECIAL CASE:
        # If all_trajectories is still None at this point, it means that no reactions at all are possible from this initial state, even unbinding reactions.
        # An example of this is if you have two non-toehold domains binding.
        # In this case, we return a single trajectory consiting just of the initial state. This will be accounted for in step 3 below.
        if all_trajectories is None:
            all_trajectories = [[0]]
        debugPrint('ALL TRAJECTORIES: '+str(all_trajectories))
        debugPrint('Total found: '+str(len(all_trajectories)))
        #####
        # Step 3: Filter the trajectories, removing those that end in a non-stable state or that end in the initial state and contain more than one state.
        #####
        filtered_trajectories = []
        for t in all_trajectories:
            if (t[-1] == 0) and (len(t) > 1):
                # Skip any trajectory that ends in the initial state and contains more than one state (to handle the SPECIAL CASE from above)
                pass
            elif (len(state_space_repr[t[-1]]['out_edge_destinations']) != 0):
                # Skip any trajectory that ends in a state with possible outward edges
                pass
            else:
                # Any other trajectory should be kept
                filtered_trajectories += [t]
        debugPrint('FILTERED TRAJECTORIES: '+str(filtered_trajectories))
        #####
        # Step 4: the end points of those that remain correspond to products from valid merged reactions (there may be duplicates to remove).
        #####
        filtered_endpoints = []
        for ft in filtered_trajectories:
            end_state = ft[-1]
            if end_state not in filtered_endpoints:
                filtered_endpoints += [end_state]
        debugPrint('FILTERED ENDPOINTS: '+str(filtered_endpoints))
        allProducts = [sorted(state_space_repr[ep]['species_list']) for ep in filtered_endpoints]
        debugPrint('ALL PRODUCTS: '+str(allProducts))
        return allProducts    
    
    # Compute all bimolecular reactions possible when "this" species is paired with "that" species, with any subsequent "fast" reactions being merged
    def bimolecularMergedReactions(self, this, that):
        assert self.settings['enumerationMode'] == 'infinite'
        allTransitions = self.allBindingTransitions(this.compose(that))
        allReactions = []
        reactants = [this, that]
        for t in allTransitions:
            thisFwdRate = 1.0 # NB: will probably want more realistic rate constants when it comes to actually running simulations.
            theseInitialProducts = [speciesFromStrandGraph(sg) for sg in t['new_strand_graph'].connectedComponents()]
            assert len(theseInitialProducts) == 1
            initialProduct = theseInitialProducts[0]
            allMergedProducts = self.unimolecularMergedReactionProducts(initialProduct)
            if allMergedProducts == []:
                # This means that no stable merged products were found.
                # An example would be a complex where branch migration is moving back and forth but the strands are all held on by long domains and thus don't unbind.
                # This could also happen if you have a regular toehold binding but display one strand and then end up in a situation with no single stable state.
                # So, this might occur in some cases. TO DO: think more about this!
                lib.error('SOME REACTION PRODUCES NO STABLE MERGED PRODUCT(S). THIS IS NOT YET SUPPORTED BY THE REACTION ENUMERATOR! I SUGGEST THAT YOU USE DETAILED MODE!')
            else:
                for theseProducts in allMergedProducts:
                    if sorted(theseProducts) == sorted(reactants):
                        # This means that the reactants were returned as products.
                        # This reaction would be "unproductive" if these were the only merged products found, so we will just ignore this reaction.
                        pass
                    else:
                        # thisMetadata = {'type':t['type'], 'edges_added':t['edges_added'], 'edges_removed':t['edges_removed'], 'all_edges_involved':t['all_edges_involved']}
                        thisMetadata = {} # NB: may need a more nuanced approach to merging metadata when merging reactions.
                        thisReaction = Reaction(reactants, thisFwdRate, theseProducts, bwdrate=None, metadata=thisMetadata)
                        if thisReaction not in allReactions:
                            allReactions += [thisReaction]
        return allReactions

    ########################################################################

    #
    # FINALLY, IMPLEMENT THE ABSTRACT METHOD FROM THE ABSTRACT SUPERCLASS!!!
    # Given a list of distinct species, enumerate all possible reactions that could occur.
    # This version encompasses both the DETAILED and INFINITE semantics for reaction enumeration.
    # Also commenting out all the debug printing in this version, to make important control flow more obvious.
    #
    def enumerateReactions(self, species_list):
        assert self.validSettings()
        #def debugPrint(x):
        #    if self.settings['debug']:
        #        print(x)
        if not isListOfSpecies(species_list):
            lib.error('In ReactionEnumerator_Original.enumerateReactions: expected list of species as argument, but found: '+str(species_list))
        if not lib.distinct(species_list): 
            lib.error('In ReactionEnumerator_Original.enumerateReactions: expected all species in argument list to be unique, but found: '+str(species_list))
        allReactions = []
        species_processed = []
        species_pairs_processed_SORTED = []
        species_to_process = list(species_list)
        #iterationcount = 1
        while (species_to_process != []):
            #msg = 'Iteration: '+str(iterationcount)
            #debugPrint(msg)
            #debugPrint('-'*len(msg))
            #debugPrint('At start of iteration:')
            #debugPrint('species_processed:')
            #for s in species_processed:
            #    debugPrint(s)
            ##debugPrint('species_pairs_processed_SORTED:')
            ##for sp in species_pairs_processed_SORTED:
            ##    debugPrint(sp)
            #debugPrint('species_to_process:')
            #for s in species_to_process:
            #    debugPrint(s)
            #debugPrint('allReactions:')
            #for r in allReactions:
            #    debugPrint(r)
            x = species_to_process.pop(0) # Remove and return first species in the list
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
                        newReactions += self.bimolecularReactions(x, y)
                    elif self.settings['enumerationMode'] == 'infinite':
                        newReactions += self.bimolecularMergedReactions(x, y)
                    else:
                        assert False
                    species_pairs_processed_SORTED += [this_sorted_pair]
            #debugPrint('newReactions found:')
            #for r in newReactions:
            #    debugPrint(r)
            possiblyNewSpecies = []
            for r in newReactions:
                assert r not in allReactions
                #if r in allReactions:
                #    debugPrint('Reaction already seen?!')
                #    debugPrint('New reaction:')
                #    debugPrint(r)
                #    debugPrint('allReactions as it currently stands:')
                #    debugPrint(len(allReactions))
                #    for q in allReactions:
                #        debugPrint(q)
                #    assert False
                allReactions += [r]
                possiblyNewSpecies += r.listOfSpeciesInvolved()
            species_processed += [x] # Do this before the next loop so we don't double-count species!
            for pns in possiblyNewSpecies:
                if (pns not in species_processed) and (pns not in species_to_process):
                    species_to_process += [pns]
            #debugPrint('At end of iteration:')
            #debugPrint('species_processed:')
            #for s in species_processed:
            #    debugPrint(s)
            ##debugPrint('species_pairs_processed_SORTED:')
            ##for sp in species_pairs_processed_SORTED:
            ##    debugPrint(sp)
            #debugPrint('species_to_process:')
            #for s in species_to_process:
            #    debugPrint(s)
            #debugPrint('allReactions:')
            #for r in allReactions:
            #    debugPrint(r)
            #iterationcount += 1
        return CRN(species_processed, allReactions)

###############################################################################################

# 
# Enumerator objects based on common settings
# 

enumeratorOriginal_Adjacent_Detailed = ReactionEnumerator_Original({'name':'adjacent_detailed', 
                                                                    'debug': False,
                                                                    'enumerationMode':'detailed',
                                                                    'maxComplexSize': math.inf,
                                                                    'threeWayMode':'adjacent',
                                                                    'unbindingMode':'adjacent'})

enumeratorOriginal_Adjacent_Infinite = ReactionEnumerator_Original({'name':'adjacent_infinite',
                                                                    'debug': False,
                                                                    'enumerationMode':'infinite',
                                                                    'maxComplexSize': math.inf,
                                                                    'threeWayMode':'adjacent',
                                                                    'unbindingMode':'adjacent'})

enumeratorOriginal_AnchoredStrgsd_Detailed = ReactionEnumerator_Original({'name':'anchored_strgsd_detailed',
                                                                          'debug': False,
                                                                          'enumerationMode':'detailed',
                                                                          'maxComplexSize': math.inf,
                                                                          'threeWayMode':'anchored_strgsd',
                                                                          'unbindingMode':'anchored_strgsd'})

enumeratorOriginal_AnchoredStrgsd_Infinite = ReactionEnumerator_Original({'name':'anchored_strgsd_infinite',
                                                                          'debug': False,
                                                                          'enumerationMode':'infinite',
                                                                          'maxComplexSize': math.inf,
                                                                          'threeWayMode':'anchored_strgsd',
                                                                          'unbindingMode':'anchored_strgsd'})

# enumeratorOriginal_AnchoredTraversal_Detailed = ReactionEnumerator_Original({'name':'anchored_traversal_detailed',
#                                                                              'debug': False,
#                                                                              'enumerationMode':'detailed',
#                                                                              'maxComplexSize': math.inf,
#                                                                              'threeWayMode':'anchored_traversal',
#                                                                              'unbindingMode':'anchored_traversal'})
# 
# enumeratorOriginal_AnchoredTraversal_Infinite = ReactionEnumerator_Original({'name':'anchored_traversal_infinite',
#                                                                              'debug': False,
#                                                                              'enumerationMode':'infinite',
#                                                                              'maxComplexSize': math.inf,
#                                                                              'threeWayMode':'anchored_traversal',
#                                                                              'unbindingMode':'anchored_traversal'})

allEnumeratorsOriginal = [enumeratorOriginal_Adjacent_Detailed,
                          enumeratorOriginal_AnchoredStrgsd_Detailed,
                          enumeratorOriginal_Adjacent_Infinite,
                          enumeratorOriginal_AnchoredStrgsd_Infinite]
#                          enumeratorOriginal_AnchoredTraversal_Infinite,
#                          enumeratorOriginal_AnchoredTraversal_Detailed]

enumeratorOriginal_Default = enumeratorOriginal_AnchoredStrgsd_Detailed

###############################################################################################
