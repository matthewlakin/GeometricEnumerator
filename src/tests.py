########################################################################
#
# tests.py
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
import sgparser
import sys
from constants import *
from timeit import default_timer as timer
from domain import *
from strand import *
from process import *
from strandgraph import *
from species import *
from reaction import *
from enumerator_geometric import *
from crn import *
from constraintchecker_sampling import *


class Skipping(Exception):
    pass

def skip():
    print('Skipping...')
    raise Skipping

def zero_nucteotide_loop_check(s):
    s1 = [(i.strip()) for i in (s[1:-1] + "").split('|')]
    res = ""
    for strand in s1:
        s4 = (strand[1:-1] + "").split()
        for k in range(0, len(s4) - 1):
            s5 = s4[k].split('!')
            s6 = s4[k+1].split('!')
            if((s5[0] + "*" == s6[0]) or (s5[0] == s6[0] + "*")):
                res += s4[k] + ' '+ s4[k+1] + ' '    
    if(len(res) > 0):
        return (True, "Zero nucteotide loop poosible domain: " + res)
    return (False, res)

def doParsingAndEnumerationTest(s, number, domainLengthStr, maxComplexSize=math.inf):
    flag_zero_nucleotide_loop_check, res = zero_nucteotide_loop_check(s)
    print(res)
    if(flag_zero_nucleotide_loop_check):
        return

    for enumerator in allEnumeratorsGeometric:
        start_time = timer()
        enumerator.settings['maxComplexSize'] = maxComplexSize
        title = 'DOING ENUMERATION TEST '+str(number)+' FOR SETTINGS '+enumerator.settings['name'].upper()+':'
        print(title)
        print('-'*len(title))
        print('Input string:')
        print(str(s))
        print(domainLengthStr)
        print("Nicked Parameters: ")
        #print("NICKED_FLAG: " + str(NICKED_FLAG))
        if(NICKED_FLAG): 
            #print("NICKED ANGLE LOWER BOUND: " + str(NICKEDANGLE_LOWER_BOUND))
            print("NICKED ANGLE UPPER BOUND: " + str(NICKEDANGLE_UPPER_BOUND))
        
        p = sgparser.parse(s)
        sg = strandGraphFromProcess(p)
        sg.displayRepresentation()
        sg.domainLength = parseDomainLength(domainLengthStr)
        speciesList = speciesListFromProcess(p)

        for spec in speciesList:
            spec.domainLength = sg.domainLength

        crn = enumeratorGeometric.enumerateReactions(speciesList)
        print('Found '+str(len(crn.species))+' species and '+str(len(crn.reactions))+' reactions in total.')
        crn.displayRepresentation()
        print()
        end_time = timer()
        elapsed_time = end_time - start_time
        print('Time taken to enumerate reactions for settings '+enumerator.settings['name']+': '+str(elapsed_time)+' seconds')
        print()

def test_branch_migration_leak():
    s = '(<x!j y x*!j> | <x>)'
    domainLengthStr = 'toeholdDomain x length 8 longDomain y length 20'
    doParsingAndEnumerationTest(s, 1, domainLengthStr)

def test003():
    s = '( <L T2^!i2 X*!i1 T1^> | <A X!i1 T2^*!i2> | <T1^* X*!j1 R> | < X!j1 A!j2 > | <A*!j2 > )'
    
    #s = '(< X!i2 A!i1 > | <A*!i1 > | <T1^*!i3 X*!i2 R> | <L T2^!i4 X*!i5 T1^!i3> | <A X!i5 T2^*!i4>)'  
    domainLengthStr = 'longDomain L length 20 toeholdDomain T2 length 5 longDomain X length 20 toeholdDomain T1 length 5 longDomain A length 20 longDomain R length 20'
    doParsingAndEnumerationTest(s, 1, domainLengthStr)

def test004():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20'
    s = '(<B!i1> | <B*!i1 A*!i2> | <A!i2 B> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)


# def test004a():
#     domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain Z length 20'
#     #doProcessTest(p, domainLengthStr)
#     s = '(<A!i1 Z> | <B!i2> | <B*!i2 A*!i1>)'

def test005():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20'
    s = '( <A!i1> | <A B!i2> | <B*!i2 A*!i1> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test005a():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain Z length 20'
    s = '( <A!i1> | <Z B!i2> | <B*!i2 A*!i1> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test006():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain C length 20 longDomain D length 20'
    s = '( <C!i1 A!i2 B> | <B!i3 D!i4> | <D*!i4 B*!i3 A*!i2 C*!i1> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)


def test006a():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain C length 20 longDomain D length 20 longDomain Z length 20'
    s = '( <C!i1 A!i2 Z> | <B!i3 D!i4> | <D*!i4 B*!i3 A*!i2 C*!i1> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test007():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain C length 20 longDomain D length 20'
    s = '( <C!i1 A!i2> | <A B!i3 D!i4> | <D*!i4 B*!i3 A*!i2 C*!i1> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)
                         
def test007a():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain C length 20 longDomain D length 20 longDomain Z length 20'
    s = '( <C!i1 A!i2> | <Z B!i3 D!i4> | <D*!i4 B*!i3 A*!i2 C*!i1> )'
    doParsingAndEnumerationTest(s, 13, domainLengthStr)
                         
def test008():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain C length 20'
    s = '( <A!i1 B!i2> | <B*!i2 C*!i3> | <C!i3 B!i4> | <B*!i4 A*!i1> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test008a():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain C length 20 longDomain Z length 20'
    s = '( <A!i1 B!i2> | <B*!i2 C*!i3> | <C!i3 Z!i4> | <Z*!i4 A*!i1> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)


def test009():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain C length 20 longDomain E length 20 longDomain F length 20'
    s = '( <E!i5 A!i1 B!i2> | <B*!i2 C*!i3 F*!i6> | <F!i6 C!i3 B!i4> | <B*!i4 A*!i1 E*!i5> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test009a():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain C length 20 longDomain E length 20 longDomain F length 20 longDomain Z length 20 '
    s = '( <E!i5 A!i1 B!i2> | <B*!i2 C*!i3 F*!i6> | <F!i6 C!i3 Z!i4> | <Z*!i4 A*!i1 E*!i5> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test010():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain C length 20 longDomain E length 20 longDomain F length 20 longDomain G length 20 longDomain H length 20'
    s = '( <E!i5 A!i1 B!i2 G!i7> | <G*!i7 B*!i2 C*!i3 F*!i6> | <F!i6 C!i3 B!i4 H!i8> | <H*!i8 B*!i4 A*!i1 E*!i5> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test010a():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain C length 20 longDomain E length 20 longDomain F length 20 longDomain G length 20 longDomain H length 20 longDomain Z length 20'
    s = '( <E!i5 A!i1 B!i2 G!i7> | <G*!i7 B*!i2 C*!i3 F*!i6> | <F!i6 C!i3 Z!i4 H!i8> | <H*!i8 Z*!i4 A*!i1 E*!i5> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test011():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20'
    s = '( <A!i1 B!i2> | <B*!i2 A*!i3> | <A!i3 B!i4> | <B*!i4 A*!i1> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test012():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20'
    s = '( <A^!i1 B^!i2> | <B^*!i2 A^*!i3> | <A^!i3 B^!i4> | <B^*!i4 A^*!i1> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test013():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain E length 20 longDomain F length 20 longDomain G length 20 longDomain H length 20'
    s = '( <E!i5 A!i1 B!i2 F!i7> | <F*!i7 B*!i2 A*!i3 G*!i6> | <G!i6 A!i3 B!i4 H!i8> | <H*!i8 B*!i4 A*!i1 E*!i5> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test014():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain E length 20 longDomain F length 20 longDomain G length 20 longDomain H length 20 longDomain I length 20 longDomain J length 20 longDomain K length 20 longDomain L length 20'
    s = '( <E A!i1 B!i2 F> | <G B*!i2 A*!i3 H> | <I A!i3 B!i4 J> | <K B*!i4 A*!i1 L> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test015():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain E length 20 longDomain F length 20 longDomain G length 20 longDomain H length 20 longDomain I length 20 longDomain J length 20 longDomain K length 20 longDomain L length 20'
    s = '( <E A!i1 B!i2 F> | <G B*!i2 A*!i3 H> | <I A!i3 B!i4 J> | <K B*!i4 A*!i1 L> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test016():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20'
    s = '( <A!i1 B> | <B!i2> | <B*!i2 A*!i1> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test017():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain C length 20'
    s = '( <A!i1 B!i2> | <B*!i2 C*!i3> | <C!i3 B!i4> | <B*!i4 A*!i1> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)
    
def test018():
    domainLengthStr = 'longDomain A length 20 longDomain B length 20'
    s = '( <A!i1 B!i2> | <B*!i2 A*!i3> | <A!i3 B!i4> | <B*!i4 A*!i1> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test019():

    domainLengthStr = 'longDomain A length 20 longDomain B length 20'
    s = '( <A B!i1> | <A!i2 B> | <A*!i2 B*!i1> )'
    p = sgparser.parse(s)
    sg = strandGraphFromProcess(p)
    sg.domainLength = parseDomainLength(domainLengthStr)
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test034():
    domainLengthStr = 'longDomain x length 20'
    s = '( <x> | <x*> )'
    doParsingAndEnumerationTest(s, 1, domainLengthStr)

def test035():
    domainLengthStr = 'toeholdDomain t length 5'
    s = '( <t^> | <t^*> )'
    doParsingAndEnumerationTest(s, 2, domainLengthStr)

def test036():
    domainLengthStr = 'toeholdDomain t length 5 longDomain x length 20'
    s = '( <t^ x> | <x!i1> | <x*!i1 t^*> )'
    doParsingAndEnumerationTest(s, 3, domainLengthStr)

def test037():
    domainLengthStr = 'toeholdDomain t length 5 longDomain x length 20 toeholdDomain u length 5'
    s = '(<t^ x> | <x!i1 u^!i2> | <u^*!i2 x*!i1 t^*>)'
    doParsingAndEnumerationTest(s, 4, domainLengthStr)
    
# def test038(): # From Groves 2015
#     try:
#         domainLengthStr = 'longDomain a length 20 toeholdDomain 2 length 5 toeholdDomain 1 length 5'
#         s = '( <a*!i1 2^> | <1^ a!i1> | <2^* a!i2> | <a*!i2 1^*> )'
#         doParsingAndEnumerationTest(s, 5, domainLengthStr, maxComplexSize=10)
#     except SystemExit:
#         print('*** CAUGHT SYS.EXIT!!! CONTINUING WITH THE TESTS...')

def test039(): # From logplc
    domainLengthStr = 'longDomain b length 20 toeholdDomain tb length 5 toeholdDomain tx length 5 longDomain x length 20 toeholdDomain to length 5'
    s = '''( <tb^ b> | <tx^ x> | <to^*!1 x*!2 tx^*!3 b*!4 tb^*>
           | <b!4 tx^!3> | <x!2 to^!1> )'''
    doParsingAndEnumerationTest(s, 6, domainLengthStr)

def test040(): # From Chatterjee 2017 (NO TETHERS THOUGH!)
    #skip()
    domainLengthStr = 'longDomain s length 20 toeholdDomain a0 length 8 toeholdDomain f length 8 toeholdDomain x length 8 longDomain y length 20'
    s = '( <s a0^> | <a0^* s*!1 f^ s!1> | <s!2 x^ s*!2 f^*> | <x^* s*!3 y^ s!3> | <s*!4 y^*> | <s!4> )'
    doParsingAndEnumerationTest(s, 7, domainLengthStr)

# def test040i(): # From Chatterjee 2017 (NO TETHERS THOUGH!)
#     #skip()
#     domainLengthStr = 'longDomain s length 20 toeholdDomain a0 length 8 toeholdDomain f length 8 toeholdDomain x length 8 longDomain y length 20'
#     s = '(<a0^* s*!1 x f^ y s!1> | <s!2 x^ s*!2 f^*> )'
#     doParsingAndEnumerationTest(s, 7, domainLengthStr)

def test_examples_spacer1():
    domainLengthStr = 'longDomain x length 20 toeholdDomain spcr1 length 6 toeholdDomain spcr2 length 6 longDomain y length 20'
    s = '(<x!i1 spcr1^ y* spcr2^ x*!i1> | <y> )'
    p = sgparser.parse(s)
    sg = strandGraphFromProcess(p)
    sg.domainLength = parseDomainLength(domainLengthStr)
    sg.displayRepresentation() 
    cc = allEnumeratorsGeometric[0].settings['constraintChecker']
    for spec in sg.connectedComponents():
        spec.domainLength = sg.domainLength
        flag, no = cc.isPlausible(spec)
        print(flag)
        print(no)
    doParsingAndEnumerationTest(s, 4, domainLengthStr)

def test_simple_loop_example():
        domainLengthStr = 'toeholdDomain s length 7 longDomain z length 20 longDomain y length 5 toeholdDomain x length 5'
        s = '(<x^!i1 s^ z y*  x^*!i1> | <y>)'
        doParsingAndEnumerationTest(s, 4, domainLengthStr)

def test40_part():
        domainLengthStr = 'toeholdDomain a length 5 longDomain s length 20 toeholdDomain f length 5 toeholdDomain x length 5'
        s = '( <a^*!i1 s*!i2 f^!i3 s!i4> | <s a^!i1> | <s!i2 x^ s*!i4 f^*!i3> )'
        #doParsingAndEnumerationTest(s, 4, domainLengthStr)
        p = sgparser.parse(s)
        sg = strandGraphFromProcess(p)
        sg.domainLength = parseDomainLength(domainLengthStr)
        sg.displayRepresentation()
        rg = regionGraphFromStrandGraph(sg)
        rg.displayRepresentation()
        cc = enumeratorGeometric.settings['constraintChecker']
        print(cc.isPlausible(sg))

# def test060():
#         s = '(<x!i1 y1!i2 z y2!i3 x*!i1> | <y1*!i2> | <y2*!i3>)'
#         domainLengthStr = 'longDomain x length 20 longDomain y1 length 20 longDomain z length 20 longDomain y2 length 20'
#         print(s)
#         print(domainLengthStr)
#         p = sgparser.parse(s)
#         sg = strandGraphFromProcess(p)
#         sg.domainLength = parseDomainLength(domainLengthStr)
#         sg.displayRepresentation()
#         rg = regionGraphFromStrandGraph(sg)
#         rg.displayRepresentation()
#         cc = enumeratorGeometric.settings['constraintChecker']
#         print("is Plausible: " + str(cc.isPlausible(sg)))

# def test061():
#         s = '(<A!i1 B!i2> | <A*!i1 B*!i2>)'
#         domainLengthStr = 'longDomain A length 20 longDomain B length 20'
#         print(s)
#         print(domainLengthStr)
#         p = sgparser.parse(s)
#         sg = strandGraphFromProcess(p)
#         sg.domainLength = parseDomainLength(domainLengthStr)
#         sg.displayRepresentation()
#         rg = regionGraphFromStrandGraph(sg)
#         rg.displayRepresentation()
#         cc = enumeratorGeometric.settings['constraintChecker']
#         print("is Plausible: " + str(cc.isPlausible(sg)))



def test_unbinding():
    s = '(<t^!i x!j> | <x u^!k y u^*!k x*!j t^*!i>)'
    domainLengthStr = 'toeholdDomain t length 8 longDomain x length 20 longDomain y length 20 toeholdDomain u length 8'
    doParsingAndEnumerationTest(s, 4, domainLengthStr)

def test_remote_toehold_TMSD():
    s = '(<t^ spcr1 y>  | <y*!i1 spcr2 t^*> | <y!i1>)'
    domainLengthStr = 'toeholdDomain t length 14 longDomain spcr1 length 5 longDomain y length 20 longDomain spcr2 length 5 longDomain x length 15'
    doParsingAndEnumerationTest(s, 4, domainLengthStr)

def test_remote_toehold_ds_ds():
    s = '(<t^ spcr1!i2 y>  | <spcr1*!i2> | <y*!i1 spcr2!i3 t^*> | <y!i1> | <spcr2*!i3>)'
    domainLengthStr = 'toeholdDomain t length 8 longDomain spcr1 length 5 longDomain y length 20 longDomain spcr2 length 5'
    doParsingAndEnumerationTest(s, 4, domainLengthStr)

def test_remote_toehold_ds_ss():
    print("double stranded spacer is on the gate")
    s = '(<t^ spcr1 y>  | <y*!i1 spcr2!i3 t^*> | <y!i1> | <spcr2*!i3>)'
    domainLengthStr = 'toeholdDomain t length 14 longDomain spcr1 length 5 longDomain y length 20 longDomain spcr2 length 5'
    doParsingAndEnumerationTest(s, 4, domainLengthStr)

def test_remote_toehold_ss_ds():
    print("double stranded spacer is on the invading strand")
    s = '(<t^ spcr1!i2 y>  | <spcr1*!i2> | <y*!i1 spcr2 t^*> | <y!i1> )'
    domainLengthStr = 'toeholdDomain t length 14 longDomain spcr1 length 1 longDomain y length 20 longDomain spcr2 length 2 longDomain x length 15'
    doParsingAndEnumerationTest(s, 4, domainLengthStr)        

def test_remote_toehold_ss_ss():
    s = '(<t^ spcr1 y> | <y*!i1 spcr2 t^*> | <y!i1> )'
    domainLengthStr = 'toeholdDomain t length 14 longDomain spcr1 length 1 longDomain y length 20 longDomain spcr2 length 1 longDomain x length 15'
    doParsingAndEnumerationTest(s, 4, domainLengthStr)        

def test_double_bond():
    print("Nicked Parameters: ")
    print("NICKED_FLAG: " + str(NICKED_FLAG)) 

    s = '(<x!i1 y!i2 t!i3> | <y*!i2> | <t*!i3 z!i4 x*!i1> | <z*!i4> )'
    s = '(<x!i1 y!i2> | <y*!i2 x*!i1> )'
    domainLengthStr = 'longDomain x length 20 longDomain y length 20 longDomain z length 20 longDomain t length 20'
    print(s)
    p = sgparser.parse(s)
    sg = strandGraphFromProcess(p)
    sg.domainLength = parseDomainLength(domainLengthStr)
    sg.displayRepresentation()
    rg = regionGraphFromStrandGraph(sg)
    rg.displayRepresentation()
    cc = enumeratorGeometric.settings['constraintChecker']
    print("is Plausible: " + str(cc.isPlausible(sg)))

def test_double_bond_juction():
    print("Nicked Parameters: ")
    print("NICKED_FLAG: " + str(NICKED_FLAG))

    s = '(<x!i1 y!i2> | <y*!i2> | <z!i4 x*!i1> | <z*!i4> )'
    domainLengthStr = 'longDomain x length 20 longDomain y length 20 longDomain z length 20 longDomain t length 20'
    print(s)
    p = sgparser.parse(s)
    sg = strandGraphFromProcess(p)
    sg.domainLength = parseDomainLength(domainLengthStr)
    sg.displayRepresentation()
    rg = regionGraphFromStrandGraph(sg)
    rg.displayRepresentation()
    cc = enumeratorGeometric.settings['constraintChecker']
    print("is Plausible: " + str(cc.isPlausible(sg)))
    print("##############################################################")
    print()

    s = '(<x!i1 y!i2 a!i3> | <y*!i2> | <a*!i3 z!i4 x*!i1> | <z*!i4> )'
    domainLengthStr = 'longDomain x length 20 longDomain y length 20 longDomain z length 20 longDomain a length 20'
    print(s)
    p = sgparser.parse(s)
    sg = strandGraphFromProcess(p)
    sg.domainLength = parseDomainLength(domainLengthStr)
    sg.displayRepresentation()
    rg = regionGraphFromStrandGraph(sg)
    rg.displayRepresentation()
    cc = enumeratorGeometric.settings['constraintChecker']
    print("is Plausible: " + str(cc.isPlausible(sg)))

    print("##############################################################")
    print()

    s = '(<x!i1 spcr1!i2 y!i3> | <y*!i3 spcr1*!i2> | <y z!i4 x*!i1> | <z*!i4> )'
    domainLengthStr = 'longDomain x length 20 longDomain y length 20 longDomain z length 20 longDomain spcr1 length 20'
    print(s)
    print(domainLengthStr)
    p = sgparser.parse(s)
    sg = strandGraphFromProcess(p)
    sg.domainLength = parseDomainLength(domainLengthStr)
    sg.displayRepresentation()
    rg = regionGraphFromStrandGraph(sg)
    rg.displayRepresentation()
    cc = enumeratorGeometric.settings['constraintChecker']
    print("is Plausible: " + str(cc.isPlausible(sg)))

def test_unimolecular_reactions():
    s = '(<t^ spcr1 y spcr2 t^* > )'
    domainLengthStr = 'toeholdDomain t length 14 longDomain spcr1 length 1 longDomain y length 20 longDomain spcr2 length 1 longDomain x length 15'
    doParsingAndEnumerationTest(s, 4, domainLengthStr)        

def test_zero_nucleotide_loops():
    s = '(<t^ t^*>)'
    domainLengthStr = 'toeholdDomain t length 14 longDomain spcr1 length 10 longDomain y length 20'
    doParsingAndEnumerationTest(s, 4, domainLengthStr)        

def getTestNames():
    all_test_names = sorted([fname for fname in globals().keys() if fname.startswith('test')])
    test_names = []
    for arg in sys.argv[1:]:
        if arg in all_test_names:
            test_names += [arg]
    if test_names == []:
        test_names = all_test_names
    return test_names


def runAllTests(subset=None):
    testsToRun = getTestNames() if subset is None else subset
    for test_name in testsToRun:
        print('')
        print('')
        print('')
        print('******************************** '+test_name+' *************************************')
        try:
            globals()[test_name]()
        except Skipping:
            pass


enumeratorGeometric = ReactionEnumerator_Geometric({'name':'adjacent_detailed',
                                                                    'debug': False,
                                                                    'enumerationMode':'detailed',
                                                                    'maxComplexSize': math.inf,
                                                                    'threeWayMode':'adjacent',
                                                                    'unbindingMode':'adjacent',
                                                                    'rate' : {'bind': 0.003, 'unbind' : 0.1, 'migrate' : 1.0, 'displace' : 1.0},
                                                                    'constraintChecker': ConstraintChecker_Sampling(seed=7) })

allEnumeratorsGeometric = [enumeratorGeometric]

runAllTests()
