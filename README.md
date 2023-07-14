# Geometric Enumerator

Code for the paper "A geometric framework for reaction enumeration in computational nucleic acid devices" 

## Description

- This package implements the geometric reaction enumeration framework for domain level DNA strand displacement reactions. 
- The code demonstrates that the geometric framework described in that paper can be used to enumerate reactions. 
- A structure sampling approach was used for constraint solving. 
- This software supports the following domain level reactions: bind, unbind, 3-way branch migration and 4-way branch migration reactions. 
- For further details please refer to the paper.

## Getting Started
- The package contains this README file and the src folder. 
- The src folder contains the code for enumerating the reactions and constraint solving using structure sampling approach.

### Dependencies

* Python 3.x
* ply
* probio\_lib
* Jupyter Notebook (optional for displaying the graphs)
* Graphviz library (optional for displaying the graphs)

### Executing program

- tests.py: defines a series of tests to test the reaction enumeration algorithms. Examples like hairpin binding,
  remote toehold reactions etc can be found here. It can be executed by following command: python3 tests.py
  This file depends directly on the following classes:

    - enumerator\_geometric.py: This class implements the geometric framework for reaction enumeration algorithms. 
      Geometric predicates plausible and same species are defined here besides other rulesets that produces the CRN objects.

    - constraintchecker\_sampling.py : This class implements the constraint solving using the structure sampling approach. 
      It also records statistical information like the number of sampling trials needed to find the plausibility of structures.

    - reaction.py: a class representing reactions between species.
      Includes a notion of "reversible" reactions, which simplifies pretty-printing and visualization.

    - strandgraph.py: a data structure class that represents strand graphs.
  	  
    - species.py: a class for individual species, which is a subclass of StrandGraph. Species must be fully connected.

    - crn.py: a class encapsulating CRNs, containing of species and their reactions.
  

- tests.ipynb: this Jupyter notebook wraps tests.py, which allows the same set of tests to be run in the graphical Jupyter environment.
  Some of the classes define graphical representations that can be used for visual debugging and development.
  The graph visualization relies on the "graphviz" library and associated command-line tool being installed.

### Input/Output format

- The input DNA strands are provided using the process calculus syntax. The lengths of the domains are provided separately. An example: 
```
# Domain length information:
domainLengthStr = 'longDomain A length 20 longDomain B length 20 longDomain C length 20 longDomain D length 20 longDomain Z length 20'

# Details of DNA strands through process calculus syntax
s = '( <C!i1 A!i2> | <Z B!i3 D!i4> | <D*!i4 B*!i3 A*!i2 C*!i1> )'
````

- Output can be produced in the following two formats:
    - Command line:
        - Lists all the possible reactions that can occur and
        - The strand graph information details of the species that are part of the reaction.
- Jupyter notebook: 
    - Lists all the possible reactions and 
    - Displays the graphical version of species and can be saved as Jupyter notebooks in a html or pdf format.


## Authors

Sarika Kumar
Matthew Lakin

## Acknowledgments

This material is based upon work supported by the National Science Foundation under Grants 1518861 and 1814906.
