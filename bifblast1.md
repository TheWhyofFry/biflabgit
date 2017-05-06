# BifBlast! Introduction

This documentation is a primer for BifBlast! Herein I will list the goals and organization of the code/modules we will develop as a group.  


# Goals
	* Develop a core backend that can:
		* Create a BLAST database from a list of sequences
		* Perform arbitrary blasts against this database
		* Parse the results
		* Keep a reference of the created BLAST databases in a SQLite database
		* Allow the module to be executed as a stand-alone application
			* e.g. ``python -m bifblast --createdb MYBLASTDB --sequences MYSEQUENCES --dbtype nt``
	* Transform this into a neat Python package
	* Continuously update the Git repository


## Core development

## Small introduction (recap?) of Python modules

Python modules allow the encapsulation of a set of functions and classes for re-use. It would be very redundant to copy all the code needed for the ``math`` module into your code to access the functions.  In the simplest case, you can have a single python file with a bunch of functions, ``myfunctions.py``.  In order for another script to use it, it has to import this file with ``import myfunctions``.  Python has a set of paths where it looks for functions and gives an error if it didn't find the functions in the predefined path.  You could, of course, add paths for python to search, but this gets very messy very quickly.  This is why when you install python modules, they are put in a central location.  This location may be accessible system-wide or in a local folder. Usually, python stores (third-party) installed modules in a ``site-packages`` folder, so when you import something it will look inside this folder for the appropriate module. So if you enter ``import Bio`` it will import the ``site-packages/Bio`` module that exists as a subfolder.  Modules may also contain ``submodules`` which can be imported using ``.`` as a separator, for instance ``import Bio.Seq``.  To avoid using the ``Bio.`` prefix, one can also use ``from Bio import Seq`` and then access the ``Seq`` modules contents with ``Seq.``. 




### Module overview

For each task, a separate Python script will be created.  This means that BLAST database generation functions will not be group together with BLAST query functions.  The reason for this is to make our library as modular as possible.  This allows for a smoother history and if bugs do creep up, we can easily isolate them to a particular module. 

The Biopython package already provides most of the functionality for our project, so most of the core functions will just act as shorthand versions.  This step is a bit redundant, but I would like, especially the Hons students, to get familiar with general Python concepts too. 


Lets look at the general structure of our core package and discuss each file in turn. 

```bash
.
├── bifblast
│   ├── blastargparser.py
│   ├── blastdbgen.py
│   ├── blastdbsqlite.py
│   ├── blastparser.py
│   ├── blastquery.py
│   ├── __init__.py
│   ├── __main__.py
│   └── setup.py
├── MAINTAINERS.md
└── README.md

```



#### blastdbgen.py

This file will handle the blastdb generation queries.  It will accept a file as input with options stating whether it is a protein or nucleotide BLAST database.  While this seems trivial, we need to take a few things into account:
	
	* The input file must be in _valid_ FASTA format
	* The headers of the FASTA entries may not include escape characters (such as tab)
	* An option should exist to automatically determine if the input sequences are protein or nucleotide sequences (within reason)

As far as I know (waiting for the egg on my face) there is no direct way to create a blast database from wihthin BioPython. To this end, we will need to create a function that interfaces with the NCBI BLAST+ suite to create a blast database for us. From a usage perspective, the BLAST database generation should happen via a single function call.

#### blastdbquery.py

This module will facilitate the execution of BLAST queries.  Essentially, a function will be written that can execute the query (via BioPython) but also handle where the raw results (before parsing) are to be stored. There should be an option for the user to get a raw string output or output the results to a file.  While there can be an option choosing what format the output should be, the default needs to be XML, since it can convert to any other format (e.g. tab-delimited).  As with ``blastdbgen``, sequences need to be validated as valid FASTA. 

#### blastparser.py

This module will parse raw BLAST results.  It will contain functions that can both parse and selectively extract info from the BLAST hits. Deviating slightly from the pure modular design, it should also be able to call the necessary BioPython functions that allow conversion between the BLAST output formats. 

The ``pandas`` Python module is great for working with (among other things) tabular data.  BLAST results can be converted to tabular format and if structured appropriately, it will make filtering BLAST results a lot easier, e.g. filtering by E-value, bit score, query hit size etc.

Ideally, this module should also include a summary statistic that includes:
	* Number of hits found per BLAST query
	* Total coverage of the BLAST query 

#### blastdbsqlite.py

Not absolutely necessary for the core module to succeed, but given that we eventually want a web interface for this project, it seems prudent to keep a structured record of all the BLAST databases that have been created.  This allows the user to access previously created BLAST databases using a key and a reference BLAST database.  This module will not be assigned to anyone per se, but if you would like to work on this, let me know.  I'd suggest you read up on the [SqlAlchemy](https://www.sqlalchemy.org) module.

#### ``__init__.py`` and ``__main__.py``

These two files are necessary for packaging our module.  The ``__init__.py`` both indicates that the ``bifblast`` directory should be treated as a module and the files therein can be imported as ``from bifblast import blastdbgen``, ``import bifblast.blastdbgen`` etc.  It also allows some module initialization (and other packaging functions). The ``__main__.py`` file helps with the execution of the module as a command.  Subsequent lectures will delve into more detail.

#### ``setup.py``

This file will contain the necessary functions to setup/install the package. 




## Small examples to get started

To help you get started, we'll look at some of the basic interfaces to BioPython as well as creating a BLAST database from the command line.  The code will generally follow the BioPython [cookbook](http://biopython.org/DIST/docs/tutorial/Tutorial.html).  These examples are merely to help you get familiar with everything.  While your code may look similar in the end, please try not to copy-and-paste it.  



### Creating a BLAST database

Creating a BLAST database is straightforward.  You will need the [NCBI BLAST+ tools](). Download the relevant file for your operating system and install.  These tools are also essential for BioPython's local BLAST functionality.



### Performing arbitrary BLASTs

### Parsing the results

### Keep a reference of the created BLAST databases in a SQLite database

### Allow the module to be executed as a stand-alone application