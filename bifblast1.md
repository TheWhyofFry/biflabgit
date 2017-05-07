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
<p>Hello</p>



## Small examples to get started

To help you get started, we'll look at some of the basic interfaces to BioPython as well as creating a BLAST database from the command line.  The code will generally follow the BioPython [cookbook](http://biopython.org/DIST/docs/tutorial/Tutorial.html).  These examples are merely to help you get familiar with everything.  While your code may look similar in the end, please try not to copy-and-paste it.  



### Creating a BLAST database

Creating a BLAST database is straightforward.  You will need the [NCBI BLAST+ tools](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). Download the relevant file for your operating system and install (Windows: ``.exe``, Lab PCs: ``.x86_64.rpm``).  These tools are also essential for BioPython's local BLAST functionality.  You would also need a [test](./data/test.fsa) file for this exercise, although you can use your own FASTA file. 

The ``makeblastdb`` command is used to create BLAST databases.  You need to specify at least an input FASTA file with the ``-in`` flag (test.fsa), the type of database (nucleotide/protein) with the ``-dbtype`` flag. The ``-out`` flag is used to name the output database, which will act as a reference when using it in a BLAST search.


* Creating a database
```bash
$ > makeblastdb -in test.fsa -dbtype nucl -out test         
... 
$ > blastn -query test.fsa -db test # Just to check. 

$ > ls -l test*

-rw-r--r-- 1 werner.local users 2199 May  4 08:13 test.fsa
-rw-r--r-- 1 werner.local users  276 May  4 08:15 test.nhr
-rw-r--r-- 1 werner.local users  104 May  4 08:15 test.nin
-rw-r--r-- 1 werner.local users  524 May  4 08:15 test.nsq

```

You'll see a bunch of files that were created by the ``makeblastdb`` command. When referencing the databse, you need to just use the name ``test``.  If you don't include a name, the ``test`` file will include the ``.fsa`` suffix (so test.nin would be test.fsa.nin) and you'll need to do a ``blastn`` query with ``-db test.fsa``.


### Installing BioPython (and packages in general)

OK, let's backtrack a bit. If you are unsure how to install Python packages, continue to read this section.  If you are comfortable installing packages, you may skip to the next section. Python packages, as stated, are a collection of files that provide functionality beyond the basic Python install.  In our case, we want to install BioPython. The ``pip`` command is used to install Python packages.  Usually, ``pip install biopython`` should suffice. However, many operating systems nowadays come with different versions of Python.  We will be using Python 2.7.  If you installed BioPython and it still appears missing, consider using ``pip2.7 install biopython``. There is another reason why it can b

If you're using the lab PCs, you may need some extra configurations in order to get an internet connection.  I can't put the details on a public guide, so feel free to ask me :-) 

Note that using ``pip install`` will install a package *system wide*.  This means that it will be accessible to all users that can log in to your PC.  In order for it to work, you will need Adminsitrative privileges.  So on your lab PCs, you would need to execute ``sudo pip install biopython``.  If you would like to not install it system wide, you may use ``sudo pip install --user biopython``.  This will install the packages in the home folder. For the time being, I would suggest installing it system wide instead of as a local user.

Generally, packages that have dependencies, i.e. other packages that need to be installed for them to work, grab and install them too.  BioPython has a little exception, which we will see after attempting to install BioPython.




```bash

$ > pip2.7 install biopython
Collecting biopython
  Downloading biopython-1.69.tar.gz (15.4MB)
    100% |████████████████████████████████| 15.4MB 43kB/s 
Building wheels for collected packages: biopython
....


$ > python
>>> import Bio.Blast
```

If all went well, there shouldn't be any errors.  Ahem, except for a moaning by BioPython about "NumPy" not being installed and that it recommends installing NumPy prior to installing BioPython.  NumPy is a large python package that contains a multitude of methods and classes that facilitate the processing of numerical data in python. So let's install ``NumPy`` and then reinstall ``BioPython``. Please note that if you _*did not get an error about numpy*_, that's OK.

```shell

$ > sudo pip install numpy #Install NumPy
Collecting numpy
  Downloading numpy-1.12.1-cp27-cp27mu-manylinux1_x86_64.whl (16.5MB)
    100% |████████████████████████████████| 16.5MB 88kB/s 
Installing collected packages: numpy

$ > sudo pip uninstall biopython
....
Proceed (y/n)? y 
Successfully uninstalled biopython-1.69

$ > pip install biopython
Collecting biopython
  Using cached biopython-1.69.tar.gz
Building wheels for collected packages: biopython
  Running setup.py bdist_wheel for biopython ... done

```


Then, install BioPython as you did before.  There should (hopefully) be no errors :-) If you do get some errors trying out the commands, check if you prefixed the ``pip`` command with ``sudo``. Alternatively, use the ``pip install --user biopython`` command.




### Performing arbitrary BLASTs using BioPython

Assuming you have BioPython installed and your BLAST database has been created, we'll delve a little into a basic usage of BioPython by BLASTing a sequences from the ``test`` BLAST database against itself.

If you followed the creation of a BLAST database, you would've noticed the default output of ``blastn`` is similar to the pretty textual outputs you would get when running BLAST on the NCBI. However, we would like to interface with these results.  Now, prior to parsing (processing the results) BioPython generates the command line similarly like what was done earlier in this guide. 

<div>Hello</div>
<iframe src="inline/blastquery.html"></iframe>

```python
>>> from Bio.Blast.Applications import NcbiblastnCommandline
>>> bc = NcbiblastnCommandline(query="test.fsa", db="test", evalue=0.001, outfmt=5)
>>> print bc
blastn -outfmt 5 -query test.fsa -db test -evalue 0.001
>>> output, error = bc()
>>> output[:100]
...nDOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.n'

```

### Parsing the results

### Keep a reference of the created BLAST databases in a SQLite database

### Allow the module to be executed as a stand-alone application

