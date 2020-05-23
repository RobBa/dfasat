# README #

flexfringe (formerly DFASAT), a flexible state-merging framework written in C++.

### What is this repository for? ###

You can issue pull requests for bug fixes and improvements. Most work will happen in the development branch while master contains a more stable version.

### How do I get set up? ###

flexfringe has one required dependency: libpopt for argument parsing. Some heuristic functions bring their own dependencies. We provide an implementation of a likelihood-based merge function for probabilistic DFAs. It needs the GNU scientific library (development) package (e.g. the libgsl-dev package in Ubuntu).
 
If you want to use the reduction to SAT and automatically invoke the SAT solver, you need to provide the path to the solver binary. flexfringe has been tested with lingeling (which you can get from http://fmv.jku.at/lingeling/ and run its build.sh).
**PLEASE NOTE:** SAT solving only works for learning plain DFAs. The current implementation is not verified to be correct. Use an older commit if you rely on SAT-solving.

You can build and compile the flexfringe project by running

$ make clean all

in the main directory to build the executable named *flexfringe*.


### How do I run it? ###

Run ./flexfringe --help to get help.

The start.sh script together with some .ini files provides a shortcut to storing 

Example:

`$ ./start.sh ini/batch-overlap.ini data/staminadata/1_training.txt.dat`

See the .ini files for documentation of parameter flags. 

#### Input files ####

The default input is formated following the Abadingo formating:

```
num_samples alphabet_size
label length sym1 sym2 ... symN
.
.
.
```
for each symbol, additional data can be attached via /, i.e. `label length sym1/data1 sym2/data2 ... symN/dataN`. These can represent outputs (e.g. for Mealy or Moore machines), or any other information needed by a custom evaluation function.

Real-valued attributes, e.g. for real-time automata, can be attached via :, i.e. `label length sym1:real1,real2,realn ...`. The number of attributes has to be specified in the header after the alphabet size, i.e. `num_samples alphabet_size:num_attributes`.

#### Output files ####

flexfringe will generate several .dot files into the specified output directory (./ by default):

* pre\:\*.dot are intermediary dot files created during the merges/search process.
* final.dot is the end result

You can plot the dot files via

`$ dot -Tpdf file.dot -o outfile.pdf`
or
`$ ./show.sh final.dot`

after installing dot from graphviz.

### Contribution guidelines ###

* Fork and implement, request pulls.
* You can find sample evaluation files in ./source/evaluation. Make sure to REGISTER your own file to be able to access it via the -h flag.

#### Writing tests ####

Unit tests are incomplete. *flexfringe* uses the Catch2 framework (see the [https://github.com/catchorg/Catch2/blob/master/docs/tutorial.md](Tutorial) and the *tests* folder for some examples.

#### Logging ####
Logging is incomplete. *flexfringe* uses the loguru framework (see the [https://github.com/emilk/loguru/blob/master/README.md](Loguru documentation)). *flexfringe* uses the stream-version. Please log using the `LOG_S(LEVEL) << "messge"` syntax to implement logging.
 
### Who do I talk to? ###

* Sicco Verwer (original author; best to reach out to for questions on batch mode, RTI+ implementation, and SAT reduction)
* Christian Hammerschmidt (author of the online/streaming mode, interactive mode, and the flexible evaluation function mechanism)
* Sofia Tsoni (scientific programmer, maintainer)

### Credits and Licences

*flexfinge* relies on a number of open source packages and libraries. You can find the respective LICENCE files in the source/utility subdirectory. Most notable, we use
* CLI11 for command line parsing
* Catch for unit testing
