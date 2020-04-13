CC	=	g++
CFLAGS	=	-g
SOURCES = 	source/*.cpp 
SOURCESPYTHON =	apta.cpp dfasat.cpp  refinement.cpp evaluation_factory.cpp random_greedy.cpp  state_merger.cpp parameters.cpp searcher.cpp stream.cpp interactive.cpp 
LFLAGS 	= -w -std=c++11 -L/opt/local/lib -I/opt/local/include -I./source -I./source/evaluation -lm -lpopt -lgsl -lgslcblas -lpthread -ldl
PYTHON_EVAL = source/evaluation/python.cpp

EVALFILES := $(wildcard source/evaluation/*.cpp)
EVALOBJS := $(addprefix source/evaluation/,$(notdir $(EVALFILES:.cpp=.o)))

ifdef WITH_PYTHON
  PYTHON_VERSION=$(shell python3 -c 'import sys; print("".join([str(v) for v in sys.version_info[:2]]))')
  PYTHON_INC=$(shell python3-config --includes)
  PYTHON_LIBS=$(shell python3-config --libs)
  BOOST_LIBS=-lboost_python-py$(PYTHON_VERSION)
else
  EVALFILES := $(filter-out $(PYTHON_EVAL), $(EVALFILES))
  EVALOBJS := $(addprefix source/evaluation/,$(notdir $(EVALFILES:.cpp=.o)))
endif


OUTDIR ?= .

.PHONY: all clean test

all: regen source/gitversion.cpp flexfringe

regen:
	sh collector.sh

debug:
	$(CC) -g $(SOURCES) -o flexfringe $(LFLAGS) $(LIBS)

flexfringe: $(EVALOBJS) source/gitversion.cpp
	$(CC) $(CFLAGS) -o $@ $(SOURCES)  $(EVALOBJS) -I./ $(LFLAGS) $(LIBS)


test: $(EVALOBJS) source/gitversion.cpp
	$(CC) $(FLAGS) -DUNIT_TESTING=1 -I./ -o runtests tests/tests.cpp tests/tail.cpp $(SOURCES) $(EVALOBJS) $(LFLAGS) $(LIBS)
	mkdir -p test-reports
	./runtests -r junit > test-reports/testresults.xml	

source/evaluation/%.o: source/evaluation/%.cpp
	$(CC) -fPIC -c -o $@ $< -I.source $(LFLAGS) $(LIBS) $(PYTHON_INC) $(PYTHON_LIBS) $(BOOST_LIBS) 

clean:
	rm -f flexfringe ./source/evaluation/*.o source/generated.cpp named_tuple.py *.dot *.json exposed_decl.pypp.txt flexfringe*.so source/gitversion.cpp

source/gitversion.cpp: 
	[ -e .git/HEAD ] && [ -e .git/index ] && echo "const char *gitversion = \"$(shell git rev-parse HEAD)\";" > $@ || echo "const char *gitversion = \"No commit info available\";" > $@

python: $(EVALOBJS) source/gitversion.cpp
	export CPLUS_INCLUDE_PATH=/usr/include/python3.5
	$(CC) -fPIC -shared $(CFLAGS)  -o flexfringe.lib.so $(SOURCESPYTHON) $^ -I./ $(LFLAGS) $(LIBS) $(PYTHON_LIBS) $(PYTHON_INC) 
	python3 generate.py
	g++ -W  -g -Wall -fPIC -shared generated.cpp flexfringe.lib.so -o flexfringe.so $(PYTHON_LIBS) $(PYTHON_INC) $(BOOST_LIBS) $(LFLAGS) -Wl,-rpath,'$$ORIGIN' -L. -l:flexfringe.lib.so
	mv flexfringe.lib.so flexfringe.so $(OUTDIR)
