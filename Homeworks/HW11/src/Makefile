#CXX=gccrec g++-code-assist
CXX=g++
CXXFLAGS=-g -O0
LDFLAGS=

SOURCES=$(wildcard *.cpp) $(wildcard *.c)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=binary

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o $@
	etags `find . -regex ".*\.[cChH]\(pp\)?" -print`
clean: 
	rm $(EXECUTABLE)
	rm $(OBJECTS)
	rm TAGS
