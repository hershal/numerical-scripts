#CXX=gccrec g++-code-assist
CXX=g++
CXXFLAGS=-g -O0
LDFLAGS=

SOURCES=gauss_elim.cpp linearcd2D.cpp matrix.cpp program11.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=binary

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o $@
	etags `find . -regex ".*\.[cChH]\(pp\)?" -print`
clean: 
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)
	rm -f TAGS
