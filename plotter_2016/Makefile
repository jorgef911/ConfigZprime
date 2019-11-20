##############################################################
# Author: Andres Florez, Universidad de los Andes, Colombia. #
##############################################################

ObjSuf = o
SrcSuf = cc
ExeSuf =
DllSuf = so
OutPutOpt = -o
HeadSuf = h

ROOTCFLAGS = $(shell root-config --cflags) -O0 -g
ROOTLIBS = $(shell root-config --libs) -O0 -g

# Linux with egcs

CXX = g++ -std=c++11
CXXFLAGS += $(ROOTCFLAGS) -I./ -g

LD = g++ -std=c++11
LDFLAGS += $(ROOTLIBS) -g

SOFLAGS = -shared
LIBS =

SRCDIR = src
OBJDIR = obj

#------------------------------------------------------------------------------
SOURCES = $(wildcard $(SRCDIR)/*.cc)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
#------------------------------------------------------------------------------

all: Plotter


Plotter: $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)	

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

%: $(OBJDIR)/%.o
	$(LD) $(LDFLAGS) -o $@ $< $(LIBS)

clean:
	@echo "Cleaning..."
	@ls $(OBJDIR)
	@rm -f $(OBJECTS)

job: Plotter
	rm *root -f
	./Plotter config/201607.config

.SUFFIXES: .$(SrcSuf) .cc .o .so
