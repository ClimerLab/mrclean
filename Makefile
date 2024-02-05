SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#---------------------------------------------------------------------------------------------------
#
# Set CPLEXDIR and CONCERTDIR to the directories where CPLEX and CONCERT are installed.
#
#---------------------------------------------------------------------------------------------------

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio221/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio221/concert

#---------------------------------------------------------------------------------------------------
# Compiler selection
#---------------------------------------------------------------------------------------------------

CXX = g++

#---------------------------------------------------------------------------------------------------
# Directories
#---------------------------------------------------------------------------------------------------

OBJDIR = build
SRCDIR = src

#---------------------------------------------------------------------------------------------------
# Executables
#---------------------------------------------------------------------------------------------------

EXE = mrclean

#---------------------------------------------------------------------------------------------------
# Object files
#---------------------------------------------------------------------------------------------------

OBJ = BinContainer.o GreedySolver.o RowColSolver.o ElementSolver.o Timer.o ConfigParser.o \
			CleanSolution.o
ALL_OBJ = $(OBJ) clean.o

#---------------------------------------------------------------------------------------------------
# Compiler options
#---------------------------------------------------------------------------------------------------

CXXFLAGS = -O3 -Wall -fPIC -fexceptions -DIL_STD -std=c++11 -fno-strict-aliasing

#---------------------------------------------------------------------------------------------------
# Link options and libraries
#---------------------------------------------------------------------------------------------------

CPLEXLIBDIR    = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR  = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CXXLNDIRS      = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CXXLNFLAGS     = -lconcert -lilocplex -lcplex

#---------------------------------------------------------------------------------------------------
# Includes
#---------------------------------------------------------------------------------------------------

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

INCLUDES = -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

#---------------------------------------------------------------------------------------------------
all: CXXFLAGS += -DNDEBUG
all: $(EXE)

debug: CXXFLAGS += -g
debug: $(EXE)

mrclean: $(addprefix $(OBJDIR)/, clean.o)
	$(CXX) $(CXXLNDIRS) -o $@ $(addprefix $(OBJDIR)/, $(ALL_OBJ)) $(CXXLNFLAGS)

$(OBJDIR)/clean.o:	$(addprefix $(SRCDIR)/, main.cpp) \
										$(addprefix $(OBJDIR)/, $(OBJ))
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJDIR)/RowColSolver.o:	$(addprefix $(SRCDIR)/, RowColSolver.cpp RowColSolver.h) \
				   							$(addprefix $(OBJDIR)/, BinContainer.o)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJDIR)/ElementSolver.o:	$(addprefix $(SRCDIR)/, ElementSolver.cpp ElementSolver.h) \
				   							$(addprefix $(OBJDIR)/, BinContainer.o)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJDIR)/GreedySolver.o:	$(addprefix $(SRCDIR)/, GreedySolver.cpp GreedySolver.h) \
													$(addprefix $(SRCDIR)/, MrCleanUtils.h) \
				   								$(addprefix $(OBJDIR)/, BinContainer.o)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJDIR)/BinContainer.o: $(addprefix $(SRCDIR)/, BinContainer.cpp BinContainer.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/Timer.o: $(addprefix $(SRCDIR)/, Timer.cpp Timer.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/ConfigParser.o: $(addprefix $(SRCDIR)/, ConfigParser.cpp ConfigParser.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/CleanSolution.o: $(addprefix $(SRCDIR)/, CleanSolution.cpp CleanSolution.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#---------------------------------------------------------------------------------------------------
.PHONY: clean
clean:
	/bin/rm -f $(OBJDIR)/*.o
#---------------------------------------------------------------------------------------------------