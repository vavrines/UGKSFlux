
#-------------------------------------------------------------------------#
#                                                                         #
#                         Cloud CFD Version 0.1                           #
#                                                                         #
#        CFD Software Package in Scientific and Engineering Research      #
#                                                                         #
# 					   Developed by Tianbai Xiao (2019)                   #    
#                                                                         #
#                         Makefile for software                           #        
#                                                                         #
#-------------------------------------------------------------------------#

SHELL = /bin/sh
.SUFFIXES: .f90 .FPP
%.o: %.mod

#--------------------------------------------------
# variables
#--------------------------------------------------

# directories
BINDIR = bin
SRCDIR = src

# compiler and flags
FCC = ifort
FCP = mpifort
OMP = off
MPI = off

# commands for example
NPS = 2
CMD = ../$(BINDIR)/CloudCFD
RUN = mpirun

FC = $(FCC)

ifeq ($(FCC),gfortran)
    FCFLAGS += -ffree-form -ffree-line-length-none -DDOUBLE -J $(BINDIR) -O3
    OMPFLAG = -fopenmp
endif

ifeq ($(FCC),ifort)
    FCFLAGS += -free -DDOUBLE -module $(BINDIR) -O3 -xHOST -no-prec-div
    OMPFLAG = -openmp -parallel
endif

ifeq ($(OMP),on)
    FCFLAGS += $(OMPFLAG)
endif

# executable and objects
EXEC = CloudCFD
SRCS = data.o grid.o preprocess.o kinetic.o math.o postprocess.o \
	   initialize.o boundary.o UGKSFlux1D.o loop.o main.o

#--------------------------------------------------
# compiling
#--------------------------------------------------

all: checkdir CloudCFD

checkdir:
	@mkdir -p $(BINDIR)

CloudCFD: checkdir $(BINDIR)/CloudCFD

$(BINDIR)/CloudCFD: $(addprefix $(BINDIR)/, $(SRCS))
	$(FC) $(FCFLAGS) -o $@ $(addprefix $(BINDIR)/, $(SRCS))

$(BINDIR)/%.o: $(SRCDIR)/%.???
	$(FC) $(FCFLAGS) -c $< -o $@

#--------------------------------------------------
# other functions
#--------------------------------------------------

clean:
	rm -f $(BINDIR)/*

veryclean:
	rm -f *.hst *.dat	
