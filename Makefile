# General variables
NAME:=AlphaMate
VERSION:= $(shell git rev-parse --short HEAD)
SUBVERSION:=0
PROGRAM:=$(NAME)$(VERSION).$(SUBVERSION)

# Set the default compiler to iFort
FC:=ifort
FFLAGS:=-O3 -DVERS=""commit-$(VERSION)""

#  If -D WEB is specified, stops will be put into alphasim.

# MS Windows
ifeq ($(OS), Windows_NT)
	SRCDIR      := src/
	BUILDDIR    :=
	TARGETDIR   :=
	OSFLAG := "OS_WIN"
	## see also https://software.intel.com/en-us/compiler_winapp_f (2014-12-03)
	FFLAGS := $(FFLAGS) /static /fpp /Qmkl /D $(OSFLAG) # /Qparallel
	obj := .obj
	MAKEDIR :=
	exe := .exe
	CC := cl
	CFLAGS := /EHsc
	DEL := del
else
	# Linux or Mac OSX
	SRCDIR      := src/
	BUILDDIR    := objs/
	TARGETDIR   := bin/
	obj := .o
	OSFLAG := "OS_UNIX"
	# TODO: can we make this generic?
	MKLROOT := /opt/intel/mkl
	# On Eddie
	# MKLROOT:=/exports/applications/apps/intel/ClusterStudio2013/mkl
	MKLLIB := -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm
	MKLINC := -I$(MKLROOT)/include
	exe :=
	FFLAGS:= $(FFLAGS) -mkl -static-intel -fpp -openmp-link=static -module $(BUILDDIR) -D $(OSFLAG) # -parallel
	uname := $(shell uname)
	MAKEDIR := @mkdir -p
	DEL := rm -rf
	# Linux only
	ifeq ($(uname), Linux)
		FFLAGS := $(FFLAGS) -static -static-libgcc -static-libstdc++
	endif
endif

# Compile everything
all: directories $(TARGETDIR)$(NAME)$(exe) $(TARGETDIR)AlphaMate$(exe)

directories:
	$(MAKEDIR)  $(TARGETDIR)
	$(MAKEDIR)  $(BUILDDIR)

# Compilation options for debugging
# With warnings about not used variables
debuglong: FFLAGS := $(FFLAGS) -traceback -g -debug all -fpp -ftrapuv -fpe0 -warn -check all

debuglong: all

# With memory checks
debug: FFLAGS := $(FFLAGS) -traceback -g -debug all -warn -check bounds -check format \
		-check output_conversion -check pointers -check uninit -fpp

debug: all

parallel: FFLAGS := $(FFLAGS) -par-report1 -guide-par

parallel: all

web: FFLAGS := $(FFLAGS) -D "WEB"

web: all

# If binary is made, intermediate files will be binary
binary: FFLAGS := $(FFLAGS) -D "BINARY"

binary: all
# Compile AlphaMate
$(TARGETDIR)AlphaMate$(exe): $(SRCDIR)AlphaMate.f90
	@echo "Compiling AlphaMate..."
	$(FC) $(SRCDIR)AlphaMate.f90 $(FFLAGS) -o $(TARGETDIR)AlphaMate$(exe)
	@echo

# Cleaning
sparklinglyclean: veryclean
	rm -rf TARGETDIR

veryclean: clean
	$(DEL) $(TARGETDIR)AlphaMate$(exe)

clean:
	$(DEL) -rf $(BUILDDIR) *$(obj) *.mod *.dwarf *.i90 *__genmod* *~

.PHONY: make veryclean all
