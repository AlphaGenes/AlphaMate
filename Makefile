# General variables
NAME:=AlphaMate
VERSION:= $(shell git rev-parse --short HEAD)
SUBVERSION:=0
PROGRAM:=$(NAME)$(VERSION).$(SUBVERSION)
ALPHAHOUSEDIR:=../AlphaHouse/

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
	FFLAGS := $(FFLAGS) /static /fpp /Qmkl /D $(OSFLAG)
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
	FFLAGS:= $(FFLAGS) -mkl -static-intel -fpp -openmp-link=static -module $(BUILDDIR) -D $(OSFLAG)
	uname := $(shell uname)
	EDDIEFLAGS := $(FFLAGS)
	MAKEDIR := @mkdir -p
	DEL := rm -rf
	# Linux only
	ifeq ($(uname), Linux)
		FFLAGS := $(FFLAGS) -static -static-libgcc -static-libstdc++
	endif
endif

MODS := $(ALPHAHOUSEDIR)AlphahouseMod.f90 \
	$(ALPHAHOUSEDIR)AlphaStatMod.f90 \
	$(ALPHAHOUSEDIR)AlphaEvolveMod.f90 \
	$(ALPHAHOUSEDIR)OrderPackMod.f90 \

# Compile everything
all: directories $(TARGETDIR)$(NAME)$(exe)

eddie: directories Makefile $(MODS) $(SRCDIR)$(NAME).f90
	$(FC) $(MODS) $(SRCDIR)$(NAME).f90 $(EDDIEFLAGS) -o $(TARGETDIR)$(NAME)$(exe)

directories:
	$(MAKEDIR) $(TARGETDIR)
	$(MAKEDIR) $(BUILDDIR)

# Compilation options for debugging
# With warnings about not used variables
debuglong: FFLAGS := $(FFLAGS) -traceback -g -debug all -ftrapuv -fpe0 -warn -check all

debuglong: all

# With memory checks
debug: FFLAGS := $(FFLAGS) -traceback -g -debug all -warn -check bounds -check format \
		-check output_conversion -check pointers -check uninit

debug: all

web: FFLAGS := $(FFLAGS) -D "WEB"

web: all

# If binary is made, intermediate files will be binary
binary: FFLAGS := $(FFLAGS) -D "BINARY"

binary: all
# Compile
$(TARGETDIR)$(NAME)$(exe): Makefile $(MODS) $(SRCDIR)$(NAME).f90
	@echo "Compiling $(NAME)..."
	$(FC) $(MODS) $(SRCDIR)$(NAME).f90 $(FFLAGS) -o $(TARGETDIR)$(NAME)$(exe)
	@echo

# Cleaning
sparklinglyclean: veryclean
	rm -rf TARGETDIR

veryclean: clean
	$(DEL) $(TARGETDIR)$(NAME)$(exe)

clean:
	$(DEL) -rf $(BUILDDIR) *$(obj) *.mod *.dwarf *.i90 *__genmod* *~

.PHONY: make veryclean all
