#=======================================================================
#                   define the compiler names
#=======================================================================
CC       = gcc 
F90      = gfortran
# F90      = ifort
PYTHON = /usr/bin/env python
PY2SO = ../../py2so.py
BUILD = .
#=======================================================================
#                     additional flags
#=======================================================================
ifeq ($(F90),gfortran)
	FPP      = gfortran -E
	FPP_F90FLAGS = -x f95-cpp-input -fPIC
	F90FLAGS = -fPIC -shared 
    FCOMP    = gfortran
    LIBS     =
endif

ifeq ($(F90),ifort)
	FPP      = gfortran -E # gfortran f90wrap temp files only. not compilation
	FPP_F90FLAGS = -x f95-cpp-input -fPIC -cpp -DMPI
	F90FLAGS = -fpscomp logicals -fPIC -fpp -DMPI # use 1 and 0 for True and False
    FCOMP    = intelem # for f2py
	MPIFC = mpiifort
    LIBS =
endif

CFLAGS = -fPIC  
ifeq (${UNAME}, Darwin)
  LIBTOOL = libtool -static -o
else
  LIBTOOL = ar src
endif

# Python module name
PYTHON_MODN = f90LegacyFingerPrint
# mapping between Fortran and C types
KIND_MAP = kindMap.json
LIBSRC_SOURCES = Kinds Math Constants BondCrt BondCharMatrix


# file names
LIBSRC_FILES = $(addsuffix .F90,${LIBSRC_SOURCES})

# object files
LIBSRC_OBJECTS = $(addsuffix .o,${LIBSRC_SOURCES})

# only used when cleaning up
LIBSRC_FPP_FILES = $(addsuffix .fpp,${LIBSRC_SOURCES})

.SUFFIXES: .F90 .fpp

.PHONY: all clean

all: LIB/${PYTHON_MODN}.*.so LIB/_${PYTHON_MODN}.*.so

LIB/${PYTHON_MODN}.*.so: ${PYTHON_MODN}.py
	sed -i'.bak' 's/import _$(PYTHON_MODN)/from . import _$(PYTHON_MODN)/g' $(PYTHON_MODN).py
	$(PYTHON) $(PY2SO) $^ LIB

${PYTHON_MODN}.py: ${LIBSRC_FPP_FILES}
	f90wrap -m ${PYTHON_MODN} \
			-k ${KIND_MAP} \
			./FPP/*
	mkdir -p F90wrap
	mv f90wrap* ./F90wrap/

LIB/_${PYTHON_MODN}.*.so: libsrc.a ${LIBSRC_FPP_FILES}
	f2py-f90wrap -c -m _${PYTHON_MODN} \
				--fcompiler=$(FCOMP) --f90exec=$(F90) --f77exec=$(F90) \
				-L. -lsrc \
				--opt=-O3 \
				--debug \
				--build-dir . \
				./F90wrap/f90wrap*.f90 || cp ./SRC/fortranobject.c ./src*/src*/. && \
	f2py-f90wrap -c -m _${PYTHON_MODN} \
				--fcompiler=$(FCOMP) --f90exec=$(F90) --f77exec=$(F90) \
				-L. -lsrc \
				--opt=-O3 \
				--debug \
				--build-dir . \
				./F90wrap/f90wrap*.f90
	mkdir -p LIB
	mv $(notdir $@) LIB/
# mv ._*.so to LIB/_*.so

libsrc.a: ${LIBSRC_OBJECTS}
	${LIBTOOL} $@ ./OBJ/*
	
${LIBSRC_OBJECTS}: %.o:./SRC/%.F90
	mkdir -p OBJ F90wrap
	${F90} ${F90FLAGS} -c $< -o ./OBJ/$@ -J ./F90wrap -O3

.c.o:
	${CC} ${CFLAGS} -c $< -o $@

${LIBSRC_FPP_FILES}: %.fpp : ./SRC/%.F90
	mkdir -p FPP
	${FPP} ${FPP_F90FLAGS} $< -o ./FPP/$@

clean:
	-rm libsrc.a ${PYTHON_MODN}.py ${PYTHON_MODN}.py.bak ${PYTHON_MODN}.c
	-rm .f2py_f2cmap
	-rm -r F90wrap FPP OBJ LIB src*
	-rm -r .libs
