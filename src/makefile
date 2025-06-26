####################################################################
# Makefile for  DNS CODE for rectangular duct                      #
# for sekimoto laptop                                              #
# atsushi sekimoto, originally from markus uhlmann                 #
####################################################################
.SUFFIXES: .o .mod .for .F .f90 .h
########################################################## compiler
#include makefile.in_ltrac-hpc-1
#include makefile.in_intel
#include makefile.in_vermeer # for intel + cuda/8.0
#include makefile.in_vermeer_gfortran    # with cuda/8.0
#include makefile.in_vermeer_gfortran_nocuda
#include makefile.in_mac_gfortran_nocuda
include makefile.in_sofia_nocuda
#include makefile.in_ito_fujitsu
#include makefile.in_ito_intel
#include makefile.in_my_ito_intel
#include makefile.in_wst_aquarius
#include makefile.in_rdb_cuda
#include makefile.in_mac_as
##########################################################cpp flags
#use the debug flags if $(DEBUG)='yes':
#ifeq ($(DEBUG),yes)
#   echo "DEBUG = yes"
#   F77FLAGS = $(F77DEBUG)
#   LFLAGS = $(LFDEBUG)
#else
#   echo "DEBUG = no"
   F77FLAGS = $(F77OPT)
   LFLAGS = $(LFOPT)
#endif
##########################################################cpp flags
CPP = /usr/bin/cpp
CPPFLAGS = -P -traditional-cpp \
-DNDEBUG -DNDEBUG_1 -DPROFILO \
-DNSHOOTING -DNMAKE_PRESSURE -DNSLONG \
-DNHELM_SHEN -DNBLOCKING_COMM \
-DTIMING -DNCONVERGENCE -DNFIXED_ZERO_MODE \
-DNFILTER -DDT_FIXED -DCMASSFLOW -DTEMPERATURE \
-DNRECTANGULAR -DNPROJECTIONBC -DNFILTERPLANE \
-DNSET_NOSLIP -DNPUFF_FORCING -DNFORCING_GLOBAL \
-DNMONITORING -DG77 -DFFTW -DNCUBLAS
#
# DEBUG            some output
#
# PROFILO          show time of stages in Runge-Kutta steps
#
# SHOOTING         shoot to find edge state between turbulent and laminar.
#
# MAKE_PRESSURE    remake pressure field from u, v, w 
#                  by solving poisson equation. 
#
#
# HELM_SHEN        use shen's method (galerkin method for large number of nodes)
#   
# LAPLACIAN        ifndef: not to use laplacian3d
#                     ---> less 'transyz, derivyz' but more memory
#
# BLOCKING_COMM    using Blocking communication in parallel.F(shuflling)
#
# TIMING           time the code
#
# DT_FIXED         use fixed time step 
#                   [else: variable, according to specified cfl value]
#
# CMASSFLOW        impose constant mass-flow rate
#                   [else: constant pressure gradient dp/dx]
#
# TEMPERATURE      solve an equation for the temperature (boussinesq)
#
# RECTANGULAR      allow for different number of modes in the cross-plane
#
# G77              use gnu fortran compiler (does not allow for "max"
#                   in parameter statements)
#
# PUFF_FORCING     add a x-localized forcing to the r.h.s.
#                  in order to generate a 'puff'
#
# FORCING_GLOBAL   add global forcing to the r.h.s. 
#                  in order to generate turbulence (may not work by sekimoto) 
# SET_NOSLIP       setting the slip-velocities to zero 
#                  before computing advection terms and so on,  
#                  i.e. before derivatives. ( by sekimoto) 
#
##################################################### objects alpha
INIT = .
SRCDIR = $(INIT)
OBJ=$(INIT)
OBJDIR = $(OBJ)
CALCDIR =$(INIT)
#
OBJECTSM1 = $(OBJ)/mod.o $(OBJ)/main1.o $(OBJECTSL) $(OBJECTSf) $(OBJECTSCUDA)
OBJECTSM2 = $(OBJ)/mod.o $(OBJ)/main1.o $(OBJECTSL) $(OBJECTSf) 
#OBJECTSCONV = $(OBJ)/mod.o $(OBJ)/convh5.o $(OBJECTSL) $(OBJECTSf) $(OBJECTSCUDA)
# most cases does not require GPU for convh5
OBJECTSCONV = $(OBJ)/mod.o $(OBJ)/convh5.o $(OBJECTSL) $(OBJECTSf) 
OBJECTSPARA1 = $(OBJ)/para1.o $(OBJECTSL) $(OBJECTSf)
OBJECTSPARA2 = $(OBJ)/para_qcri.o $(OBJECTSL) $(OBJECTSf)
OBJECTSPUFF = $(OBJ)/puff.o $(OBJECTSL) $(OBJECTSf)
OBJECTSCORR = $(OBJ)/corr.o $(OBJECTSL) $(OBJECTSf)
#
OBJECTSL = \
$(OBJ)/parallel.o\
$(OBJ)/pointer.o\
$(OBJ)/fourier.o\
$(OBJ)/tchebw.o\
$(OBJ)/chebyshev.o\
$(OBJ)/ininewc.o\
$(OBJ)/readwrite.o\
$(OBJ)/readwrite_hdf5.o\
$(OBJ)/operators.o\
$(OBJ)/mult.o\
$(OBJ)/cheb_collocation.o\
$(OBJ)/cpp.o\
$(OBJ)/helm3d.o\
$(OBJ)/int2d3d.o\
$(OBJ)/pseudomat.o\
$(OBJ)/solvediag.o\
$(OBJ)/helm_shen.o\
$(OBJ)/preshen.o\
$(OBJ)/legendre.o\
$(OBJ)/solveshen.o\
$(OBJ)/pseudopois.o\
$(OBJ)/helm3d_temp.o\
$(OBJ)/int2d3d_temp.o\
$(OBJ)/pseudomat_temp.o\
$(OBJ)/solvediag_temp.o


#$(OBJ)/helm3dsv.o\
#$(OBJ)/int2d3dsv.o\
#$(OBJ)/pseudomatsv.o\
#$(OBJ)/solvediagsv.o\
# $(OBJ)/solvepois.o
#
# tcheb.f the extension should be .f
#
OBJECTSf = $(OBJ)/concat.o
#
# parallel_xor.F uses more effective and scalable 
# exclusive-or communication (-DNBLOCKING_COMM is ignored) 
#                       2012/06/05 by sekimoto 
#
#! see /opt/cuda/8.0/src/
OBJECTSCUDA = \
$(OBJ)/fortran.o \
$(OBJ)/cusparse_fortran.o \
$(OBJ)/mytranspose_cuda.o
#
#################################################### compile 
#
#.for.o:
#	cp $< ./tmp/$<
#	$(F77) $(F77FLAGS)  $<      
#
.F.o:
	$(CPP) $(CPPFLAGS)  $< > ./tmp/$*.for 
	$(F77) $(F77FLAGS)  ./tmp/$*.for
##
.f90.o:
	$(CPP) $(CPPFLAGS)  $< > ./tmp/$*.f90 
	$(F77) $(F77FLAGS)  ./tmp/$*.f90
.f.o:
	cp $< ./tmp/$<
	$(F77) $(F77FLAGS)  ./tmp/$<
fortran.o:
	cp fortran.c ./tmp/.
	cp fortran.h ./tmp/.
	cp fortran_common.h ./tmp/.
	$(CC) -c $(CCFLAGS) $(GPUOPT) ./tmp/fortran.c
#       nvcc -c $(GPUFLAGS) $(GPUOPT) ./tmp/fortran.c
cusparse_fortran.o:
	cp cusparse_fortran.c ./tmp/.
	cp cusparse_fortran.h ./tmp/.
	cp cusparse_fortran_common.h ./tmp/.
	$(CC) -c $(CCFLAGS) $(GPUOPT) ./tmp/cusparse_fortran.c
transpose_cuda.o:
	cp transpose_cuda.c ./tmp/.
	cp transpose_cuda.h ./tmp/.
	$(CC) -c $(CCFLAGS) $(GPUOPT) ./tmp/transpose_cuda.c
mytranspose_cuda.o:
	cp mytranspose_cuda.cu ./tmp/.
	cp transpose_cuda.h ./tmp/.
	nvcc -c $(GPUFLAGS) $(GPUOPT) ./tmp/mytranspose_cuda.cu
############################################################ build
#
#
MAIN1   : printmsgA DIM $(OBJECTSM1) 
	@echo Linking... 
	$(F77) -o $(CALCDIR)/MAIN1 $(OBJECTSM1) $(LFLAGS) 
	@echo MAIN1 with CUDA built, enjoy.
##
MAIN2   : printmsgA DIM $(OBJECTSM2) 
	@echo Linking... 
	$(F77) -o $(CALCDIR)/MAIN2 $(OBJECTSM2) $(LFLAGS) 
	@echo MAIN2 without CUDA built, enjoy.
##
convh5	: printmsgA DIM $(OBJECTSCONV)
	@echo Linking...
	$(F77) -o $(CALCDIR)/convh5 $(OBJECTSCONV) $(LFLAGS)
	@echo convh5 built.
#

AVS1    : printmsgA DIM $(OBJECTSAVS)
	@echo Linking...
	$(F77) -o $(CALCDIR)/AVS1 $(OBJECTSAVS) $(LFLAGS)
	@echo AVS1 built.
#
PARA1    : printmsgA DIM $(OBJECTSPARA1)
	@echo Linking...
	$(F77) -o $(CALCDIR)/PARA1 $(OBJECTSPARA1) $(LFLAGS)
	@echo PARA1 built for visualizing all in paraview hdf5 format.
#
PARA2    : printmsgA DIM $(OBJECTSPARA2)
	@echo Linking...
	$(F77) -o $(CALCDIR)/PARA2 $(OBJECTSPARA2) $(LFLAGS)
	@echo PARA2 built for visualizing lap.P in paraview hdf5 format.
###
PUFF    : printmsgA DIM $(OBJECTSPUFF)
	@echo Linking...
	$(F77) -o $(CALCDIR)/PUFF $(OBJECTSPUFF) $(LFLAGS)
	@echo PUFF built.
#
CORR    : printmsgA DIM $(OBJECTSCORR)
	@echo Linking...
	$(F77) -o $(CALCDIR)/CORR $(OBJECTSCORR) $(LFLAGS)
	@echo CORR built.
#
BLAS    : 
	( cd $(BLASDIR); \
	make -f makefile_blas clean; \
	make -f makefile_blas LIB DEBUG=$(DEBUG); )
#
#
DIM   : printmsgA
#	$(CPP) $(CPPFLAGS) ctes3D_in > ctes3D
#
clean:
	rm -f *.o *.mod *.for TEST* AVS* ctes3D
	rm -r ./tmp/*
	cp concat.F concat.for
########################################################## message
printmsgA :
	@echo
	@echo Building POST_PROCESSOR for $(SYSTEM)
	@echo Compiler flags  : $(F77FLAGS)
	@echo Compiler        : $(F77)
	@echo Linker   flags  : $(LFLAGS)
	@echo cpp flags       : $(CPPFLAGS)
#
###################################################### end of file

