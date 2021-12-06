#$Id: Makefile 23 2012-06-27 10:06:21Z fguillaume $
#=FG===========================================================
# Makefile
#==============================================================
#--------------------------------------------------------------
#MAIN variables
#--------------------------------------------------------------
# Here declare the program name, compilator, Binary and Manual directories are declared
PROG=GLASCOW
#Fortran compiler 
COMP=ifort
#/!\ No space after compiler name Prefer ifort if possible 
#NB :gfortran work also, but ideally if ifort is available we strongly
#recommend you to use it
BINDIR=../BIN
dirMan=../MAN/

#--------------------------------------------------------------------
#Compilations options + Library path
#--------------------------------------------------------------------
#Comment and/or eventually modify if you want to use lapack 
#Define LAPACK and openmp flag depending on compiler
ifeq ($(COMP),gfortran)
	LAPACK=-DLAPACK -llapack 
	OMP= -fopenmp 
else
	LAPACK=-DLAPACK  -mkl
	OMP=-openmp 	
endif
#NB : LAPACK is USE by default, but if it is not set-up on your
#computer just comment (add a # ) before the expression starting with -DLAPACK

FF= -O3 $(OMP) $(LAPACK) -static

# /!\ with ifort -ipo flag should be avoided as it can interfer with openMP barrier
LibMod=
#--------------------------------------------------------------------
#Main rules
#--------------------------------------------------------------------
$(PROG): $(PROG).F90	 	
	$(COMP) $(FF) -o $(PROG) $(PROG).F90
	mv $(PROG) $(BINDIR)/.

#-Cleaning-----------------------------------------------------------
clean:
	rm -f *.o *.a *.mod *~ TAGS

cleanall: clean	
	rm -f ../EX/EX*/*.out
	rm -f ../MAN/*.log ../MAN/*.aux ../MAN/*.blg ../MAN/*.bbg 
	rm -f ../MAN/*.toc ../MAN/*.out 

#-TAGS---------------------------------------------------------------
TAGS:
	etags *.F90 $(LibMod)*.F90  $(LibMod)*.f90 
#-Documentation------------------------------------------------------
doc:
	(cd $(dirMan); pdflatex $(PROG).tex ;pdflatex $(PROG).tex) 
	rm -fr $(dirMan)*.log $(dirMan)*.aux $(dirMan)*.bbg $(dirMan)*.bbg $(dirMan)*.toc 
	rm -fr $(dirMan)*.blg $(dirMan)*.bbl $(dirMan)*.out

docview:
	acroread $(dirMan)$(PROG).pdf

#-Checking for correct compilation and results correctness-----------
check:
	../.CHK/ScriptCHK.sh

#-Archive------------------------------------------------------------
ARCH: clean
	mv ../../trunk ../../$(PROG) #trick to give a correct name to the archive!
	tar -zcvf ../../../$(PROG).tar.gz ../../$(PROG)  --exclude=MAN/$(PROG).tex --exclude=*.svn 
	mv ../../$(PROG)  ../../trunk 

#-Save-a-time-stamped-version-of-the-binary--------------------------
SAVE:
	mv $(BINDIR)/$(PROG) $(BINDIR)/$(PROG)`date +%d%m%y`.bck
