#
# testje: ifc -O3 -> user    4m50.642s
#         g95 -O3 -> user    4m56.303s
#
F90 = f95    # latest vs g95 (usr15)
#F90 = ../g95.1  # versie jun 06 (usr16)
#F90 = ifort


F77FILES = \
 solver_sparsekit2.f\
 lapacks.f  
# solver_bcg2.f
# solver_ma27.f

F90FILES=\
 modules.f90\
 patches.f90\
 particles.f90\
 initialfield.f90\
 saverestart.f90\
 artificial.f90\
 solverinterface.f90\
 readcontrolfile.f90\
 patchscalars.f90\
 collectcells.f90\
 user.f90\
 tools.f90\
 opendx.f90\
 vtk.f90\
 tecplt.f90\
 gmv.f90\
 dolfyn.f90

F77OBJS = ${F77FILES:.f=.o}
F90OBJS = ${F90FILES:.f90=.o}

##################################################################
# Targets
##################################################################

dolfyn:	${F90OBJS} ${F77OBJS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} ${F77OBJS} ${LFLAGS} ${INCS}  -o dolfyn 
clean:
	@rm  *.o

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 .f 


.f90.o:
#	${F90} -O0 -g -Wall -c ${INCS} ${@:.o=.f90}
#	${F90} -O2 -warn alignments,argument_checking,errors,uninitialized,usage -traceback -c ${INCS} ${@:.o=.f90}
	${F90} -O2  -c ${INCS} ${@:.o=.f90}
.f.o:
	${F90} -O2 -c ${INCS} ${@:.o=.f}
