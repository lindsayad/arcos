# This Makefile builds the FISH90 library.
# FISH90 is a modified F90 version of FISHPACK

TOPDIR		:= ../..

include 	$(TOPDIR)/Makefile.config
include 	$(TOPDIR)/Makefile.inc

LIBFISH		= ../lib/libfish90.a

.PHONY:		all

all:
	@echo
	@echo "###############################"
	@echo "### BUILDING FISH90 library ###"
	@echo "###############################"
	@echo
	make lib

SRCF		= genbunal.f90 gnbnaux.f90 hst_procs.f90 poisson.f90
SRCC		= hstcrt_wrap.c hstcyl_wrap.c

OBJF		= $(subst .f90,.o,$(SRCF))
OBJC		= $(subst .c,.o,$(SRCC))
OBJ		:= fish.o $(OBJF) $(OBJC) 

hst_procs.o :	fish.o genbunal.o poisson.o
genbunal.o :	gnbnaux.o poisson.o
poisson.o :	gnbnaux.o

lib :		$(LIBFISH)

$(LIBFISH) :	$(OBJ)
	$(AR) $@ $? 
	mv fish.mod ../lib

clean:
	rm -f $(LIBFISH) $(OBJ) *~ ../lib/fish.mod ../lib/FISH.mod fish.f
