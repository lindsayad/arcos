# Alejandro Luque Estepa, 2006
# Margreet Nool 2013

TOPDIR		:= ..

include		$(TOPDIR)/Makefile.config
include		$(TOPDIR)/Makefile.inc

.PHONY:		example2 all clean

all:
	@echo
	@echo "#################################"
	@echo "### BUILDING arcos executable ###"
	@echo "#################################"
	@echo
	make arcos
	mv arcos ..

# ---------------------------------------------------------------------------
DEBUG_LEVEL	= 1

all:	arcos
OBJ = cdr.o poisson.o rz_array.o interpol2.o misc.o \
      cstream.o dft.o grid.o reaction.o mapper.o \
      photo.o main.o sprites.o react_table.o rt.o configuration.o	

LIBDIR1		= $(TOPDIR)/fish90/lib
LIBDIR2		= $(TOPDIR)/arcos_f90/lib
LIBDIR3		= $(TOPDIR)/libconfig-1.4.9/lib
LIBS		= -L./$(LIBDIR1) -L./$(LIBDIR2) $(LIBDIR3)/.libs/libconfig.so -lfish90 -larcos -lm -lfftw3 
INCDIRS1	= ../include
INCDIRS2	= $(TOPDIR)/fish90/include

$(OBJ) : $(LIBDIR1)/libfish90.a	 $(LIBDIR1)/fish.mod \
	$(LIBDIR2)/libarcos.a $(LIBDIR3)/.libs/libconfig.la


%.o:	%.c
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o $@ -c $< 

example2.o :	$(INCDIRS1)/libconfig.h $(INCDIRS1)/configuration.h

example2: configuration.o example2.o
	$(CC) $(CCFLAGS) example2.o configuration.o -o $@ $(LIBS_new)

cdr.o	:	$(INCDIRS1)/cdr.h $(INCDIRS1)/cstream.h \
		$(INCDIRS1)/grid.h $(INCDIRS1)/interpol2.h \
		$(INCDIRS1)/mapper.h $(INCDIRS1)/parameters.h \
		$(INCDIRS1)/poisson.h $(INCDIRS1)/proto.h \
		$(INCDIRS1)/react_table.h $(INCDIRS1)/rz_array.h \
		$(INCDIRS1)/species.h $(INCDIRS1)/tree.h

cstream.o :	$(INCDIRS1)/parameters.h $(INCDIRS1)/proto.h \
		$(INCDIRS1)/species.h

dft.o	:	$(INCDIRS1)/cdr.h $(INCDIRS1)/cstream.h \
		$(INCDIRS1)/grid.h $(INCDIRS1)/parameters.h \
		$(INCDIRS1)/proto.h $(INCDIRS1)/rz_array.h \
		$(INCDIRS1)/species.h

grid.o	:	$(INCDIRS1)/grid.h $(INCDIRS1)/parameters.h \
		$(INCDIRS1)/proto.h $(INCDIRS1)/species.h

interpol2.o :	$(INCDIRS1)/grid.h $(INCDIRS1)/interpol2.h \
		$(INCDIRS1)/parameters.h $(INCDIRS1)/proto.h \
		$(INCDIRS1)/rz_array.h $(INCDIRS1)/species.h

main.o :	$(INCDIRS1)/cdr.h $(INCDIRS1)/cstream.h \
		$(INCDIRS1)/grid.h $(INCDIRS1)/parameters.h \
	       	$(INCDIRS1)/proto.h $(INCDIRS1)/species.h \
	        $(INCDIRS1)/configuration.h	

mapper.o :	$(INCDIRS1)/cstream.h $(INCDIRS1)/grid.h \
		$(INCDIRS1)/interpol2.h $(INCDIRS1)/mapper.h \
		$(INCDIRS1)/proto.h $(INCDIRS1)/species.h

photo.o :	$(INCDIRS1)/cdr.h $(INCDIRS1)/interpol2.h \
		$(INCDIRS1)/mapper.h $(INCDIRS1)/parameters.h \
		$(INCDIRS1)/photo.h $(INCDIRS1)/poisson.h \
		$(INCDIRS1)/proto.h $(INCDIRS1)/rz_array.h \
		$(INCDIRS1)/species.h

poisson.o :	$(INCDIRS1)/cdr.h $(INCDIRS1)/cstream.h \
		$(INCDIRS1)/grid.h $(INCDIRS1)/interpol2.h \
		$(INCDIRS1)/mapper.h $(INCDIRS1)/parameters.h \
		$(INCDIRS1)/poisson.h $(INCDIRS1)/proto.h \
		$(INCDIRS1)/rz_array.h $(INCDIRS1)/species.h \
		$(INCDIRS2)/fishpack.h 

reaction.o :	$(INCDIRS1)/cdr.h $(INCDIRS1)/grid.h \
		$(INCDIRS1)/parameters.h $(INCDIRS1)/photo.h \
		$(INCDIRS1)/proto.h $(INCDIRS1)/react_table.h \
		$(INCDIRS1)/species.h

react_table.o : $(INCDIRS1)/react_table.h

rt.o :		$(INCDIRS1)/parameters.h $(INCDIRS1)/proto.h \
		$(INCDIRS1)/species.h

rz_array.o :	$(INCDIRS1)/parameters.h $(INCDIRS1)/proto.h \
		$(INCDIRS1)/rz_array.h $(INCDIRS1)/species.h

sprites.o :	$(INCDIRS1)/cdr.h $(INCDIRS1)/parameters.h \
		$(INCDIRS1)/photo.h $(INCDIRS1)/poisson.h \
		$(INCDIRS1)/rz_array.h $(INCDIRS1)/species.h

arcos: $(OBJ) 
	$(F90) $(F90FLAGS) $(OBJ) -o $@ $(LIBS)

clean:
	rm -f *.o *~ arcos 2> /dev/null
