# This file defines common Makefile variables and rules that are used
# to compile the FISHPACK libraries. Normally, this
# file does not have to be modified; see the file Makefile.conf.

# The names of the fishpack libraries.

FISHPACK_LIB	= libfishpack.a
FISH90_LIB	= libfish90.a

# Commands to compile and link C and Fortran 90 code.

Compile		= $(CC) $(CFLAGS) $(DFLAGS) $(addprefix -I,$(INCDIRS)) \
		        -o $@ -c $<

LinkExe		= $(CC) $(addprefix -L,$(LIBDIRS)) -o $@ $^ \
		        $(F90FLAGS) $(addprefix -l,$(LIBS))

CompileF90	= $(F90) $(F90FLAGS) $(addprefix -I,$(INCDIRS)) \
		         -o $@ -c $<

LinkF90Exe	= $(F90) $(addprefix -L,$(LIBDIRS)) -o $@ $^ \
		         $(F90FLAGS) $(addprefix -l,$(LIBS))

CompileF	= $(F90) $(F90FLAGS) $(addprefix -I,$(INCDIRS)) \
		         -o $@ -c $<

Compilef	= $(F90) $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
		         -o $@ -c $<

# Implicit rules for compiling C and Fortran source files.

%.o:		%.c
		$(Compile)

%.o:		%.f90
		$(CompileF90)

%.o:		%.F
		$(CompileF)

%.o:		%.f
		$(Compilef)

