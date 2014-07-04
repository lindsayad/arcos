
TOPDIR	:= ../..

include	$(TOPDIR)/Makefile.config
include	$(TOPDIR)/Makefile.inc

.PHONY:	all

all:	lib

SRCF	= cdr_advect_diffu_vec.f90 
SRCC	= cdr_advect_diffu_wrap.c

INCDIRS	= ../../include

OBJF	=$(subst .f90,.o,$(SRCF))
OBJC	=$(subst .c,.o,$(SRCC))
OBJ := $(OBJF) $(OBJC)

LIB	=../lib/libarcos.a

lib : $(LIB)

$(LIB) : $(OBJ)
	$(AR) $@ $? 

clean:
	rm -f $(LIB) $(OBJ) *~
