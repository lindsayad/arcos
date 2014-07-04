include 	Makefile.inc
include 	Makefile.config

SOURCE_DIR	= src
FISH90_DIR	= fish90/src
F90_DIR		= arcos_f90/src
CONFIG_DIR	= libconfig-1.4.9/lib

SUBDIRS		= $(F90_DIR) $(FISH90_DIR) $(SOURCE_DIR) $(CONFIG_DIR)
CLEANDIRS	= $(SUBDIRS:%=clean-%)

.PHONY:		all allclean $(SUBDIRS) $(CLEANDIRS)

# Default Target is to build everything
all: $(SUBDIRS)
allclean: $(CLEANDIRS) clean

$(SUBDIRS):
	$(MAKE) -C $@

$(CLEANDIRS):
	$(MAKE) -C $(@:clean-%=%) clean

clean:
	$(RM) $(SOURCE_DIR)/*.o ./arcos


