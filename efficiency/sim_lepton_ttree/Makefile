APP_NAME      := emc_analysis

SOURCE_FILES  := emc_analysis.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := no

include $(HADDIR)/hades.def.mk
LIBS := $(patsubst -ltirpc,-l:libtirpc.so.3,$(LIBS))  # Debian 12 and 13

HYDRA_LIBS    += -lDst 


.PHONY:  default
default: clean build install

include $(HADDIR)/hades.app.mk

