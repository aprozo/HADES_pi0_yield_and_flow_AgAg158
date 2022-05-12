APP_NAME      := emc_analysis

SOURCE_FILES  := emc_analysis.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := no

include $(HADDIR)/hades.def.mk

HYDRA_LIBS    += -lDst 


.PHONY:  default
default: clean build install

include $(HADDIR)/hades.app.mk

