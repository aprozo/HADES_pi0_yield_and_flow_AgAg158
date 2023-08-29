APP_NAME      := emc_diphoton_ana

SOURCE_FILES  := emc_diphoton_ana.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := no

include $(HADDIR)/hades.def.mk

#HYDRA_LIBS    += -lDst 
CPP_FLAGS     += -std=c++11

.PHONY:  default
default: clean build install

include $(HADDIR)/hades.app.mk

