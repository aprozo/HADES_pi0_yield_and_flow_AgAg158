APP_NAME      := exp_sim_to_photons

SOURCE_FILES  := main.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := no

include $(HADDIR)/hades.def.mk

HYDRA_LIBS    += -lDst
CPP_FLAGS     += -std=c++14

.PHONY:  default
default: clean build install

include $(HADDIR)/hades.app.mk
