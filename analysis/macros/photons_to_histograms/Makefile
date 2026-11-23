APP_NAME      := photons_to_histograms

CPP_FLAGS     += -std=c++17

SOURCE_FILES  := main.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := no

include $(HADDIR)/hades.def.mk

HYDRA_LIBS    += -lDst 


.PHONY:  default
default: clean build install

include $(HADDIR)/hades.app.mk

