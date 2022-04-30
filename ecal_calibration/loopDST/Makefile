APP_NAME      := loopDST
SOURCE_FILES  := loopDST.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := no

include $(HADDIR)/hades.def.mk

HYDRA_LIBS    += -lDst
CPP_FLAGS     += -std=c++14

.PHONY:  default
default: clean build install

#HYDRA_LIBS    +=  -lMyLib

include $(HADDIR)/hades.app.mk
