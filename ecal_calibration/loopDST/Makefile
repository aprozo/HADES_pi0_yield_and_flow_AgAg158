APP_NAME      := loopDST
SOURCE_FILES  := loopDST.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := no

include $(HADDIR)/hades.def.mk
LIBS := $(patsubst -ltirpc,-l:libtirpc.so.3,$(LIBS))  # Debian 12 and 13

HYDRA_LIBS    += -lDst
CPP_FLAGS     += -std=c++17

.PHONY:  default
default: clean build install

#HYDRA_LIBS    +=  -lMyLib

include $(HADDIR)/hades.app.mk
