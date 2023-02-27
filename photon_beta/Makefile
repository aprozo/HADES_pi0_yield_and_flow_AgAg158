APP_NAME      := photons_to_histograms

SOURCE_FILES  := emc_diphoton_ana.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := no

include $(HADDIR)/hades.def.mk

.PHONY:  default
default: clean build install

include $(HADDIR)/hades.app.mk
