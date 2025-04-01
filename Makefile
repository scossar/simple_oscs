lib.name = oscillators

class.sources = src/triangle~.c src/simple_osc~.c src/cubic_osc~.c src/fold_osc~.c src/simple_phasor~.c src/tri_phase~.c

PDLIBBUILDER_DIR=pd-lib-builder/
include ${PDLIBBUILDER_DIR}/Makefile.pdlibbuilder
