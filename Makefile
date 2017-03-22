DEMO_DEP:=src/astronomy.cc src/ets_demo.cc src/ets_assigners.cc src/ets_tools.cc src/*.h src/external/*
DEMO_SRC:=src/astronomy.cc src/ets_demo.cc src/ets_assigners.cc src/ets_tools.cc src/external/*.cc
PYETS_DEP:=src/astronomy.cc src/pyETS.cc src/ets_assigners.cc src/ets_tools.cc src/*.h src/external/*
PYETS_SRC:=src/astronomy.cc src/pyETS.cc src/ets_assigners.cc src/ets_tools.cc src/external/*.cc

PYBIND_INC:=src
PYFLAGS:=$(filter-out -Wstrict-prototypes,$(shell python-config --cflags --ldflags))

COMMON_FLAGS:=-O3 -g -fmax-errors=1 -std=c++14 -I src/external -W -Wall

ifdef ORTOOLS
OR_INC:=-I$(ORTOOLS)/include -DHAVE_ORTOOLS
OR_LIB:=-L$(ORTOOLS)/lib -lortools
endif

ets_demo: $(DEMO_DEP)
	g++ $(COMMON_FLAGS) $(OR_INC) $(DEMO_SRC) -o ets_demo $(OR_LIB)

pyETS: $(PYETS_DEP)
	g++ $(COMMON_FLAGS) -shared -fPIC -Isrc $(OR_INC) $(PYETS_SRC) -o pyETS.so $(OR_LIB) $(PYFLAGS)

clean:
	rm -f ./ets_demo pyETS.so
