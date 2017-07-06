CXX:=g++

PYETS_DEP:=src/ets.cc src/pyETS.cc src/*.h src/external/*
PYETS_SRC:=src/ets.cc src/pyETS.cc src/external/*.cc
PYCCONV_DEP:=src/cconv.cc src/pycconv.cc src/*.h src/external/*
PYCCONV_SRC:=src/cconv.cc src/pycconv.cc src/external/*.cc

PYBIND_INC:=src
PYFLAGS:=$(filter-out -Wstrict-prototypes,$(shell python-config --cflags --ldflags))

COMMON_FLAGS:=-O3 -g -fmax-errors=1 -std=c++14 -I src/external -W -Wall

pyETS: $(PYETS_DEP)
	$(CXX) $(COMMON_FLAGS) -shared -fPIC -Isrc $(PYETS_SRC) -o pyETS.so $(PYFLAGS)

pycconv: $(PYCCONG_DEP)
	$(CXX) $(COMMON_FLAGS) -shared -fPIC -Isrc $(PYCCONV_SRC) -o pycconv.so $(PYFLAGS)

clean:
	rm -f pyETS.so
