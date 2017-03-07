SRC:=src/*.cc src/*.h src/external/*

ifdef ORTOOLS
OR_INC:=-I$(ORTOOLS)/include -DHAVE_ORTOOLS
OR_LIB:=-L$(ORTOOLS)/lib -lortools
endif

ets_demo: $(SRC)
	g++ -g -std=c++14 -I src/external $(OR_INC) -O3 -W -Wall src/*.cc src/external/*.cc -o ets_demo $(OR_LIB)

clean:
	rm -f ./ets_demo
