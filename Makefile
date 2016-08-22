SRC:=src/*.cc src/*.h src/external/*

ets_demo: $(SRC)
	g++ -g -std=c++14 -I src/external -O3 -W -Wall src/*.cc src/external/*.cc -o ets_demo

clean:
	rm -f ./ets_demo
