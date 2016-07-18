EXTSRC:=src/external/*

ets_demo: src/ets_demo.cc $(EXTSRC)
	g++ -g -std=c++14 -I src/external -O3 -W -Wall src/ets_demo.cc src/external/*.cc -o ets_demo

clean:
	rm -f ./ets_demo
