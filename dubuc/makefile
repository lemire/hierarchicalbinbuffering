#
.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h
#
#

VPATH = .:./lemurcore

all: regression benchmark

regression: virtualarray.h externalarray.h transform.cpp dubuccoefficients.h olabuffer.h
	g++ -o regression transform.cpp -g3 -Wall -Winline -I../function

benchmark: virtualarray.h externalarray.h benchmark.cpp dubuccoefficients.h olabuffer.h
	g++ -o benchmark benchmark.cpp -g3 -Wall -Winline -I../function


benchmark1: virtualarray.h externalarray.h benchmark.cpp dubuccoefficients.h olabuffer.h
	g++ -o benchmark1 benchmark.cpp -O2 -g3 -DUSE_EXTERNAL -Wall  -I../function ../lemurcore/lemurcore.a

papibenchmark: virtualarray.h externalarray.h benchmark.cpp dubuccoefficients.h olabuffer.h
	g++ -DDO_PAPI -O2 -o papibenchmark benchmark.cpp -g3 -Wall  -I../function -lpapi -lperfctr

toy: virtualarray.h externalarray.h test.cpp dubuccoefficients.h olabuffer.h
	g++ -o toy test.cpp -g3 -Wall -Winline -I../function


test: regression
	./regression

release: regressionrelease benchmarkrelease

regressionrelease: virtualarray.h externalarray.h transform.cpp dubuccoefficients.h olabuffer.h
	g++ -o regression transform.cpp -O2 -Wall -Winline -I../function

benchmarkrelease: virtualarray.h externalarray.h benchmark.cpp dubuccoefficients.h olabuffer.h
	g++ -o benchmark benchmark.cpp  -O2 -Wall -Winline -I../function #-DNDEBUG

testrelease: regressionrelease
	./regression


.PHONY : clean tags


clean :
	rm -f *.o
	rm -f regression
	rm -f transform
	rm -f benchmark1
	rm -f benchmark
	rm -f *.a
	rm -f tags
