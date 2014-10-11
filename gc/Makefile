all: wolo wolum

wolo: maxflow.o LinkedBlockList.o GCoptimization.o graph.o letsgo.cpp
	g++ -o wolo GCoptimization.o LinkedBlockList.o graph.o maxflow.o letsgo.cpp

wolum: maxflow.o LinkedBlockList.o GCoptimization.o graph.o example.cpp
	g++ -o wolum GCoptimization.o LinkedBlockList.o graph.o maxflow.o example.cpp

maxflow.o: maxflow.cpp
	g++ -g -O -c maxflow.cpp

GCoptimization.o: GCoptimization.cpp 
	g++ -g -O -c GCoptimization.cpp

LinkedBlockList.o: LinkedBlockList.cpp
	g++ -g -O -c LinkedBlockList.cpp

graph.o: graph.cpp
	g++ -g -O -c graph.cpp

clean:
	rm *.o

mathtest: mathtest.cpp
	g++ mathtest.cpp -o testmath
	./testmath