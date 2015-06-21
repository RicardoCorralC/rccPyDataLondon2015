all: TMalign

TMalign: TMalign.cpp basic_fun.h  Kabsch.h  NW.h TMalign.h

	g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp 
clean:
	rm -f TMalign
