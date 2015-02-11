
CXX      = g++   

BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)

CXXFLAGS = -Wall -lm -O3 -lz -I${LIBGAB} -I${LIBGAB}/gzstream/  -c
LDFLAGS  = -lz


all: nucular  nucularMCMC

nucular.o:	nucular.cpp libnuc.o
	${CXX} ${CXXFLAGS} nucular.cpp

nucularMCMC.o:	nucularMCMC.cpp libnuc.o
	${CXX} ${CXXFLAGS} nucularMCMC.cpp

libnuc.o:	libnuc.cpp
	${CXX} ${CXXFLAGS} libnuc.cpp


nucular:	nucular.o ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

nucularMCMC:	nucularMCMC.o ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f nucular.o nucular nucularMCMC.o nucularMCMC

