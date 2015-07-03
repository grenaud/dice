CXX      = g++ -g 


BAMTOOLS= ../bamtools/
LIBGAB=  ../libgab/
#LIBTABIX = $(realpath tabix/)

CXXFLAGS = -Wall -lm -O3 -lz -std=c++0x -I${LIBGAB} -I${LIBGAB}/gzstream/  -I${BAMTOOLS}/include/ -Itabix -I${BAMTOOLS}/src/  -c
LDFLAGS  = tabix/libtabix.a ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o  ${BAMTOOLS}/lib/libbamtools.a  -lpthread -lm -lz



all: dice diceBAM BAM2DICE


tabix/tabix.h:
	rm -rf tabix/*
	git clone --recursive https://github.com/samtools/tabix.git

tabix/libtabix.a: tabix/tabix.h
	make -C tabix

${LIBGAB}/gzstream/libgzstream.a: ${LIBGAB}/utils.o

${LIBGAB}/utils.h:
	rm -rf ${LIBGAB}/
	git clone --recursive https://github.com/grenaud/libgab.git ${LIBGAB}


${LIBGAB}/utils.o: ${BAMTOOLS}/lib/libbamtools.a  ${LIBGAB}/utils.h
	make -C ${LIBGAB}

${BAMTOOLS}/src/api/BamAlignment.h:
	rm -rf ${BAMTOOLS}/
	git clone --recursive https://github.com/pezmaster31/bamtools.git ${BAMTOOLS}/


${BAMTOOLS}/lib/libbamtools.a: ${BAMTOOLS}/src/api/BamAlignment.h
	cd ${BAMTOOLS}/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..


miscfunc.o:	${LIBGAB}/utils.o miscfunc.cpp
	${CXX} ${CXXFLAGS} miscfunc.cpp


diceBAM.o:	diceBAM.cpp 
	${CXX} ${CXXFLAGS} diceBAM.cpp


diceBAM:	diceBAM.o miscfunc.o libnuc.o MistarParser.o  AlleleRecords.o ReadTabix.o SingleAllele.o mistarOperations.o tabix/libtabix.a ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

BAM2DICE.o:     BAM2DICE.cpp 
	${CXX} ${CXXFLAGS} BAM2DICE.cpp

BAM2DICE:       BAM2DICE.o miscfunc.o libnuc.o MistarParser.o  AlleleRecords.o ReadTabix.o SingleAllele.o mistarOperations.o tabix/libtabix.a ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

libnuc.o:	libnuc.cpp  ${LIBGAB}/utils.o ${LIBGAB}/gzstream/libgzstream.a
	${CXX} ${CXXFLAGS} libnuc.cpp


dice.o:	dice.cpp libnuc.o
	${CXX} ${CXXFLAGS} dice.cpp

dice:	dice.o tabix/libtabix.a ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f *.o dice  diceBAM BAM2DICE

