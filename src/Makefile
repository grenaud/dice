CXX      = g++ #-g 


BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)
LIBTABIX = $(realpath tabix/)

CXXFLAGS = -Wall -lm -O3 -lz -std=c++0x -I${LIBGAB} -I${LIBGAB}/gzstream/  -I${BAMTOOLS}/include/ -I${LIBTABIX} -I${BAMTOOLS}/src/  -c
LDFLAGS  = ${LIBTABIX}/libtabix.a ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o  ${BAMTOOLS}/lib/libbamtools.a  -lpthread -lm -lz



all: nucular  nucularMCMC nucularMCMCBAM


${LIBTABIX}/tabix.h:
	rm -rf ${LIBTABIX}/*
	git clone --recursive https://github.com/samtools/tabix.git

${LIBTABIX}/libtabix.a: ${LIBTABIX}/tabix.h
	make -C tabix


${LIBGAB}/utils.h:
	rm -rf ${LIBGAB}/
	git clone --recursive https://github.com/grenaud/${LIBGAB}.git


${LIBGAB}/utils.o: ${BAMTOOLS}/lib/libbamtools.a  ${LIBGAB}/utils.h
	make -C ${LIBGAB}

${BAMTOOLS}/src/api/BamAlignment.h:
	rm -rf ${BAMTOOLS}/
	git clone --recursive https://github.com/pezmaster31/${BAMTOOLS}.git


${BAMTOOLS}/lib/libbamtools.a: ${BAMTOOLS}/src/api/BamAlignment.h
	cd ${BAMTOOLS}/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..


miscfunc.o:	${LIBGAB}/utils.o miscfunc.cpp
	${CXX} ${CXXFLAGS} miscfunc.cpp


nucular.o:	nucular.cpp  libnuc.o
	${CXX} ${CXXFLAGS} nucular.cpp

nucularMCMC.o:	nucularMCMC.cpp libnuc.o
	${CXX} ${CXXFLAGS} nucularMCMC.cpp


nucularMCMCBAM.o:	nucularMCMCBAM.cpp 
	${CXX} ${CXXFLAGS} nucularMCMCBAM.cpp



nucularMCMCBAM:	nucularMCMCBAM.o ${LIBTABIX}/libtabix.a miscfunc.o libnuc.o MistarParser.o  AlleleRecords.o ReadTabix.o SingleAllele.o mistarOperations.o ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a

nucularMCMCBAM:	nucularMCMCBAM.o miscfunc.o libnuc.o MistarParser.o  AlleleRecords.o ReadTabix.o SingleAllele.o mistarOperations.o tabix/libtabix.a ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

libnuc.o:	libnuc.cpp
	${CXX} ${CXXFLAGS} libnuc.cpp


nucular:	nucular.o  ${LIBTABIX}/libtabix.a ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

nucularMCMC:	nucularMCMC.o ${LIBTABIX}/libtabix.a ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f nucular.o nucular nucularMCMC.o nucularMCMC libnuc.o MistarParser.o  AlleleRecords.o ReadTabix.o SingleAllele.o  mistarOperations.o

