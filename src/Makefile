CXX      = g++ -g 


BAMTOOLS= $(realpath ../bamtools/)
#LIBGAB= $(realpath ../libgab/)
#LIBTABIX = $(realpath tabix/)

CXXFLAGS = -Wall -lm -O3 -lz -std=c++0x -Ilibgab -Ilibgab/gzstream/  -I${BAMTOOLS}/include/ -Itabix -I${BAMTOOLS}/src/  -c
LDFLAGS  = tabix/libtabix.a ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o  ${BAMTOOLS}/lib/libbamtools.a  -lpthread -lm -lz



all: nucular  nucularMCMC nucularMCMCBAM


tabix/tabix.h:
	rm -rf tabix/*
	git clone --recursive https://github.com/samtools/tabix.git

tabix/libtabix.a: tabix/tabix.h
	make -C tabix


libgab/utils.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git


libgab/utils.o: ${BAMTOOLS}/lib/libbamtools.a  libgab/utils.h
	make -C libgab

${BAMTOOLS}/src/api/BamAlignment.h:
	rm -rf ${BAMTOOLS}/
	git clone --recursive https://github.com/pezmaster31/${BAMTOOLS}.git


${BAMTOOLS}/lib/libbamtools.a: ${BAMTOOLS}/src/api/BamAlignment.h
	cd ${BAMTOOLS}/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..


miscfunc.o:	libgab/utils.o miscfunc.cpp
	${CXX} ${CXXFLAGS} miscfunc.cpp


nucular.o:	nucular.cpp  libnuc.o
	${CXX} ${CXXFLAGS} nucular.cpp

nucularMCMC.o:	nucularMCMC.cpp libnuc.o
	${CXX} ${CXXFLAGS} nucularMCMC.cpp


nucularMCMCBAM.o:	nucularMCMCBAM.cpp 
	${CXX} ${CXXFLAGS} nucularMCMCBAM.cpp



nucularMCMCBAM:	nucularMCMCBAM.o tabix/libtabix.a miscfunc.o libnuc.o MistarParser.o  AlleleRecords.o ReadTabix.o SingleAllele.o mistarOperations.o libgab/utils.o  libgab/gzstream/libgzstream.a

nucularMCMCBAM:	nucularMCMCBAM.o miscfunc.o libnuc.o MistarParser.o  AlleleRecords.o ReadTabix.o SingleAllele.o mistarOperations.o tabix/libtabix.a libgab/utils.o  libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

libnuc.o:	libnuc.cpp
	${CXX} ${CXXFLAGS} libnuc.cpp


nucular:	nucular.o  tabix/libtabix.a libgab/utils.o  libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

nucularMCMC:	nucularMCMC.o tabix/libtabix.a ${BAMTOOLS}/lib/libbamtools.a libgab/utils.o  libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f *.o nucular  nucularMCMC  nucularMCMCBAM

