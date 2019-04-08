CXX      = g++ -g 



BAMTOOLS= $(realpath ../bamtools/)
LIBGAB= $(realpath ../libgab/)

#LIBTABIX = $(realpath tabix/)

CXXFLAGS = -Wall -lm -O3 -lz -std=c++11 -I${LIBGAB} -I${LIBGAB}/gzstream/  -I${BAMTOOLS}/src/ -Itabix -I${BAMTOOLS}/src/  -c
LDFLAGS  = tabix/libtabix.a  ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o  ${BAMTOOLS}/build/src/api/libbamtools.a ${LIBGAB}/utils.o  -lpthread -lm -lz
LDFLAGSNBT  = tabix/libtabix.a ${LIBGAB}/utils.o  -lpthread -lm -lz



all: dice diceBAM BAM2DICE


tabix/tabix.h:
	rm -rf tabix/*
	git clone --recursive https://github.com/samtools/tabix.git

tabix/libtabix.a: tabix/tabix.h
	make -C tabix

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

dice:	dice.o tabix/libtabix.a  ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGSNBT)

clean :
	rm -f *.o dice  diceBAM BAM2DICE


allelefreq :
	mkdir -p ../alleleFreqNuc
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/ACB.mst.gz     -O ../alleleFreqNuc/ACB.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/ACB.mst.gz.tbi -O ../alleleFreqNuc/ACB.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/ASW.mst.gz     -O ../alleleFreqNuc/ASW.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/ASW.mst.gz.tbi -O ../alleleFreqNuc/ASW.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/BEB.mst.gz     -O ../alleleFreqNuc/BEB.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/BEB.mst.gz.tbi -O ../alleleFreqNuc/BEB.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/CDX.mst.gz     -O ../alleleFreqNuc/CDX.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/CDX.mst.gz.tbi -O ../alleleFreqNuc/CDX.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/CEU.mst.gz     -O ../alleleFreqNuc/CEU.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/CEU.mst.gz.tbi -O ../alleleFreqNuc/CEU.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/CHB.mst.gz     -O ../alleleFreqNuc/CHB.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/CHB.mst.gz.tbi -O ../alleleFreqNuc/CHB.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/CHS.mst.gz     -O ../alleleFreqNuc/CHS.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/CHS.mst.gz.tbi -O ../alleleFreqNuc/CHS.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/CLM.mst.gz     -O ../alleleFreqNuc/CLM.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/CLM.mst.gz.tbi -O ../alleleFreqNuc/CLM.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/ESN.mst.gz     -O ../alleleFreqNuc/ESN.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/ESN.mst.gz.tbi -O ../alleleFreqNuc/ESN.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/FIN.mst.gz     -O ../alleleFreqNuc/FIN.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/FIN.mst.gz.tbi -O ../alleleFreqNuc/FIN.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/GBR.mst.gz     -O ../alleleFreqNuc/GBR.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/GBR.mst.gz.tbi -O ../alleleFreqNuc/GBR.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/GIH.mst.gz     -O ../alleleFreqNuc/GIH.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/GIH.mst.gz.tbi -O ../alleleFreqNuc/GIH.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/GWD.mst.gz     -O ../alleleFreqNuc/GWD.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/GWD.mst.gz.tbi -O ../alleleFreqNuc/GWD.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/IBS.mst.gz     -O ../alleleFreqNuc/IBS.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/IBS.mst.gz.tbi -O ../alleleFreqNuc/IBS.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/ITU.mst.gz     -O ../alleleFreqNuc/ITU.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/ITU.mst.gz.tbi -O ../alleleFreqNuc/ITU.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/JPT.mst.gz     -O ../alleleFreqNuc/JPT.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/JPT.mst.gz.tbi -O ../alleleFreqNuc/JPT.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/KHV.mst.gz     -O ../alleleFreqNuc/KHV.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/KHV.mst.gz.tbi -O ../alleleFreqNuc/KHV.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/LWK.mst.gz     -O ../alleleFreqNuc/LWK.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/LWK.mst.gz.tbi -O ../alleleFreqNuc/LWK.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/MSL.mst.gz     -O ../alleleFreqNuc/MSL.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/MSL.mst.gz.tbi -O ../alleleFreqNuc/MSL.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/MXL.mst.gz     -O ../alleleFreqNuc/MXL.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/MXL.mst.gz.tbi -O ../alleleFreqNuc/MXL.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/PEL.mst.gz     -O ../alleleFreqNuc/PEL.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/PEL.mst.gz.tbi -O ../alleleFreqNuc/PEL.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/PJL.mst.gz     -O ../alleleFreqNuc/PJL.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/PJL.mst.gz.tbi -O ../alleleFreqNuc/PJL.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/PUR.mst.gz     -O ../alleleFreqNuc/PUR.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/PUR.mst.gz.tbi -O ../alleleFreqNuc/PUR.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/STU.mst.gz     -O ../alleleFreqNuc/STU.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/STU.mst.gz.tbi -O ../alleleFreqNuc/STU.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/TSI.mst.gz     -O ../alleleFreqNuc/TSI.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/TSI.mst.gz.tbi -O ../alleleFreqNuc/TSI.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/YRI.mst.gz     -O ../alleleFreqNuc/YRI.mst.gz
	wget https://bioinf.eva.mpg.de/dice/alleleFreqNuc/YRI.mst.gz.tbi -O ../alleleFreqNuc/YRI.mst.gz.tbi

testdatasimulated :
	mkdir -p ../testData/simulated/
	wget https://bioinf.eva.mpg.de/dice/testData/simulated/all.bam           -O ../testData/simulated/all.bam
	wget https://bioinf.eva.mpg.de/dice/testData/simulated/all.bam.bai       -O ../testData/simulated/all.bam.bai
	wget https://bioinf.eva.mpg.de/dice/testData/simulated/all.med.regions   -O ../testData/simulated/all.med.regions
	wget https://bioinf.eva.mpg.de/dice/testData/simulated/all.mst.gz    	 -O ../testData/simulated/all.mst.gz
	wget https://bioinf.eva.mpg.de/dice/testData/simulated/all.mst.gz.tbi 	 -O ../testData/simulated/all.mst.gz.tbi
	wget https://bioinf.eva.mpg.de/dice/testData/simulated/all.ref.fa    	 -O ../testData/simulated/all.ref.fa
	wget https://bioinf.eva.mpg.de/dice/testData/simulated/all.ref.fa.fai 	 -O ../testData/simulated/all.ref.fa.fai




