#!/bin/bash

PROG=$(readlink -f $0)
GOATDIR=`echo $(dirname "${PROG}") | sed 's/\/scripts//'`

NAME=$1
MAIN=`echo $NAME | sed 's/-/_/'`
CODE=`echo $NAME | sed 's/-/ /g' | sed -e "s/\b\(.\)/\u\1/g" | sed 's/ //g'`

# Create new CMakeLists.txt

LINE=`grep -n pi0-example ${GOATDIR}/CMakeLists.txt | tail -1 | sed 's/:/ /' | awk '{print $1+1}'`
head -n $LINE ${GOATDIR}/CMakeLists.txt > ${GOATDIR}/CMakeLists.new

echo "add_executable(${NAME}" >> ${GOATDIR}/CMakeLists.new
echo "   \${GOAT_BASE}" >> ${GOATDIR}/CMakeLists.new
echo "   \${GOAT_GHIST}" >> ${GOATDIR}/CMakeLists.new
echo "   inc/PPhysics.h" >> ${GOATDIR}/CMakeLists.new
echo "   src/PPhysics.cc" >> ${GOATDIR}/CMakeLists.new
echo "   inc/P${CODE}.h" >> ${GOATDIR}/CMakeLists.new
echo "   src/P${CODE}.cc" >> ${GOATDIR}/CMakeLists.new
echo "   src/${MAIN}_main.cc" >> ${GOATDIR}/CMakeLists.new
echo "   )" >> ${GOATDIR}/CMakeLists.new
echo "target_link_libraries(${NAME} \${LIBS} \${ROOT_LIBRARIES})" >> ${GOATDIR}/CMakeLists.new

tail -n +$LINE ${GOATDIR}/CMakeLists.txt >> ${GOATDIR}/CMakeLists.new

mv ${GOATDIR}/CMakeLists.new ${GOATDIR}/CMakeLists.txt

# Create new main

cat ${GOATDIR}/src/pi0_example_main.cc | sed "s/Pi0Example/${CODE}/g" > ${GOATDIR}/src/${MAIN}_main.cc

# Create new header

cat ${GOATDIR}/inc/PPi0Example.h | sed "s/Pi0Example/${CODE}/g" > ${GOATDIR}/inc/P${CODE}.h

# Create new source

cat ${GOATDIR}/src/PPi0Example.cc | sed "s/Pi0Example/${CODE}/g" > ${GOATDIR}/src/P${CODE}.cc

echo "Finished. Now just edit:"
echo "     ${GOATDIR}/inc/P${CODE}.h"
echo "     ${GOATDIR}/src/P${CODE}.cc"
echo "Then recompile and run ${NAME}"
