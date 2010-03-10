# Make sure the compiler can find include files from our cminpack library.
include_directories (${CMINPACK_SOURCE_DIR})

# Make sure the linker can find the cminpack library once it is built.
link_directories (${CMINPACK_BINARY_DIR}/Hello)

SET(FPGM tchkder thybrd thybrd1 thybrj thybrj1 tlmder tlmder1 tlmdif 	tlmdif1 tlmstr tlmstr1)
SET(PGM tchkder_ thybrd_ thybrd1_ thybrj_ thybrj1_ tlmder_ tlmder1_ tlmdif_ tlmdif1_ tlmstr_ tlmstr1_ tfdjac2_)
SET(CPGM tchkderc thybrdc thybrd1c thybrjc thybrj1c tlmderc tlmder1c tlmdifc tlmdif1c tlmstrc tlmstr1c tfdjac2c)

#enable_language(Fortran)
#FOREACH(source ${FPGM})
#  ADD_EXECUTABLE(${source} ${source}.f)
#  TARGET_LINK_LIBRARIES (${source} cminpack) 
#ENDFOREACH(source)

FOREACH(source ${PGM})
  ADD_EXECUTABLE(${source} ${source}.c)
  TARGET_LINK_LIBRARIES (${source} cminpack) 
ENDFOREACH(source)

FOREACH(source ${CPGM})
  ADD_EXECUTABLE(${source} ${source}.c)
  TARGET_LINK_LIBRARIES (${source} cminpack) 
ENDFOREACH(source)


