INCLUDE(CheckIncludeFile)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckLibraryExists)

#INCLUDE_DIRECTORIES(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")

ADD_EXECUTABLE(ibd
  Main.cpp Pop.cpp
  rexp.cpp rnormal.cpp
  Main.h Pop.h disperse.h
  xorshift.h
  rexp.h rnormal.h
)

