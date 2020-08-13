project(lpopc_test)
include_directories(../src/Common)
include_directories(../src/Core)
include_directories(../src/Core/spdlog)
include_directories(../src/SparseMatrix)
include_directories(../../ThirdParty/ipopt/include/coin)
include_directories(../../ThirdParty/arma/include)

include_directories(../example/hypersensitive)

find_library(IPOPT_LIB ipopt ../../ThirdParty/ipopt/lib)
find_library(OPENBLAS_LIB openblas ../../ThirdParty/openblas/lib)
message("Path to Ipopt: ${IPOPT_LIB}")
message("Path to OpenBLAS: ${OPENBLAS_LIB}")
link_directories(../../ThirdParty/openblas/lib)
link_directories(../../ThirdParty/ipopt/lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
#You can replace the source file with your problem or try other examples
set(SOURCE_FILES ../example/hypersensitive/HyperSensitive.cpp)
add_executable(lpopc_test ${SOURCE_FILES})


target_link_libraries(lpopc_test ipopt coinmetis coinmumps lpopc openblas)