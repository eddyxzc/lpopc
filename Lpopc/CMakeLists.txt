project(liblpopc)
include_directories(src/Common)
include_directories(src/Core)
include_directories(src/Core/spdlog)
include_directories(src/SparseMatrix)
include_directories(../ThirdParty/ipopt/include/coin)
include_directories(../ThirdParty/arma/include)

link_directories(../ThirdParty/ipopt/lib)
link_directories(../ThirdParty/openblas/lib)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SRC_LIST
src/Common/LpDebug.cpp
src/Common/LpOption.cpp
src/Common/LpOptionList.cpp
src/Common/LpReporter.cpp
src/Common/LpUtils.cpp
src/Core/LpANDeriveChecker.cpp
src/Core/LpBoundsChecker.cpp
src/Core/LpDerivDependciesChecker.cpp
src/Core/LpFiniteDifferenceDerive.cpp
src/Core/LpGuessChecker.cpp
src/Core/LpHessian.cpp
src/Core/LpLpopcAlgorithm.cpp
src/Core/LpLpopcApplication.cpp
src/Core/LpMeshRefiner.cpp
src/Core/LpNLPSolver.cpp
src/Core/LpNLPWrapper.cpp
src/Core/LpopcIpopt.cpp
src/Core/LpOptimalProblem.cpp
src/Core/LpPhMeshRefineAlg.cpp
src/Core/LpSacleOCP.cpp
src/Core/LpSizeChecker.cpp
src/Core/Nlp2OPConverter.cpp
src/Core/RPMGenerator.cpp
src/Core/LpSolutionError.cpp
        src/Core/LpLiuHpMeshRefineAlg.cpp
src/SparseMatrix/LpSparseArray.cpp
src/SparseMatrix/LpSparseMatrix.cpp
)
add_library(lpopc ${SRC_LIST})
target_link_libraries(lpopc ipopt coinmetis coinmumps openblas)