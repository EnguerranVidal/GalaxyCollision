target_sources(nbody_engine PRIVATE
        Calculator.cpp
        Calculator.h
        NewtonCalculator.cpp
        NewtonCalculator.h
        BarnesHutCalculator.cpp
        BarnesHutCalculator.h
        BarnesHutTree.cpp
        BarnesHutTree.h
)

if(GALAXY_ENABLE_CUDA)
    add_subdirectory(cuda)
endif()