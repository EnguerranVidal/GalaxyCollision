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

add_subdirectory(cuda)