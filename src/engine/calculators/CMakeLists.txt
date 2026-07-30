target_sources(engine PRIVATE
        Calculator.cpp
        Calculator.h
        NewtonCalculator.cpp
        NewtonCalculator.h
        BarnesHutCalculator.cpp
        BarnesHutCalculator.h
)

add_subdirectory(cuda)