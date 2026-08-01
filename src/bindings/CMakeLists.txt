pybind11_add_module(
        engine
        bindings.cpp
        bindings.h
        bindMath.cpp
        bindParticles.cpp
        bindParameters.cpp
        bindCalculator.cpp
        bindIntegrator.cpp
        bindDistribution.cpp
)


target_link_libraries(
        engine
        PRIVATE
        nbody_engine
)