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


target_link_libraries(engine PRIVATE nbody_engine)
set_target_properties(engine PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
)
add_custom_command(TARGET engine POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:engine>"
        "${CMAKE_SOURCE_DIR}/"
        COMMENT "Copying engine module next to main.py"
)