add_library(nbody_engine)

target_include_directories(nbody_engine
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

add_subdirectory(particles)
add_subdirectory(calculators)
add_subdirectory(distributions)
add_subdirectory(integrators)
add_subdirectory(parameters)
add_subdirectory(math)

if(MSVC)
    target_compile_options(nbody_engine PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:/W4>
    )
else()
    target_compile_options(nbody_engine PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic>
    )
endif()