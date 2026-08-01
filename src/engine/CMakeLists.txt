add_library(engine)

target_include_directories(engine
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

add_subdirectory(particles)
add_subdirectory(calculators)
add_subdirectory(distributions)
add_subdirectory(integrators)
add_subdirectory(parameters)
add_subdirectory(utils)
add_subdirectory(math)

if(MSVC)
    target_compile_options(engine PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:/W4>
    )
else()
    target_compile_options(engine PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic>
    )
endif()