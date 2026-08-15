target_sources(nbody_engine PRIVATE
        NewtonKernel.cu
        NewtonKernel.cuh
        BarnesHutKernel.cu
        BarnesHutKernel.cuh
)

set_source_files_properties(
        NewtonKernel.cu
        NewtonKernel.cuh
        BarnesHutKernel.cu
        BarnesHutKernel.cuh
        PROPERTIES
        LANGUAGE CUDA
)

target_compile_definitions(nbody_engine PRIVATE GALAXY_ENABLE_CUDA=1)
set_property(TARGET nbody_engine PROPERTY CUDA_ARCHITECTURES 75 80 86 89)