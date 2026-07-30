target_sources(engine PRIVATE
        BarnesHutTree.cu
        BarnesHutKernel.cu
)

set_source_files_properties(
        BarnesHutTree.cu
        BarnesHutKernel.cu
        PROPERTIES
        LANGUAGE CUDA
)