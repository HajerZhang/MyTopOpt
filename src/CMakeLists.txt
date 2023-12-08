set(dir ${CMAKE_CURRENT_SOURCE_DIR})

list(
    APPEND
    srcs
    ${dir}/filter.cpp
    ${dir}/FiniteElement.cpp
    ${dir}/mma.cpp
)

set(srcs ${srcs} PARENT_SCOPE)