set(dir ${CMAKE_CURRENT_SOURCE_DIR})

list(
    APPEND
    srcs
    ${dir}/filter.cpp
    ${dir}/FiniteElement.cpp
    ${dir}/mma.cpp
    ${dir}/Topopt.cpp
)

set(srcs ${srcs} PARENT_SCOPE)