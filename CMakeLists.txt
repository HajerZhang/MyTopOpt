cmake_minimum_required(VERSION 3.15)

project(
    myTopOpt
    LANGUAGES CXX
    VERSION 2.0.1
)
set(CMAKE_CXX_COMPILER "clang++")

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 11)
add_definitions("-Wall -Wextra")

include_directories(
    "include"
)
set(srcs)
add_subdirectory("src")

add_library(
    ${PROJECT_NAME}
    ${srcs}
)

add_executable(
    main
    "src/main.cpp"
)

target_link_libraries(
    main
    ${PROJECT_NAME}
)

set_target_properties(
    ${PROJECT_NAME} main
    PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
)

install(
    TARGETS
    main
    "${PROJECT_NAME}"
    EXPORT
    "${PROJECT_NAME}-targets"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
    EXPORT
    "${PROJECT_NAME}-targets"
    NAMESPACE
    "${PROJECT_NAME}::"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
    DIRECTORY
    "${PROJECT_BINARY_DIR}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)