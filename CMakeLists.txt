#===================================================================================================== PROJECT SETUP ===
cmake_minimum_required(VERSION 3.19)
cmake_policy(SET CMP0110 NEW)
project(tardigrade_cmml)

string(TOUPPER "${PROJECT_NAME}" project_name_upper)

# Set common project paths relative to project root directory
set(CPP_SRC_PATH "src/cpp")
set(CPP_TEST_PATH "${CPP_SRC_PATH}/tests")
set(CMAKE_SRC_PATH "src/cmake")
set(DOXYGEN_SRC_PATH "docs/doxygen")
set(SPHINX_SRC_PATH "docs/sphinx")

# Add a flag for whether the python bindings should be built or not
set(TARDIGRADE_CMML_BUILD_PYTHON_BINDINGS
    ON
    CACHE BOOL
    "Flag for whether the python bindings should be built for hydra"
)
set(TARDIGRADE_CMML_USE_LIBXSMM OFF CACHE BOOL "Flag for whether to use libxsmm for matrix math")
set("${project_name_upper}_BUILD_AS_SHARED"
    OFF
    CACHE BOOL
    "Flag for whether to build the library as a shared object library or header only"
)

if(${project_name_upper}_BUILD_AS_SHARED)
    add_compile_definitions(TARDIGRADE_CMML_BUILD_AS_SHARED)
endif()

# Add a flag for if a full build of all tardigrade repositories should be performed
set(TARDIGRADE_FULL_BUILD
    OFF
    CACHE BOOL
    "Flag for whether a full build of Tardigrade should be performed (i.e., all repos pulled from git and built)"
)

# Set the material library names
set(HEADER_ONLY_MATERIAL_MODELS "tardigrade_BasicSolid" "tardigrade_DefinedPlasticEvolution")
set(HEADER_ONLY_MATERIAL_MODEL_UPSTREAM_REQUIREMENTS "tardigrade_hydra")

# Set the internal support libraries
set(INTERNAL_SUPPORT_LIBRARIES)
set(ADDITIONAL_HEADER_ONLY_LIBRARIES "tardigrade_CMMLMultiphaseMaterialBase")
set(ADDITIONAL_HEADER_ONLY_UPSTREAM_LIBRARIES)

set(ADDITIONAL_HEADER_ONLY_LIBRARIES ${ADDITIONAL_HEADER_ONLY_LIBRARIES} ${HEADER_ONLY_MATERIAL_MODELS})
set(ADDITIONAL_HEADER_ONLY_UPSTREAM_LIBRARIES
    ${ADDITIONAL_HEADER_ONLY_UPSTREAM_LIBRARIES}
    ${HEADER_ONLY_MATERIAL_MODEL_UPSTREAM_REQUIREMENTS}
)
set(PROJECT_SOURCE_FILES ${PROJECT_NAME}.cpp ${PROJECT_NAME}.h ${PROJECT_NAME}_header_only_materials.h)
# For each support library you need to set the upstream libraries required via
#set(TARDIGRADE_CMMLMYSUPPORTLIBRARY_UPSTREAM_LIBRARIES "tardigrade_constitutive_tools")

# Set some Abaqus integration variables
set(MISSING_INTEGRATION_DEPS_MESSAGE "Abaqus integration tests will be skipped.")

# Add the cmake folder to locate project CMake module(s)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}" ${CMAKE_MODULE_PATH})

# Set build type checks
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lower)
set(upstream_required "")
set(not_conda_test "true")
if(cmake_build_type_lower STREQUAL "release")
    set(upstream_required "REQUIRED")
elseif(cmake_build_type_lower STREQUAL "conda-test")
    set(upstream_required "REQUIRED")
    set(not_conda_test "false")
    # Find the installed project umat
endif()

set(PROJECT_PRIVATE_HEADERS "")
foreach(package ${INTERNAL_SUPPORT_LIBRARIES})
    set(PROJECT_SOURCE_FILES ${PROJECT_SOURCE_FILES} ${package}.cpp ${package}.h)
    set(PROJECT_PRIVATE_HEADERS ${PROJECT_PRIVATE_HEADERS} ${package}.h)
endforeach(package)
if(${not_conda_test} STREQUAL "true")
    foreach(package ${ADDITIONAL_HEADER_ONLY_LIBRARIES})
        add_library(${package} INTERFACE "${CPP_SRC_PATH}/${package}.cpp" "${CPP_SRC_PATH}/${package}.h")
    endforeach(package)
endif()
if(TARDIGRADE_CMML_USE_LIBXSMM)
    message(STATUS "Using libxsmm for matrix math")
    find_library(libxsmm_DIR xsmm REQUIRED PATHS "$ENV{CONDA_PREFIX}/lib")
    add_definitions(-DTARDIGRADE_CMML_USE_LLXSMM)
    set(PROJECT_LINK_LIBRARIES ${PROJECT_LINK_LIBRARIES} "xsmm" "blas")
endif()

if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

# Get version number from Git
if(${not_conda_test} STREQUAL "true")
    if(DEFINED ENV{PREFIX})
        set(Python ROOT_DIR "$ENV{PREFIX}/bin")
    endif()
    set(Python_FIND_STRATEGY LOCATION)
    find_package(Python COMPONENTS Interpreter REQUIRED)
    execute_process(
        COMMAND ${Python_EXECUTABLE} -m setuptools_scm
        OUTPUT_VARIABLE ${PROJECT_NAME}_VERSION_STRING_FULL
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(${${PROJECT_NAME}_VERSION_STRING_FULL} STREQUAL "")
        set(${PROJECT_NAME}_VERSION 0.0.0)
    else()
        string(REGEX MATCH "[0-9]+\.[0-9]+\.[0-9]+" ${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_VERSION_STRING_FULL})
    endif()
    project(${PROJECT_NAME} VERSION ${${PROJECT_NAME}_VERSION})
endif()

# Add installation directory variables
include(GNUInstallDirs)

# Make the code position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set the c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt -O3 -DNDEBUG -ggdb")
endif()

# Enable CTest
enable_testing()

#================================================================================================= FIND DEPENDENCIES ===

#=============================================================================================== FIND PROGRAMS ===

find_program(CLANG_FORMAT_EXECUTABLE NAMES clang-format REQUIRED)
set_property(GLOBAL PROPERTY CLANG_FORMAT_SOURCE_FILES)

#=============================================================================================== ADD PROJECT TARGETS ===

# Find eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
if(EIGEN3_FOUND)
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
endif()

# Find bash (Required for abaqus integration tests)
find_program(BASH_PROGRAM bash)
if(BASH_PROGRAM)
    message(STATUS "Found bash: ${BASH_PROGRAM}")
else()
    message(STATUS "bash not found. ${MISSING_INTEGRATION_DEPS_MESSAGE}")
endif()

# Find related, but external, projects in installed environment
include(FetchContent)
set(upstream_packages "tardigrade_error_tools" "tardigrade_vector_tools" "tardigrade_hydra")
if(${not_conda_test} STREQUAL "false")
    # During conda-build testing, we must find the installed project files as if they were an external project
    set(upstream_packages ${upstream_packages} ${PROJECT_NAME})
endif()
set(LOCAL_BUILD_INCLUDE_DIRECTORIES "")
foreach(package ${upstream_packages})
    string(TOUPPER "${package}" package_upper)
    set(${package_upper}_BUILD_PYTHON_BINDINGS
        ${TARDIGRADE_CMML_BUILD_PYTHON_BINDINGS}
        CACHE INTERNAL
        "Setting ${package}'s python binding flag to the global value"
    )
    if(TARDIGRADE_FULL_BUILD)
        message("Not attempting to find ${package}")
    else()
        find_package(${package} ${upstream_required} CONFIG)
    endif()
    if(${package}_FOUND)
        message(STATUS "Found ${package}: ${${package}_DIR}")
    else()
        # Find related, but external, projects using FetchContent and building locally
        message(WARNING "Did not find an installed ${package} package. Attempting local build with FetchContent.")
        if(NOT DEFINED ${package_upper}_FETCHCONTENT_VERSION)
            set(${package_upper}_FETCHCONTENT_VERSION "origin/dev")
        endif()
        message("${package_upper} is being built with version ${${package_upper}_FETCHCONTENT_VERSION}")
        FetchContent_Declare(
            ${package}
            GIT_REPOSITORY https://github.com/UCBoulder/${package}.git
            GIT_TAG ${${package_upper}_FETCHCONTENT_VERSION}
        )
        FetchContent_MakeAvailable(${package})
        set(LOCAL_BUILD_INCLUDE_DIRECTORIES ${LOCAL_BUILD_INCLUDE_DIRECTORIES} ${${package}_SOURCE_DIR}/${CPP_SRC_PATH})
    endif()
endforeach(package)

foreach(package ${upstream_packages})
    if(NOT cmake_build_type_lower STREQUAL "release")
        string(TOLOWER "${package}" lower_package)
        foreach(sub_package ${${upper_package}_UPSTREAM_LIBRARIES})
            target_include_directories(${package} PUBLIC ${${sub_package}_SOURCE_DIR}/${CPP_SRC_PATH})
        endforeach(sub_package)
    endif()
endforeach(package)

# Only add tests and documentation for current project builds. Protects downstream project builds.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Add test-only external projects
    set(test_upstream_packages)

    foreach(package ${test_upstream_packages})
        string(TOUPPER "${package}" package_upper)
        set(${package_upper}_BUILD_PYTHON_BINDINGS
            ${TARDIGRADE_CMML_BUILD_PYTHON_BINDINGS}
            CACHE INTERNAL
            "Setting ${package}'s python binding flag to the global value"
        )
        if(TARDIGRADE_FULL_BUILD)
            message("Not attempting to find ${package}")
        else()
            find_package(${package} ${upstream_required} CONFIG)
        endif()
        if(${package}_FOUND)
            message(STATUS "Found ${package}: ${${package}_DIR}")
        else()
            # Find related, but external, projects using FetchContent and building locally
            message(WARNING "Did not find an installed ${package} package. Attempting local build with FetchContent.")
            if(NOT DEFINED ${package_upper}_FETCHCONTENT_VERSION)
                set(${package_upper}_FETCHCONTENT_VERSION "origin/dev")
            endif()
            message("${package_upper} is being built with version ${${package_upper}_FETCHCONTENT_VERSION}")
            FetchContent_Declare(
                ${package}
                GIT_REPOSITORY https://github.com/UCBoulder/${package}.git
                GIT_TAG ${${package_upper}_FETCHCONTENT_VERSION}
            )
            FetchContent_MakeAvailable(${package})
            set(LOCAL_BUILD_INCLUDE_DIRECTORIES
                ${LOCAL_BUILD_INCLUDE_DIRECTORIES}
                ${${package}_SOURCE_DIR}/${CPP_SRC_PATH}
            )
        endif()
    endforeach(package)
endif()

#=============================================================================================== ADD PROJECT TARGETS ===
# MUST COME AFTER DEPENDENCY LOCATING
# Add project source directories
if(${not_conda_test} STREQUAL "true")
    include_directories(${CPP_SRC_PATH})
    add_subdirectory(${CPP_SRC_PATH})
endif()

# Only add tests and documentation for current project builds. Protects downstream project builds.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Find Boost. Required for tests
    find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)
    # Add c++ tests and docs
    add_subdirectory(${CPP_TEST_PATH})
    if(${not_conda_test} STREQUAL "true")
        add_subdirectory(${DOXYGEN_SRC_PATH})
        add_subdirectory(${SPHINX_SRC_PATH})
    endif()
endif()

#==================================================================================== SETUP INSTALLATION CMAKE FILES ===
if(${not_conda_test} STREQUAL "true")
    foreach(package ${PROJECT_NAME})
        include(CMakePackageConfigHelpers)
        write_basic_package_version_file(
            "${package}ConfigVersion.cmake"
            VERSION ${PROJECT_VERSION}
            COMPATIBILITY SameMajorVersion
        )
        configure_package_config_file(
            "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}/Config.cmake.in"
            "${PROJECT_BINARY_DIR}/${package}Config.cmake"
            INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake
        )

        # CMake won't build the targets for local builds of upstream projects
        if(cmake_build_type_lower STREQUAL release)
            install(
                EXPORT ${package}_Targets
                FILE ${package}Targets.cmake
                DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake
            )
        endif()

        install(
            FILES "${PROJECT_BINARY_DIR}/${package}Config.cmake" "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
            DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake
        )
    endforeach(package)
endif()

#================================================================================== SETUP CLANG FORMATTING TARGETS ===

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    get_property(CLANG_FORMAT_SOURCE_FILES_VALUE GLOBAL PROPERTY CLANG_FORMAT_SOURCE_FILES)

    add_custom_target(
        cpp-format-check
        COMMAND ${CLANG_FORMAT_EXECUTABLE} -style=file --dry-run -Werror -i ${CLANG_FORMAT_SOURCE_FILES_VALUE}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Check the C++ code format with clang-format"
    )
    set_target_properties(cpp-format-check PROPERTIES EXCLUDE_FROM_ALL TRUE)

    add_custom_target(
        cpp-format
        COMMAND ${CLANG_FORMAT_EXECUTABLE} -style=file -i ${CLANG_FORMAT_SOURCE_FILES_VALUE}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Format the C++ code with clang-format"
    )
    set_target_properties(cpp-format PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()
