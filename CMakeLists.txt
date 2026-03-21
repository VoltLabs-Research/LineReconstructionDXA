#
# Copyright (C) Rodolfo Herrera Hernandez. All rights reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(line-reconstruction-dxa VERSION 1.0.0 LANGUAGES CXX)

option(LINE_RECONSTRUCTION_DXA_BOOTSTRAP_LOCAL_DEPS "Bootstrap local VoltLabs dependencies when packaged configs are unavailable" ON)
set(VOLTLABS_ROOT "${CMAKE_SOURCE_DIR}/../voltlabs" CACHE PATH "Path to the local VoltLabs workspace")
set(CORETOOLKIT_SOURCE_DIR "${VOLTLABS_ROOT}/CoreToolkit" CACHE PATH "Local CoreToolkit source directory")
set(STRUCTURE_IDENTIFICATION_SOURCE_DIR "${VOLTLABS_ROOT}/StructureIdentification" CACHE PATH "Local StructureIdentification source directory")
set(CORETOOLKIT_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_deps/install/coretoolkit" CACHE PATH "Install prefix for local CoreToolkit bootstrap")
set(STRUCTURE_IDENTIFICATION_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_deps/install/structure-identification" CACHE PATH "Install prefix for local StructureIdentification bootstrap")

if(LINE_RECONSTRUCTION_DXA_BOOTSTRAP_LOCAL_DEPS)
    list(PREPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)

find_package(coretoolkit REQUIRED)
find_package(structure-identification REQUIRED)
find_package(common-neighbor-analysis REQUIRED)
find_package(polyhedral-template-matching REQUIRED)
find_package(opendxa REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)

get_target_property(OPENDXA_INCLUDE_DIRS opendxa::opendxa INTERFACE_INCLUDE_DIRECTORIES)
if(NOT OPENDXA_INCLUDE_DIRS OR OPENDXA_INCLUDE_DIRS STREQUAL "OPENDXA_INCLUDE_DIRS-NOTFOUND")
    set(OPENDXA_INCLUDE_DIRS)
endif()
set(LOCAL_OPENDXA_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../OpenDXA/include")
if(EXISTS "${LOCAL_OPENDXA_INCLUDE_DIR}")
    list(PREPEND OPENDXA_INCLUDE_DIRS "${LOCAL_OPENDXA_INCLUDE_DIR}")
endif()

function(append_imported_library out_var target_name)
    if(NOT TARGET "${target_name}")
        return()
    endif()

    get_target_property(imported_location "${target_name}" IMPORTED_LOCATION_RELEASE)
    if(NOT imported_location OR imported_location STREQUAL "imported_location-NOTFOUND")
        get_target_property(imported_location "${target_name}" IMPORTED_LOCATION)
    endif()

    if(imported_location AND NOT imported_location STREQUAL "imported_location-NOTFOUND")
        set(${out_var} ${${out_var}} "${imported_location}" PARENT_SCOPE)
    endif()
endfunction()

file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/*.cpp)
list(FILTER LIB_SOURCES EXCLUDE REGEX ".*main\.cpp$")

add_library(${PROJECT_NAME}_lib STATIC ${LIB_SOURCES})
set_target_properties(${PROJECT_NAME}_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(${PROJECT_NAME}_lib BEFORE PUBLIC
    ${OPENDXA_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/include
)
target_precompile_headers(${PROJECT_NAME}_lib PRIVATE <volt/core/volt.h>)
target_link_libraries(${PROJECT_NAME}_lib PRIVATE
    opendxa::opendxa
    coretoolkit::coretoolkit
    structure-identification::structure-identification
    common-neighbor-analysis::common-neighbor-analysis
    polyhedral-template-matching::polyhedral-template-matching
    spdlog::spdlog
    nlohmann_json::nlohmann_json
)

add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_precompile_headers(${PROJECT_NAME} PRIVATE <volt/core/volt.h>)

set(LRDXA_MANUAL_LINK_LIBS)
append_imported_library(LRDXA_MANUAL_LINK_LIBS "CONAN_LIB::opendxa_opendxa_lib_RELEASE")
append_imported_library(LRDXA_MANUAL_LINK_LIBS "CONAN_LIB::structure-identification_structure-identification_RELEASE")
append_imported_library(LRDXA_MANUAL_LINK_LIBS "CONAN_LIB::common-neighbor-analysis_common-neighbor-analysis_lib_RELEASE")
append_imported_library(LRDXA_MANUAL_LINK_LIBS "CONAN_LIB::polyhedral-template-matching_polyhedral-template-matching_lib_RELEASE")
append_imported_library(LRDXA_MANUAL_LINK_LIBS "CONAN_LIB::coretoolkit_coretoolkit_RELEASE")
append_imported_library(LRDXA_MANUAL_LINK_LIBS "CONAN_LIB::coretoolkit_ptm_RELEASE")
append_imported_library(LRDXA_MANUAL_LINK_LIBS "CONAN_LIB::coretoolkit_mwm_csp_RELEASE")
append_imported_library(LRDXA_MANUAL_LINK_LIBS "CONAN_LIB::coretoolkit_geogram_RELEASE")

target_link_libraries(${PROJECT_NAME} PRIVATE
    ${PROJECT_NAME}_lib
    opendxa::opendxa
    structure-identification::structure-identification
    common-neighbor-analysis::common-neighbor-analysis
    polyhedral-template-matching::polyhedral-template-matching
    coretoolkit::coretoolkit
    spdlog::spdlog
    nlohmann_json::nlohmann_json
    ${LRDXA_MANUAL_LINK_LIBS}
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
install(TARGETS ${PROJECT_NAME}_lib DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)
