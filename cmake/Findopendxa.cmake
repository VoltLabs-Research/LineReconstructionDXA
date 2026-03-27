if(TARGET opendxa::opendxa)
    set(opendxa_FOUND TRUE)
    return()
endif()

find_package(coretoolkit REQUIRED)
find_package(structure-identification REQUIRED)
find_package(Threads REQUIRED)
find_package(TBB REQUIRED)
find_package(Boost REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)

if(NOT DEFINED OPENDXA_SOURCE_DIR OR NOT EXISTS "${OPENDXA_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "opendxa not found. Set OPENDXA_SOURCE_DIR to a valid source tree.")
endif()

if(NOT DEFINED OPENDXA_INSTALL_PREFIX OR OPENDXA_INSTALL_PREFIX STREQUAL "")
    set(OPENDXA_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_deps/install/opendxa")
endif()

set(_opendxa_library "${OPENDXA_INSTALL_PREFIX}/lib/libopendxa_lib.a")
set(_opendxa_build_dir "${CMAKE_BINARY_DIR}/_deps/opendxa-build")

if(NOT EXISTS "${_opendxa_library}")
    execute_process(
        COMMAND cmake -S "${OPENDXA_SOURCE_DIR}" -B "${_opendxa_build_dir}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${OPENDXA_INSTALL_PREFIX} -DCMAKE_MODULE_PATH=${CMAKE_SOURCE_DIR}/cmake -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=FALSE -DCORETOOLKIT_SOURCE_DIR=${CORETOOLKIT_SOURCE_DIR} -DCORETOOLKIT_INSTALL_PREFIX=${CORETOOLKIT_INSTALL_PREFIX} -DSTRUCTURE_IDENTIFICATION_SOURCE_DIR=${STRUCTURE_IDENTIFICATION_SOURCE_DIR} -DSTRUCTURE_IDENTIFICATION_INSTALL_PREFIX=${STRUCTURE_IDENTIFICATION_INSTALL_PREFIX}
        RESULT_VARIABLE _opendxa_configure_result
    )
    if(NOT _opendxa_configure_result EQUAL 0)
        message(FATAL_ERROR "Failed to configure OpenDXA from ${OPENDXA_SOURCE_DIR}")
    endif()

    execute_process(
        COMMAND cmake --build "${_opendxa_build_dir}"
        RESULT_VARIABLE _opendxa_build_result
    )
    if(NOT _opendxa_build_result EQUAL 0)
        message(FATAL_ERROR "Failed to build OpenDXA from ${OPENDXA_SOURCE_DIR}")
    endif()

    execute_process(
        COMMAND cmake --install "${_opendxa_build_dir}"
        RESULT_VARIABLE _opendxa_install_result
    )
    if(NOT _opendxa_install_result EQUAL 0)
        message(FATAL_ERROR "Failed to install OpenDXA into ${OPENDXA_INSTALL_PREFIX}")
    endif()
endif()

set(_opendxa_boost_target "")
if(TARGET Boost::headers)
    set(_opendxa_boost_target Boost::headers)
elseif(TARGET Boost::boost)
    set(_opendxa_boost_target Boost::boost)
elseif(TARGET boost::headers)
    set(_opendxa_boost_target boost::headers)
elseif(TARGET boost::boost)
    set(_opendxa_boost_target boost::boost)
else()
    message(FATAL_ERROR "OpenDXA requires a Boost headers target")
endif()

add_library(opendxa STATIC IMPORTED GLOBAL)
set_target_properties(opendxa PROPERTIES
    IMPORTED_LOCATION "${_opendxa_library}"
    INTERFACE_INCLUDE_DIRECTORIES "${OPENDXA_INSTALL_PREFIX}/include"
    INTERFACE_LINK_LIBRARIES "${_opendxa_boost_target};Threads::Threads;TBB::tbb;structure-identification::structure-identification;coretoolkit::coretoolkit;spdlog::spdlog;nlohmann_json::nlohmann_json"
)

add_library(opendxa::opendxa ALIAS opendxa)

set(opendxa_FOUND TRUE)
