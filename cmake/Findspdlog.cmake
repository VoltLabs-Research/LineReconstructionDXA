if(TARGET spdlog::spdlog_header_only)
    set(spdlog_FOUND TRUE)
    return()
endif()

find_package(Threads REQUIRED)
find_package(fmt REQUIRED)

set(_spdlog_include_dirs
    "/usr/include"
    "/usr/local/include"
)

find_path(SPDLOG_INCLUDE_DIR spdlog/spdlog.h PATHS ${_spdlog_include_dirs})

if(NOT SPDLOG_INCLUDE_DIR)
    message(FATAL_ERROR "spdlog headers not found")
endif()

add_library(spdlog::spdlog_header_only INTERFACE IMPORTED)
set_target_properties(spdlog::spdlog_header_only PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SPDLOG_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "Threads::Threads;fmt::fmt-header-only"
    INTERFACE_COMPILE_DEFINITIONS "SPDLOG_FMT_EXTERNAL;SPDLOG_HEADER_ONLY;SPDLOG_FWRITE_UNLOCKED"
)

if(NOT TARGET spdlog::spdlog)
    add_library(spdlog::spdlog INTERFACE IMPORTED)
    set_target_properties(spdlog::spdlog PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${SPDLOG_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "spdlog::spdlog_header_only"
        INTERFACE_COMPILE_DEFINITIONS "SPDLOG_FMT_EXTERNAL;SPDLOG_HEADER_ONLY;SPDLOG_FWRITE_UNLOCKED"
    )
endif()

set(spdlog_FOUND TRUE)
