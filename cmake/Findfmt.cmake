if(TARGET fmt::fmt-header-only)
    set(fmt_FOUND TRUE)
    return()
endif()

set(_fmt_include_dirs
    "/usr/include"
    "/usr/local/include"
)

find_path(FMT_INCLUDE_DIR fmt/core.h PATHS ${_fmt_include_dirs})

if(NOT FMT_INCLUDE_DIR)
    message(FATAL_ERROR "fmt headers not found")
endif()

add_library(fmt::fmt-header-only INTERFACE IMPORTED)
set_target_properties(fmt::fmt-header-only PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${FMT_INCLUDE_DIR}"
    INTERFACE_COMPILE_DEFINITIONS "FMT_HEADER_ONLY=1"
)

if(NOT TARGET fmt::fmt)
    add_library(fmt::fmt INTERFACE IMPORTED)
    set_target_properties(fmt::fmt PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${FMT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "fmt::fmt-header-only"
        INTERFACE_COMPILE_DEFINITIONS "FMT_HEADER_ONLY=1"
    )
endif()

set(fmt_FOUND TRUE)
