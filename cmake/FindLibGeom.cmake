if (NOT LibGeom_FOUND)
    find_path(LIBGEOM_INCLUDE_DIR NAMES geometry.hh)
    find_library(LIBGEOM_LIBRARY_RELEASE
        NAMES geom
        HINTS ${LIBGEOM_INCLUDE_DIR}
        PATH_SUFFIXES build release
        )
    find_library(LIBGEOM_LIBRARY_DEBUG
        NAMES geom
        HINTS ${LIBGEOM_INCLUDE_DIR}
        PATH_SUFFIXES build debug
        )
    if (LIBGEOM_LIBRARY_RELEASE)
        if (LIBGEOM_LIBRARY_DEBUG)
            set(LIBGEOM_LIBRARY optimized ${LIBGEOM_LIBRARY_RELEASE} debug ${LIBGEOM_LIBRARY_DEBUG})
        else()
            set(LIBGEOM_LIBRARY ${LIBGEOM_LIBRARY_RELEASE})
        endif()
    endif()
endif()

include(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(LibGeom
    REQUIRED_VARS LIBGEOM_INCLUDE_DIR LIBGEOM_LIBRARY
    FOUND_VAR LibGeom_FOUND
    )
