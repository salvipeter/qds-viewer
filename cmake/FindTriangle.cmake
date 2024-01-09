if (NOT Triangle_FOUND)
  find_path(TRIANGLE_DIR NAMES triangle.h triangle.o)
endif()

include(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(Triangle
  REQUIRED_VARS TRIANGLE_DIR
  FOUND_VAR Triangle_FOUND
)
