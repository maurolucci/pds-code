include(FindPackageHandleStandardArgs)
find_path(
  Gurobi_INCLUDE_DIRS
  NAMES gurobi_c.h
  HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
  PATH_SUFFIXES include)

find_library(
  Gurobi_LIBRARY
  NAMES gurobi gurobi81 gurobi90 gurobi95
  HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
  PATH_SUFFIXES lib)

if(MSVC)
  # determine Visual Studio year
  if(MSVC_TOOLSET_VERSION EQUAL 142)
    set(MSVC_YEAR "2019")
  elseif(MSVC_TOOLSET_VERSION EQUAL 141)
    set(MSVC_YEAR "2017")
  elseif(MSVC_TOOLSET_VERSION EQUAL 140)
    set(MSVC_YEAR "2015")
  endif()

  if(MT)
    set(M_FLAG "mt")
  else()
    set(M_FLAG "md")
  endif()

  find_library(
    Gurobi_CXX_LIBRARY
    NAMES gurobi_c++${M_FLAG}${MSVC_YEAR}
    HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
    PATH_SUFFIXES lib)
  find_library(
    Gurobi_CXX_DEBUG_LIBRARY
    NAMES gurobi_c++${M_FLAG}d${MSVC_YEAR}
    HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
    PATH_SUFFIXES lib)
else()
    message("not msvc")
  find_library(
    Gurobi_CXX_LIBRARY
    NAMES gurobi_c++
    HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
    PATH_SUFFIXES lib)
  set(Gurobi_CXX_DEBUG_LIBRARY ${Gurobi_CXX_LIBRARY})
endif()

if (Gurobi_CXX_LIBRARY)
  set(Gurobi_CXX_FOUND YES)
endif()
message("Libraries: c " ${Gurobi_LIBRARY} " cxx " ${Gurobi_CXX_LIBRARY})

find_package_handle_standard_args(Gurobi REQUIRED_VARS Gurobi_LIBRARY
        Gurobi_INCLUDE_DIRS Gurobi_CXX_LIBRARY HANDLE_COMPONENTS)

message(${Gurobi_LIBRARY} " " \"${Gurobi_LIBRARY_FOUND}\")
if (Gurobi_FOUND)
  mark_as_advanced(Gurobi_LIBRARY)
endif()
if (Gurobi_CXX_FOUND)
  mark_as_advanced(Gurobi_CXX_LIBRARY)
endif()


if (Gurobi_FOUND AND NOT TARGET gurobi::gurobi)
  add_library(gurobi::gurobi SHARED IMPORTED)
  set_property(TARGET gurobi::gurobi PROPERTY IMPORTED_LOCATION ${Gurobi_LIBRARY})
  target_include_directories(gurobi::gurobi INTERFACE ${Gurobi_INCLUDE_DIRS})
endif()

if (Gurobi_CXX_FOUND AND NOT TARGET gurobi::gurobi_cxx)
  add_library(gurobi::gurobi_cxx STATIC IMPORTED)
  set_property(TARGET gurobi::gurobi_cxx PROPERTY IMPORTED_LOCATION ${Gurobi_CXX_LIBRARY})
  get_target_property(GUROBI_CXX_IMPORTED_LOC gurobi::gurobi_cxx IMPORTED_LOCATION)
  message("gurobi cxx location " ${GUROBI_CXX_IMPORTED_LOC})
  target_link_libraries(gurobi::gurobi_cxx INTERFACE gurobi::gurobi)
endif()
