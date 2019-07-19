# - Try to find GUROBI

#  GUROBI_BASE - The libraries needed to use Gurobi

# Once done this will define
#  GUROBI_FOUND - System has Gurobi
#  GUROBI_INCLUDE_DIRS - The Gurobi include directories
#  GUROBI_LIBRARIES - The libraries needed to use Gurobi


set (GUROBI_BASE "c:" CACHE PATH "Base path of your gurobi installation")

if (GUROBI_INCLUDE_DIR)
  # in cache already
  set(GUROBI_FOUND TRUE)
  set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}" )
  set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}" )
else (GUROBI_INCLUDE_DIR)

  

find_path(GUROBI_INCLUDE_DIR 
          NAMES gurobi_c++.h
          PATHS "$ENV{GUROBI_HOME}/include"
                 "C:\\gurobi752\\include"
                 "${GUROBI_BASE}/include"
          )

find_library( GUROBI_LIBRARY 
              NAMES gurobi
                    gurobi75
			        gurobi60
                    gurobi56
                    gurobi55
                    gurobi51
                    gurobi50 
        		    gurobi46
				    gurobi45

              PATHS "$ENV{GUROBI_HOME}/lib" 
                    "C:\\gurobi752\\lib"
				    "${GUROBI_BASE}/lib"
              )

  if ( CMAKE_GENERATOR MATCHES "^Visual Studio 14.*")
    SET(GUROBI_LIB_NAME "gurobi_c++md2015")
  endif()
  
find_library( GUROBI_CXX_LIBRARY 
              NAMES gurobi_c++ 
			        ${GUROBI_LIB_NAME}
              PATHS "$ENV{GUROBI_HOME}/lib" 
                    "C:\\gurobi752\\lib"
					"${GUROBI_BASE}/lib"
              )
			  
# Binary dir for DLLs			
find_path(GUROBI_BIN_DIR 
                NAMES "gurobi75.dll" 
                PATHS "${GUROBI_INCLUDE_DIR}/../bin"
				      "${GUROBI_BASE}/bin"
                DOC "Directory containing the GUROBI DLLs"
               ) 		  

set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}" )
set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}" )

# use c++ headers as default
# set(GUROBI_COMPILER_FLAGS "-DIL_STD" CACHE STRING "Gurobi Compiler Flags")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBCPLEX_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GUROBI  DEFAULT_MSG
                                  GUROBI_CXX_LIBRARY GUROBI_LIBRARY GUROBI_INCLUDE_DIR)

mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY GUROBI_CXX_LIBRARY GUROBI_BIN_DIR )

endif(GUROBI_INCLUDE_DIR)