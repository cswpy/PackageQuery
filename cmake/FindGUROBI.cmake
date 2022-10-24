#set (GUROBI_HOME "C:/gurobi912/win64")
set(GUROBI_HOME "/home/pw1287/downloads/gurobi950/linux64")

message("CUSTOM GUROBI HOME: " ${GUROBI_HOME})

find_path(GUROBI_INCLUDE_DIR
    NAMES gurobi_c.h
    HINTS ${GUROBI_DIR} ${GUROBI_HOME}
    PATH_SUFFIXES include)

find_library(GUROBI_LIBRARY
    NAMES gurobi gurobi95
    HINTS ${GUROBI_DIR} ${GUROBI_HOME}
    PATH_SUFFIXES lib)

if(CXX)
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
        
        find_library(GUROBI_CXX_LIBRARY
            NAMES gurobi_c++${M_FLAG}${MSVC_YEAR}
            HINTS ${GUROBI_DIR} ${GUROBI_HOME}
            PATH_SUFFIXES lib)
        find_library(GUROBI_CXX_DEBUG_LIBRARY
            NAMES gurobi_c++${M_FLAG}d${MSVC_YEAR}
            HINTS ${GUROBI_DIR} ${GUROBI_HOME}
            PATH_SUFFIXES lib)
    else()
        find_library(GUROBI_CXX_LIBRARY
            NAMES gurobi_c++
            HINTS ${GUROBI_DIR} ${GUROBI_HOME}
            PATH_SUFFIXES lib)
        set(GUROBI_CXX_DEBUG_LIBRARY ${GUROBI_CXX_LIBRARY})
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY)
