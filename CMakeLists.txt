cmake_minimum_required(VERSION 3.16)
project(virtual_casing)
#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

OPTION (USE_PROFILER "Activate code profiling." OFF)
if (USE_PROFILER)
    add_definitions (-DSCTL_PROFILE=5 -DSCTL_VERBOSE)
endif ()
set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON")

find_package (BLAS)
if (BLAS_FOUND)
    add_definitions (-DSCTL_HAVE_BLAS)
    if (NOT ${BLA_VENDOR} MATCHES "^Intel")
        find_package (LAPACK)
        if(LAPACK_FOUND)
            add_definitions (-DSCTL_HAVE_LAPACK)
        endif()
    else()
        add_definitions (-DSCTL_HAVE_LAPACK)
    endif()
endif ()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(FFTW COMPONENTS FLOAT_LIB DOUBLE_LIB LONGDOUBLE_LIB DOUBLE_OPENMP_LIB)
if(FFTW_FOUND)
    add_definitions (-DSCTL_HAVE_FFTW)
    if (FFTW_FLOAT_LIB_FOUND)
        add_definitions (-DSCTL_HAVE_FFTWF)
    endif()
    if (FFTW_LONGDOUBLE_LIB_FOUND)
        add_definitions (-DSCTL_HAVE_FFTWL)
    endif()
    if (FFTW_FLOAT_OPENMP_FOUND)
        add_definitions (-DSCTL_FFTWF_THREADS)
    endif()
endif()

set(Python_FIND_STRATEGY LOCATION)
find_package(Python 3 COMPONENTS Interpreter Development.Module NumPy)

message(status "Python executable is ${Python_EXECUTABLE}")
message(status "Python Development Module found value is ${Python_Development.Module_FOUND}")
message(status "Python header dirs are ${Python_INCLUDE_DIRS}")
message(status "Python library dirs are ${Python_LIBRARY_DIRS}")
message(status "Python version is ${Python_VERSION}")


find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

IF(DEFINED ENV{CI})
    message(STATUS "CI environment detected. Set compilation flags accordingly (target ivybridge which supports avx).")
    set(CMAKE_CXX_FLAGS "-O3 -march=ivybridge -mfma -ffp-contract=fast")
elseif(DEFINED ENV{CONDA_BUILD})
    message(STATUS "conda build environment detected. Let conda set compilation flags accordingly.")
    # set(CMAKE_CXX_FLAGS "-O3 -march=ivybridge -mfma -ffp-contract=fast")
else()
    message(STATUS "Local build detected. Set compilation flags accordingly (march=native).")
    set(CMAKE_CXX_FLAGS "-O3 -march=native -mfma -ffp-contract=fast")
endif()

find_package(OpenMP REQUIRED)

add_subdirectory(thirdparty/pybind11)
# add_subdirectory(thirdparty/fmt)
# set(XTENSOR_USE_OPENMP 0)
# set(XTENSOR_USE_TBB 0)


pybind11_add_module(${PROJECT_NAME} src/python.cpp)

set_target_properties(${PROJECT_NAME}
    PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON)
add_executable(vc_testing EXCLUDE_FROM_ALL test/virtual-casing.cpp)
set_target_properties(vc_testing
    PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON)

#target_include_directories(${PROJECT_NAME} PRIVATE 
    # "thirdparty/xtensor/include" "thirdparty/xtensor-python/include" "thirdparty/xsimd/include" "thirdparty/xtl/include" "thirdparty/eigen"
#    ${Python_NumPy_INCLUDE_DIRS} "include" "BIEST/include")
# target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt)

target_include_directories(vc_testing PRIVATE  "include" "BIEST/include")
target_include_directories(${PROJECT_NAME} PRIVATE  "include" "BIEST/include")
target_include_directories(${PROJECT_NAME} PRIVATE  ${Python_NumPy_INCLUDE_DIRS})

target_link_libraries(vc_testing PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)

if(FFTW_FOUND)
    target_link_libraries (vc_testing PRIVATE FFTW::Double)
    target_link_libraries (${PROJECT_NAME} PRIVATE FFTW::Double)
    if (FFTW_FLOAT_LIB_FOUND)
        target_link_libraries (vc_testing PRIVATE FFTW::Float)
        target_link_libraries (${PROJECT_NAME} PRIVATE FFTW::Float)
    endif()
    if (FFTW_LONGDOUBLE_LIB_FOUND)
        target_link_libraries (vc_testing PRIVATE FFTW::LongDouble)
        target_link_libraries (${PROJECT_NAME} PRIVATE FFTW::LongDouble)
    endif()
    if (FFTW_FLOAT_OPENMP_FOUND)
        target_link_libraries (vc_testing PRIVATE FFTW::DoubleOpenMP)
        target_link_libraries (${PROJECT_NAME} PRIVATE FFTW::DoubleOpenMP)
    endif()
endif()
if (BLAS_FOUND)
	target_link_libraries (vc_testing PRIVATE ${BLAS_LIBRARIES})
	target_link_libraries (${PROJECT_NAME} PRIVATE ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (vc_testing PRIVATE ${LAPACK_LIBRARIES})
	target_link_libraries (${PROJECT_NAME} PRIVATE ${LAPACK_LIBRARIES})
endif ()


#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

#install(TARGETS ${PROJECT_NAME}
#        #LIBRARY
#        DESTINATION src/${PROJECT_NAME})
#install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${PROJECT_NAME})
