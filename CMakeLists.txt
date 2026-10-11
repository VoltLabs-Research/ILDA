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

project(ilda VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)

if(NOT TARGET coretoolkit::coretoolkit)
    find_package(coretoolkit REQUIRED)
endif()
if(NOT TARGET structure-identification::structure-identification)
    find_package(structure-identification REQUIRED)
endif()
if(NOT TARGET common-neighbor-analysis::common-neighbor-analysis)
    find_package(common-neighbor-analysis REQUIRED)
endif()
if(NOT TARGET polyhedral-template-matching::polyhedral-template-matching)
    find_package(polyhedral-template-matching REQUIRED)
endif()
if(NOT TARGET opendxa::opendxa)
    find_package(opendxa REQUIRED)
endif()
if(NOT TARGET spdlog::spdlog)
    find_package(spdlog REQUIRED)
endif()
if(NOT TARGET nlohmann_json::nlohmann_json)
    find_package(nlohmann_json REQUIRED)
endif()
if(NOT TARGET Boost::headers AND NOT TARGET Boost::boost AND NOT TARGET boost::headers AND NOT TARGET boost::boost)
    find_package(Boost REQUIRED)
endif()
if(NOT EMSCRIPTEN)
    if(NOT TARGET TBB::tbb)
        find_package(TBB REQUIRED)
    endif()
    if(NOT TARGET duckdb::duckdb)
        find_package(duckdb REQUIRED)
    endif()
endif()

set(ILDA_BOOST_TARGET "")
if(TARGET Boost::headers)
    set(ILDA_BOOST_TARGET Boost::headers)
elseif(TARGET Boost::boost)
    set(ILDA_BOOST_TARGET Boost::boost)
elseif(TARGET boost::headers)
    set(ILDA_BOOST_TARGET boost::headers)
elseif(TARGET boost::boost)
    set(ILDA_BOOST_TARGET boost::boost)
else()
    message(FATAL_ERROR "ILDA requires a Boost headers target")
endif()

file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
list(FILTER LIB_SOURCES EXCLUDE REGEX ".*main\\.cpp$")

add_library(${PROJECT_NAME}_lib STATIC ${LIB_SOURCES})
set_target_properties(${PROJECT_NAME}_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(${PROJECT_NAME}_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_precompile_headers(${PROJECT_NAME}_lib PRIVATE <volt/core/volt.h>)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC
    ${ILDA_BOOST_TARGET}
)
if(NOT EMSCRIPTEN)
    target_link_libraries(${PROJECT_NAME}_lib PUBLIC TBB::tbb)
endif()
target_link_libraries(${PROJECT_NAME}_lib PRIVATE
    opendxa::opendxa
    coretoolkit::coretoolkit
    structure-identification::structure-identification
    common-neighbor-analysis::common-neighbor-analysis
    polyhedral-template-matching::polyhedral-template-matching
    spdlog::spdlog
    nlohmann_json::nlohmann_json
)
if(NOT EMSCRIPTEN)
    target_link_libraries(${PROJECT_NAME}_lib PRIVATE duckdb::duckdb)
endif()

add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_precompile_headers(${PROJECT_NAME} PRIVATE <volt/core/volt.h>)

target_link_libraries(${PROJECT_NAME} PRIVATE
    ${PROJECT_NAME}_lib
    opendxa::opendxa
    structure-identification::structure-identification
    common-neighbor-analysis::common-neighbor-analysis
    polyhedral-template-matching::polyhedral-template-matching
    coretoolkit::coretoolkit
    spdlog::spdlog
    nlohmann_json::nlohmann_json
)

if(EMSCRIPTEN)
    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".js")
    target_compile_options(${PROJECT_NAME} PRIVATE -fexceptions)
    target_link_options(${PROJECT_NAME} PRIVATE
        "SHELL:-sALLOW_MEMORY_GROWTH=1"
        "SHELL:-sMAXIMUM_MEMORY=4GB"
        "SHELL:-sEXIT_RUNTIME=0"
        "SHELL:-sMODULARIZE=0"
        "SHELL:-sENVIRONMENT=web,node"
        "SHELL:-sFORCE_FILESYSTEM=1"
        "SHELL:-sDISABLE_EXCEPTION_CATCHING=0"
        "SHELL:-sEXPORTED_RUNTIME_METHODS=['callMain','FS']"
    )
    install(FILES "$<TARGET_FILE_DIR:${PROJECT_NAME}>/${PROJECT_NAME}.wasm" DESTINATION bin OPTIONAL)
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
install(TARGETS ${PROJECT_NAME}_lib DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include)
