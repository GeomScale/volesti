# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.16.3)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(lpsolve-populate NONE)

include(ExternalProject)
ExternalProject_Add(lpsolve-populate
                     "UPDATE_DISCONNECTED" "False" "URL" "https://cran.r-project.org/src/contrib/lpSolve_5.6.15.tar.gz" "URL_HASH" "MD5=3b8d780f703e0da2e4863939add164bc"
                    SOURCE_DIR          "/home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build/_deps/lpsolve-src"
                    BINARY_DIR          "/home/agentperry/Documents/CodeFiles2/GeomScale/volesti_fork/test/build/_deps/lpsolve-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)
