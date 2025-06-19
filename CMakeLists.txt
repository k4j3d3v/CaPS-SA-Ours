
# =============================================================================
# Minimimum CMake version.
cmake_minimum_required(VERSION 3.14)
include(CheckCXXSourceCompiles)

# Project name and languages used.
set(PROJECT_NAME caps_sa)
project(${PROJECT_NAME}
        VERSION 0.1.0
        LANGUAGES CXX
)

set(CAPSA_LIB_PROJECT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CAPSA_LIB_PROJECT_PATH}/cmake/Modules/")

# Language standards and hard requirements for such.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


if(ASAN OR UBSAN)
    set(CLANG TRUE)
endif()

if(CLANG)
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)
    # set(CMAKE_LINKER_TYPE LLD)
endif()

# Accumulate the compile- and link-flags.
set(COMPILE_FLAGS "")
set(LINK_FLAGS "")

if(ASAN)
    list(APPEND COMPILE_FLAGS ${ASAN_FLAGS} ${DEBUG_FLAGS})
    list(APPEND LINK_FLAGS ${ASAN_LINK_FLAGS} ${ASAN_FLAGS})
endif()

if(UBSAN)
    list(APPEND COMPILE_FLAGS ${UBSAN_FLAGS} ${DEBUG_FLAGS})
    list(APPEND LINK_FLAGS ${UBSAN_LINK_FLAGS} ${UBSAN_FLAGS})
endif()

if(DEBUG_INFO)
    list(APPEND COMPILE_FLAGS ${DEBUG_FLAGS})
endif()


# Add the required preprocessor definitions (`#define`s) to pass on.
add_compile_definitions(CMAKE_EXPORT_COMPILE_COMMANDS)  # For better `clangd` scanning.

# Detect OS
if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(OS_NAME "macOS")
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    set(OS_NAME "Linux")
else()
    set(OS_NAME "Unknown")
endif()

# Detect architecture
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^[xX]86_64$")
    set(ARCHITECTURE "x86_64")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)$")
    if(${OS_NAME} STREQUAL "macOS")
        set(ARCHITECTURE "AppleSilicon")
    else()
        set(ARCHITECTURE "ARM")
    endif()
else()
    set(ARCHITECTURE "Unknown")
endif()

if (ARCHITECTURE STREQUAL "x86_64")
    # Check for AVX2 support by compiling a small test program
    set(CMAKE_REQUIRED_FLAGS "-mavx2")
    check_cxx_source_compiles("
    #include <immintrin.h>
    int main() {
        __m256i x = _mm256_set1_epi32(0);
        (void)x;
        return 0;
    }
    " HAS_AVX2)

    # Reset flags to avoid affecting other checks
    unset(CMAKE_REQUIRED_FLAGS)

    # Report result
    if (HAS_AVX2)
    message(STATUS "AVX2 is supported by the compiler and system.")
    else()
        message(FATAL_ERROR "You are compiling on an x86_64 system without AVX2 support, this is not supported.")
    endif()
endif()


# Print results
message(STATUS "Operating System: ${OS_NAME}")
message(STATUS "Architecture: ${ARCHITECTURE}")

# Set other compile-time definitions.
if(OS_NAME STREQUAL "Linux")
    execute_process(
        COMMAND getconf LEVEL1_DCACHE_LINESIZE
        COMMAND tr -d '\n'
        OUTPUT_VARIABLE L1_CACHE_LINE_SIZE
    )
elseif(OS_NAME STREQUAL "macOS")
    if(ARCHITECTURE STREQUAL "x86_64")
        execute_process(
            COMMAND sysctl machdep.cpu.cache.linesize
            COMMAND awk "{print $2}"
            COMMAND tr -d '\n'
            OUTPUT_VARIABLE L1_CACHE_LINE_SIZE
        )
    elseif(ARCHITECTURE STREQUAL "AppleSilicon")
        execute_process(
            COMMAND sysctl hw.cachelinesize
            COMMAND awk "{print $2}"
            COMMAND tr -d '\n'
            OUTPUT_VARIABLE L1_CACHE_LINE_SIZE
        )
    else()
        message(FATAL_ERROR "Unable to identify Apple Silicon processor architecture")
    endif()
endif()

add_compile_definitions(L1_CACHE_LINE_SIZE=${L1_CACHE_LINE_SIZE})

set(NEED_SIMDE FALSE)
if ((ARCHITECTURE STREQUAL "ARM") OR (ARCHITECTURE STREQUAL "AppleSilicon"))
    add_compile_definitions(USE_SIMDE)
    set(NEED_SIMDE TRUE)
else()
    set(NEED_SIMDE FALSE)
endif()

if(CAPSSA_QUIET)
    message("Compiling CaPS-SA with silenced logging.")
    add_compile_definitions(CAPSSA_QUIET)
endif()

### set compiler flags based on what we learned about the machine!
# Warning flags.
set(WARNING_FLAGS -Wall -Wextra -Wpedantic)
set(SUPPRESS_WARNING_FLAGS -Wno-self-assign-overloaded)

# Bundle debugging information flags that we want to conditionally
# pass on to the compiler.
# Reference: https://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html.
set(DEBUG_FLAGS -ggdb)

# Bundle address sanitization flags: https://github.com/google/sanitizers/wiki/AddressSanitizer.
set(ASAN_FLAGS -fsanitize=address -fno-omit-frame-pointer)
# set(ASAN_LINK_FLAGS -static-libasan)

# Bundle undefined behavior sanitization flags: https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html.
set(UBSAN_FLAGS -fsanitize=undefined
                -fsanitize=float-divide-by-zero #-fsanitize=unsigned-integer-overflow
                -fsanitize=implicit-conversion -fsanitize=local-bounds -fsanitize=nullability)
# set(UBSAN_LINK_FLAGS -fuse-ld=lld)

set(OPTIMIZE_FLAGS -Ofast)

# Extra optimization flags not associated to the `-O` levels.
if (ARCHITECTURE STREQUAL "x86_64")
  list(APPEND OPTIMIZE_FLAGS -mavx2)
endif()

list(APPEND OPTIMIZE_FLAGS -funroll-loops -march=native)

message("OPTIMIZE FLAGS = ${OPTIMIZE_FLAGS}")

# Specify the warnings and the extra optimization flags to the compiler for the target library.
list(APPEND COMPILE_FLAGS ${WARNING_FLAGS} ${SUPPRESS_WARNING_FLAGS} ${OPTIMIZE_FLAGS})



# Determine the thread library of the system.
include(FindThreads)
if(NOT Threads_FOUND)
    message(FATAL_ERROR "A supported threads library is required. Aborting.")
endif()

# External projects installer.
include(ExternalProject)

set(EXT_INCLUDE ${CAPSA_LIB_PROJECT_PATH}/external/include)
file(MAKE_DIRECTORY ${EXT_INCLUDE})


# Prepare the `parlay` library—a parallel algorithms programming toolkit.
message("Build system will fetch and install parlaylib")
ExternalProject_Add(prj_parlaylib
    DOWNLOAD_DIR        ${CAPSA_LIB_PROJECT_PATH}/external
    DOWNLOAD_COMMAND    git clone https://github.com/cmuparlay/parlaylib.git
    SOURCE_DIR          ${CAPSA_LIB_PROJECT_PATH}/external/parlaylib
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CAPSA_LIB_PROJECT_PATH}/external
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       ""
    INSTALL_COMMAND     cp -rf include/parlay ${EXT_INCLUDE}
)

if (NEED_SIMDE)
# Prepare the `simde` library—provides implementations of SIMD instructions for
# systems without native SIMD support.
set(SIMDE_MIN_VERSION "0.8.2")
message("Build system will fetch and install simde")
ExternalProject_Add(prj_simde
    DOWNLOAD_DIR        ${CAPSA_LIB_PROJECT_PATH}/external
    DOWNLOAD_COMMAND    curl -k -L https://github.com/simd-everywhere/simde/releases/download/v${SIMDE_MIN_VERSION}/simde-amalgamated-${SIMDE_MIN_VERSION}.tar.xz -o simde-amalgamated-${SIMDE_MIN_VERSION}.tar.xz &&
                        tar -xf simde-amalgamated-${SIMDE_MIN_VERSION}.tar.xz &&
                        rm simde-amalgamated-${SIMDE_MIN_VERSION}.tar.xz

    SOURCE_DIR          ${CAPSA_LIB_PROJECT_PATH}/external/simde-amalgamated-${SIMDE_MIN_VERSION}
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CAPSA_LIB_PROJECT_PATH}/external
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       ""
    INSTALL_COMMAND     cp -rf x86 ${EXT_INCLUDE}
)
endif()

# Default build type.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build Type" FORCE)
endif()


# Add `src` to the build.
add_subdirectory(src)
