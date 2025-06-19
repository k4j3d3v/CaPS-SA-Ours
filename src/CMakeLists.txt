
# =============================================================================

# Core source files.
set(PROJECT_SRC
        Suffix_Array.cpp
        Genomic_Text.cpp
)


# Core library.
add_library(core STATIC ${PROJECT_SRC})


# Include paths for the compiler.
set(INCLUDE_DIR ${CAPSA_LIB_PROJECT_PATH}/include)
target_include_directories(core PUBLIC ${INCLUDE_DIR} ${EXT_INCLUDE})


# Executable for the driver program.
add_executable(${PROJECT_NAME} main.cpp)
add_executable(test test.cpp)


# Add dependency of the core to parlaylib.
add_dependencies(core prj_parlaylib)
if (NEED_SIMDE)
  add_dependencies(core prj_simde)
endif()

# Compiler flags.
target_compile_options(core PRIVATE ${COMPILE_FLAGS})
target_compile_options(${PROJECT_NAME} PRIVATE ${COMPILE_FLAGS})
target_compile_options(test PRIVATE ${COMPILE_FLAGS})


# Link the core lib to the threads lib.
target_link_libraries(core PRIVATE Threads::Threads)

# Link the executable to the core.
target_link_libraries(${PROJECT_NAME} PRIVATE core)
target_link_libraries(test PRIVATE core)

# Linker flags.
target_link_libraries(core PRIVATE ${LINK_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${LINK_FLAGS})
target_link_libraries(test PRIVATE ${LINK_FLAGS})


# Installation directory.
install(TARGETS ${PROJECT_NAME}
        DESTINATION bin
)

