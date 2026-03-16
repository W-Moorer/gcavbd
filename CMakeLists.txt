cmake_minimum_required(VERSION 3.13)

# Project Name
project(avbd_demo3d)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(AVBD_CORE_SOURCES
    source/collide.cpp
    source/force.cpp
    source/joint.cpp
    source/manifold.cpp
    source/rigid.cpp
    source/solver.cpp
    source/spring.cpp
)

add_executable(${PROJECT_NAME} source/main.cpp ${AVBD_CORE_SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE source)

add_executable(avbd_benchmark source/benchmark.cpp ${AVBD_CORE_SOURCES})
target_include_directories(avbd_benchmark PRIVATE source)

if(MSVC)
    target_compile_options(avbd_benchmark PRIVATE /MT$<$<CONFIG:Debug>:d>)
    set_target_properties(${PROJECT_NAME} PROPERTIES VS_GLOBAL_VcpkgEnabled false)
    set_target_properties(avbd_benchmark PROPERTIES VS_GLOBAL_VcpkgEnabled false)
endif()

# Include SDL2 first (ImGui needs it)
add_subdirectory(external/SDL)

# Add ImGui as a separate project
set(IMGUI_PROJECT_NAME "imgui")

file(GLOB IMGUI_SRC
    external/imgui/*.cpp
    external/imgui/backends/imgui_impl_sdl2.cpp
    external/imgui/backends/imgui_impl_opengl3.cpp
)

add_library(${IMGUI_PROJECT_NAME} STATIC ${IMGUI_SRC})

# Set include directories for ImGui
target_include_directories(${IMGUI_PROJECT_NAME} PUBLIC
    external/imgui
    external/imgui/backends
    external/SDL/include
)

# Link ImGui with SDL2
target_link_libraries(${IMGUI_PROJECT_NAME} PUBLIC SDL2::SDL2 SDL2::SDL2main)

# Link the main project with SDL2 and ImGui
target_include_directories(${PROJECT_NAME} PRIVATE external/SDL/include)
target_include_directories(${PROJECT_NAME} PRIVATE external/SDL/include/SDL2)
target_link_libraries(${PROJECT_NAME} PRIVATE SDL2::SDL2 SDL2::SDL2main ${IMGUI_PROJECT_NAME})

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -s USE_SDL=2 -s USE_WEBGL2=1 -s ALLOW_MEMORY_GROWTH=1 -s SINGLE_FILE=1 -s LEGACY_GL_EMULATION=1 --shell-file ../source/shell.html")
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
else()
    # Find and Link OpenGL
    find_package(OpenGL REQUIRED)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenGL::GL)
endif()

# Copy SDL2.dll to the output folder on Windows
if(WIN32)
    add_custom_command(
        TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:SDL2::SDL2>    # Path to the built SDL2.dll
            $<TARGET_FILE_DIR:${PROJECT_NAME}> # Output directory of the executable
    )
endif()

include(CTest)

if(BUILD_TESTING)
    set(BENCHMARK_TEST_CSV "${CMAKE_BINARY_DIR}/test/benchmark_smoke.csv")

    add_test(
        NAME benchmark_smoke
        COMMAND avbd_benchmark --scene 9 --steps 10 --output ${BENCHMARK_TEST_CSV}
    )

    set_tests_properties(benchmark_smoke PROPERTIES FIXTURES_SETUP benchmark_csv_ready)

    add_test(
        NAME benchmark_csv_validate
        COMMAND ${CMAKE_COMMAND}
            -DCSV_PATH=${BENCHMARK_TEST_CSV}
            -P ${CMAKE_SOURCE_DIR}/test/validate_benchmark_csv.cmake
    )

    set_tests_properties(benchmark_csv_validate PROPERTIES FIXTURES_REQUIRED benchmark_csv_ready)
endif()
