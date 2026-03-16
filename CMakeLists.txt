cmake_minimum_required(VERSION 3.24)

project(gcavbd LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

# Keep third-party dependencies minimal in this scaffold.
set(BUILD_TESTING OFF CACHE BOOL "Disable third-party tests" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "Disable Eigen docs" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "Disable Eigen pkg-config" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "Disable gtest install target" FORCE)
set(gtest_force_shared_crt ON CACHE BOOL "Use shared CRT on MSVC" FORCE)

FetchContent_Declare(
  eigen
  URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)

FetchContent_MakeAvailable(eigen googletest)

# Re-enable tests for this project.
set(BUILD_TESTING ON CACHE BOOL "Enable gcavbd tests" FORCE)

add_library(gcavbd STATIC
  src/core/types.cpp
)

target_include_directories(gcavbd
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(gcavbd
  PUBLIC
    Eigen3::Eigen
)

include(CTest)
enable_testing()
if(BUILD_TESTING)
  add_executable(test_math
    tests/test_math.cpp
  )

  target_link_libraries(test_math
    PRIVATE
      gcavbd
      GTest::gtest_main
  )

  add_test(NAME test_math COMMAND test_math)
endif()
