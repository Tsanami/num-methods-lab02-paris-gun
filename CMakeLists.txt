cmake_minimum_required(VERSION 3.10)

# Определение проекта и стандартов
project(ParishGunModel VERSION 1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Сборка библиотек и исполняемого файла
add_library(geodesy geodesy.cpp)
add_library(atmosphere atmosphere.cpp)
add_executable(run_test main.cpp)

target_link_libraries(run_test PRIVATE geodesy atmosphere)

# --- Тестирование (опционально) ---
# Для включения тестов GoogleTest выполните следующие шаги:
# 1. Добавьте в корень проекта папку googletest с исходниками GoogleTest.
#    Пример: https://github.com/google/googletest
# 2. Раскомментируйте строки ниже:
enable_testing()
add_subdirectory(googletest)
add_executable(geodesy_test tests/geodesy_test.cpp)
target_link_libraries(geodesy_test PRIVATE geodesy gtest_main)
add_test(NAME GeodesyTest COMMAND geodesy_test)
# 3. Повторите шаги сборки: mkdir build && cd build && cmake .. && cmake --build .
# 4. Запустите: ctest   или   ./geodesy_test