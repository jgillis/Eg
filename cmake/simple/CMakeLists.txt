cmake_minimum_required(VERSION 2.8)

project(simple)
MESSAGE(STATUS "\n\n")
MESSAGE(STATUS "binary dir: ${simple_BINARY_DIR}   : ${PROJECT_BINARY_DIR}")
MESSAGE(STATUS "source dir: ${simple_SOURCE_DIR} : ${PROJECT_SOURCE_DIR}")


find_file(LOCATION_TREASURE
treasure.txt
${PROJECT_SOURCE_DIR}/hidden
)

find_file(LOCATION_CHEST
chest.txt
${PROJECT_SOURCE_DIR}
)

IF(LOCATION_TREASURE)
MESSAGE(STATUS "Treasure found:${LOCATION_TREASURE} ")
ENDIF()

IF(NOT LOCATION_CHEST)
MESSAGE(STATUS "Chest not found: ${LOCATION_CHEST}")
ENDIF()

MESSAGE(STATUS "\n\n")