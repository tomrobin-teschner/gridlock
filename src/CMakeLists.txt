target_sources(${CMAKE_PROJECT_NAME}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/gridlock.cpp
  PUBLIC FILE_SET HEADERS
  BASE_DIRS ${PROJECT_SOURCE_DIR}
)

add_subdirectory(infrastructure)
add_subdirectory(fieldArray)
add_subdirectory(timeInfo)
add_subdirectory(postProcessing)
add_subdirectory(residuals)
add_subdirectory(boundaryConditions)
add_subdirectory(mesh)
# add_subdirectory(governingEquations)