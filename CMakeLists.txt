SET(PROJECT_NAME SophusDistributions)

PROJECT(${PROJECT_NAME})
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET( CMAKE_VERBOSE_MAKEFILE ON)

################################################################################
# Release by default
# Turn on Debug with "-DCMAKE_BUILD_TYPE=Debug"
IF( NOT CMAKE_BUILD_TYPE )
   SET( CMAKE_BUILD_TYPE Release )
ENDIF()

IF( CMAKE_COMPILER_IS_GNUCXX )
   SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")
   SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

   ADD_DEFINITIONS("-Wall -Werror")
ENDIF()

################################################################################
# Create variables used for exporting (from build tree) in SophusConfig.cmake
SET( ${PROJECT_NAME}_LIBRARIES "" )
SET( ${PROJECT_NAME}_INCLUDE_DIR ${PROJECT_SOURCE_DIR} )

################################################################################
FIND_PACKAGE( Sophus REQUIRED )
LIST( APPEND ${PROJECT_NAME}_INCLUDE_DIR ${Sophus_INCLUDE_DIR} )

################################################################################
SET( SOURCE_DIR "sophus")


SET( SOURCES  ${SOURCE_DIR}/SophusAdaptors.hpp
              ${SOURCE_DIR}/distributions/Distributions.hpp
              ${SOURCE_DIR}/distributions/NormalDistributionOn.hpp
              ${SOURCE_DIR}/distributions/NormalDistributionSample.hpp
              ${SOURCE_DIR}/distributions/NormalDistributionConfidenceOperations.hpp
              ${SOURCE_DIR}/distributions/NormalDistributionConversions.hpp
              ${SOURCE_DIR}/distributions/MixtureOf.hpp
              ${SOURCE_DIR}/distributions/MixtureOfSample.hpp
)


# Add sources as custom target so that they are shown in IDE's
ADD_CUSTOM_TARGET( lib${PROJECT_NAME} SOURCES ${SOURCES} )


################################################################################
# Export package for use from the build tree
EXPORT( PACKAGE ${PROJECT_NAME} )

# Create the SophusConfig.cmake file for other cmake projects.
# ... for the build tree
SET( CONFIG_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET( CONFIG_DIR ${CMAKE_CURRENT_BINARY_DIR})
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake @ONLY )
# ... for the install tree
SET( CMAKECONFIG_INSTALL_DIR lib/cmake/${PROJECT_NAME} )
FILE( RELATIVE_PATH REL_INCLUDE_DIR
    "${CMAKE_INSTALL_PREFIX}/${CMAKECONFIG_INSTALL_DIR}"
    "${CMAKE_INSTALL_PREFIX}/include" )

SET( Sophus_INCLUDE_DIR "\${${PROJECT_NAME}_CMAKE_DIR}/${REL_INCLUDE_DIR}" )
SET( CONFIG_SOURCE_DIR )
SET( CONFIG_DIR )
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}Config.cmake @ONLY )

INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}Config.cmake"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR} )

# Install headers
INSTALL(DIRECTORY sophus DESTINATION ${CMAKE_INSTALL_PREFIX}/include
        FILES_MATCHING PATTERN "*.hpp" )
