if(CMAKE_BUILD_TYPE MATCHES "Release")
  ADD_SUBDIRECTORY(Wrapping)
endif()

set(GRAPHCUT_PROGRAM graphCut3D)

set( CMAKE_EXECUTABLE_SUFFIX ${CMAKE_EXECUTABLE_SUFFIX})
set( CMAKE_SHARED_LIBRARY_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX} )
set( CMAKE_STATIC_LIBRARY_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX} )
set( CMAKE_DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX} )
set( CMAKE_RELEASE_POSTFIX ${CMAKE_RELEASE_POSTFIX} )

set(GRAPHCUT_HFILES
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/graphs/src/graphDisjointForest.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/graphs/src/graphEdge.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/image.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/image_utils.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/image2D.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/image2D_utils.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/image3D.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/image3D_utils.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/image4D.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/image4D_utils.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/mhdImage.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/io/src/io_commons.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/common/src/std_typedefs.h 
	${CMAKE_CURRENT_SOURCE_DIR}/src/GraphCut3D.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/main.hpp 
	${ADDITIONAL_HEADERS}
)
set(GRAPHCUT_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/graphs/src/graphDisjointForest.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/graphs/src/graphEdge.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/image.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/image2D.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/image3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/image4D.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/imaging/src/mhdImage.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/GraphCut3D.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp 
	${ADDITIONAL_SOURCES}
)

include_directories(${INCLUDE_DIRS} ${CMAKE_INCLUDE_PATH})
link_directories(${graphCut3D_SOURCE_DIR} ${ADDITIONAL_LIBRARY_PATHS} ${CMAKE_BINARY_DIR})

add_executable(${GRAPHCUT_PROGRAM} ${GRAPHCUT_SOURCES})
target_link_libraries(${GRAPHCUT_PROGRAM} ${LIBRARIES})

install(TARGETS ${GRAPHCUT_PROGRAM} RUNTIME DESTINATION bin)
