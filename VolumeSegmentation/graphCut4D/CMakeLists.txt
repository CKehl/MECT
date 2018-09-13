if(CMAKE_BUILD_TYPE MATCHES "Release")
  ADD_SUBDIRECTORY(Wrapping)
endif()
set(GRAPHCUT_USHRT_PROGRAM graphCut4D_ushort)
set(GRAPHCUT_7_PROGRAM graphCut4D_N7)
set(GRAPHCUT_27_PROGRAM graphCut4D_N27)
set(GRAPHCUT_27w_PROGRAM graphCut4D_N27w)

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
	${CMAKE_CURRENT_SOURCE_DIR}/src/GraphCut4D.hpp 
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
	${CMAKE_CURRENT_SOURCE_DIR}/src/GraphCut4D.cpp 
	${ADDITIONAL_SOURCES}
)
set(GRAPHCUT_7_SOURCES
	${GRAPHCUT_SOURCES} 
	${CMAKE_CURRENT_SOURCE_DIR}/src/main_7.cpp
)
set(GRAPHCUT_27_SOURCES
        ${GRAPHCUT_SOURCES} 
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main_27.cpp
)
set(GRAPHCUT_27w_SOURCES
        ${GRAPHCUT_SOURCES} 
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main_27w.cpp
)
set(GRAPHCUT_USHRT_SOURCES
        ${GRAPHCUT_SOURCES} 
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main_ushort.cpp
)

include_directories(${INCLUDE_DIRS} ${CMAKE_INCLUDE_PATH})
link_directories(${graphCut4D_SOURCE_DIR} ${ADDITIONAL_LIBRARY_PATHS} ${CMAKE_BINARY_DIR})

add_executable(${GRAPHCUT_7_PROGRAM} ${GRAPHCUT_7_SOURCES})
target_link_libraries(${GRAPHCUT_7_PROGRAM} ${LIBRARIES})
install(TARGETS ${GRAPHCUT_7_PROGRAM} RUNTIME DESTINATION bin)

add_executable(${GRAPHCUT_27_PROGRAM} ${GRAPHCUT_27_SOURCES})
target_link_libraries(${GRAPHCUT_27_PROGRAM} ${LIBRARIES})
install(TARGETS ${GRAPHCUT_27_PROGRAM} RUNTIME DESTINATION bin)

add_executable(${GRAPHCUT_27w_PROGRAM} ${GRAPHCUT_27w_SOURCES})
target_link_libraries(${GRAPHCUT_27w_PROGRAM} ${LIBRARIES})
install(TARGETS ${GRAPHCUT_27w_PROGRAM} RUNTIME DESTINATION bin)

add_executable(${GRAPHCUT_USHRT_PROGRAM} ${GRAPHCUT_USHRT_SOURCES})
target_link_libraries(${GRAPHCUT_USHRT_PROGRAM} ${LIBRARIES})
install(TARGETS ${GRAPHCUT_USHRT_PROGRAM} RUNTIME DESTINATION bin)
