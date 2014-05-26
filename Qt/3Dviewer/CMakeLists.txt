PROJECT(3DViewer)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

IF(COMMAND CMAKE_POLICY)
	   CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
	INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
	MESSAGE(FATAL_ERROR "Cannot build InsightApplications without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE (VTK)
IF (VTK_FOUND)
   INCLUDE(${USE_VTK_FILE})
ENDIF(VTK_FOUND)

# Find QT 
FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE( ${QT_USE_FILE} )

# support for out-of-source build
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR} #this is where ui_3DViewerUI.h is generated
  ${CMAKE_CURRENT_SOURCE_DIR}
)
 
# Set your files and resources here
SET(3DViewer_Srcs 3Dviewer.cpp imagewidget.cpp main.cpp)
    
SET(3DViewer_Headers 3Dviewer.h imagewidget.h)
    
SET(3DViewer_UI mainwindow.ui)
 
# for generate qt additional files 
QT4_WRAP_UI(UISrcs ${3DViewer_UI})
QT4_WRAP_CPP(MOCSrcs ${3DViewer_Headers} )
 
# Packaging source code in same files 
###SOURCE_GROUP("Resources" FILES
SOURCE_GROUP("Forms" FILES
  ${3DViewer_UI}
)
 
SOURCE_GROUP("Generated" FILES
  ${UISrcs}
  ${MOCSrcs}
)
 
ADD_EXECUTABLE(3DViewer ${3DViewer_Srcs} ${3DViewer_Headers} ${UISrcs} ${MOCSrcs})

TARGET_LINK_LIBRARIES(3DViewer QVTK ${VTK_LIBRARIES} ${ITK_LIBRARIES})

