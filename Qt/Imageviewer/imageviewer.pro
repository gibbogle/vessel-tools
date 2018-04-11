CONFIG += uitools
CONFIG += release

INCLUDEPATH  += "C:/Program Files/ITK4.8.0/include/ITK-4.8"
INCLUDEPATH  += C:/VTK-VS10-x64/include/vtk-5.10

HEADERS       = imageviewer.h imagewidget.h
               
SOURCES       = main.cpp imageviewer.cpp imagewidget.cpp
              
LIBS += -LC:/VTK-VS10-x64/lib/vtk-5.10 -lQVTK -lvtkRendering -lvtkGraphics -lvtkImaging -lvtkIO -lvtkFiltering -lvtkCommon \
-lvtkpng -lvtktiff -lvtkjpeg -lvtkexpat -lvfw32 -lopengl32  \
-lwsock32 -lvtksys -lws2_32 -lvtkexoIIc -lvtkNetCDF \
-lvtklibxml2 -lvtkzlib -lvtkalglib \
-lgdi32 -lkernel32 -luser32 -lgdi32 -lwinspool -lshell32 -lole32 -loleaut32 -luuid -lcomdlg32 -ladvapi32

LIBS += -L"C:/Program Files/ITK4.8.0/lib" -lITKcommon-4.8 -lITKIOImageBase-4.8 -lITKIOTIFF-4.8 -lITKIOVTK-4.8 -lITKPath-4.8 \
        -litksys-4.8 -litktiff-4.8 -lITKVTK-4.8 -lITKVtkGlue-4.8 -litkzlib-4.8 -litkvnl_algo-4.8 -litkvnl-4.8 -lITKVNLInstantiation-4.8 \
        -litkNetlibSlatec-4.8 -litkv3p_netlib-4.8

#QT           += network

# install
target.path = .
sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS imageviewer.pro
sources.path = .
INSTALLS += target sources

FORMS += mainwindow.ui

QMAKE_LFLAGS += /OPT:NOREF
