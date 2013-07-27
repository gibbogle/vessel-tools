rem Notes:
rem It is necessary to ensure that Qt Creator finds the 64-bit Qt dlls.  Put them at the front of the PATH list.

set PATH=C:\Qt64\4.8.1\bin;C:\Qt64\4.8.1\lib;%PATH%
C:\QtSDK\QtCreator\bin\qtcreator CMakeLists.txt

