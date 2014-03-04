// System:     Linux
// *******************************************************************

#ifdef __linux__
#include <stdio.h>
#include <unistd.h>		// comment this for Windows testing
#include <stdlib.h>
#include <string.h>
#include <errno.h>

// Abstract:   split a path into its parts
// Parameters: Path:      Object to split
//             Drive:     Logical drive , only for compatibility , not considered
//             Directory: Directory part of path
//             Filename:  File part of path minus extension
//             Extension: Extension part of path (includes the leading point)
// Returns:    Directory Filename and Extension are changed

void _splitpath(const char* Path, char* Drive, char* Directory, char* Filename, char* Extension)
{
  char* CopyOfPath = (char*) Path;	// We do not want to change Path
  char buffer[2048];
  char *wholeFilename;
  int Counter = 0;
  int Last = 0;
  int Rest = 0;
  wholeFilename = buffer;

  // no drives available in linux
  Drive[0] = '\0';
  while (*CopyOfPath != '\0') {
      // search for the last slash
      while (*CopyOfPath != '/' && *CopyOfPath != '\0') {
          CopyOfPath++;
          Counter++;
      }
      if (*CopyOfPath == '/') {
          CopyOfPath++;
          Counter++;
          Last = Counter;
	  } else {
          Rest = Counter - Last;
	  }
  }
  // directory is the first part of the path until the last slash appears
  strncpy(Directory,Path,Last);
  // strncpy doesnt add a '\0'
  Directory[Last] = '\0';
  // Filename is the part behind the last slash
  strcpy(wholeFilename,CopyOfPath -= Rest);
  // get extension if there is any
  while (*wholeFilename != '\0') {
    // The '.' is used as part of the extension
      if (*wholeFilename == '.') {
          while (*wholeFilename != '\0') {
              *Extension = *wholeFilename;
              Extension++;
              wholeFilename++;
          }
		  break;
	  } else {
		  *Filename = *wholeFilename;
		  Filename++;
	  }
      if (*wholeFilename != '\0')
          wholeFilename++;
  }
  *Filename = '\0';
  *Extension = '\0';
  return;
}
#endif

/*
// For testing

int main(int argv, char **argc)
{
	char path[] = "/home/gib/zzz.txt";
	char drive[16], directory[1024], filename[256], ext[32];
	char newpath[1024];

	_mysplitpath(path,drive,directory,filename,ext);
	printf("drive: %s\n",drive);
	printf("directory: %s\n",directory);
	printf("filename: %s\n",filename);
	printf("ext: %s\n",ext);
	strcpy(newpath,drive);
	strcat(newpath,directory);
	strcat(newpath,filename);
	strcat(newpath,ext);
	printf("newpath: %s\n",newpath);

	strcpy(path,"C:/dir/zzz.txt");
	_splitpath(path,drive,directory,filename,ext);
	printf("drive: %s\n",drive);
	printf("directory: %s\n",directory);
	printf("filename: %s\n",filename);
	printf("ext: %s\n",ext);
	strcpy(newpath,drive);
	strcat(newpath,directory);
	strcat(newpath,filename);
	strcat(newpath,ext);
	printf("newpath: %s\n",newpath);
}

*/