//You seem to be lucky since i had the same problem some time ago
//and written an equivalent for gcc . Note that this is not tested much
//and maybe does not has exactly the same functionality as the
//dos/windows function .
//I have included the other functions that i had to replace too ,
//maybe they will be of some use for you later.

// System:     Linux
// *******************************************************************

#ifdef __linux__
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

// Abstract:   split a path into its parts
// Parameters: Path: Object to split
//             Drive: Logical drive , only for compatibility , not considered
//             Directory: Directory part of path
//             Filename: File part of path
//             Extension: Extension part of path (includes the leading point)
// Returns:    Directory Filename and Extension are changed
// Comment:    Note that the concept of an extension is not available in Linux,
//             nevertheless it is considered

_splitpath(const char* Path, char* Drive, char* Directory, char* Filename, char* Extension)
{
  char* CopyOfPath = (char*) Path;
  int Counter = 0;
  int Last = 0;
  int Rest = 0;

  // no drives available in linux .
  // extensions are not common in linux
  // but considered anyway
  Drive = NULL;

  while(*CopyOfPath != '\0')
    {
      // search for the last slash
      while(*CopyOfPath != '/' && *CopyOfPath != '\0')
        {
          CopyOfPath++;
          Counter++;
        }
      if(*CopyOfPath == '/')
        {
          CopyOfPath++;
         Counter++;
          Last = Counter;
        }
      else
          Rest = Counter - Last;
    }
  // directory is the first part of the path until the
  // last slash appears
  strncpy(Directory,Path,Last);
  // strncpy doesnt add a '\0'
  Directory[Last] = '\0';
  // Filename is the part behind the last slahs
  strcpy(Filename,CopyOfPath -= Rest);
  // get extension if there is any
  while(*Filename != '\0')
  {
    // the part behind the point is called extension in windows systems
    // at least that is what i thought apperantly the '.' is used as part
    // of the extension too .
    if(*Filename == '.')
      {
        while(*Filename != '\0')
        {
          *Extension = *Filename;
          Extension++;
          Filename++;
        }
      }
      if(*Filename != '\0')
        {Filename++;}
  }
  *Extension = '\0';
  return;
}

#endif

