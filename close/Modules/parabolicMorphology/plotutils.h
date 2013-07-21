#ifndef __plotutils_h
#define __plotutils_h
#include "itkLineConstIterator.h"
#include <iostream>
#include <fstream>

// some simple profile extraction utilities for production of plots of
// image transects

template <typename TImage>
void extractProfile(typename TImage::Pointer Im, typename TImage::IndexType first, 
		    typename TImage::IndexType last, std::string outputfile)
{
  typedef typename itk::LineConstIterator<TImage> ItType;

  ItType it(Im, first, last);

  std::ofstream outfile;
  outfile.open(outputfile.c_str());
  unsigned pos = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++pos)
    {
    // write out the pixel value
    outfile << pos << "\t" << (float)it.Get() << std::endl;
    }
  outfile.close();

}

#endif
