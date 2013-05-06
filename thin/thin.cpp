#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkConnectedThresholdImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkBinaryThinningImageFilter3D.h"

#include <iostream>
#include <stdlib.h>   // for atoi()
using namespace std;

int main(int argc, char* argv[])
{
  	FILE *fp;
	char errfile[] = "error.log";
	fp = fopen(errfile,"w");
	if( argc != 3 )
	{
		printf("Usage: thin  inputImageFile outputImageFile" );
		fprintf(fp,"Usage: thin  inputImageFile outputImageFile" );
 		fprintf(fp,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fp,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fp);
		return 1;
	}
  char* infilename  = argv[1];
  char* outfilename = argv[2];

  const   unsigned int Dimension = 3;
//  typedef signed short PixelType;   // must be signed for CT since Hounsfield units can be < 0
  typedef unsigned char PixelType;   // no CT images
  typedef itk::Image< PixelType, Dimension > ImageType;

  // Read image
  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( infilename );
  try
  {
  	printf("Loading input TIFF: %s\n",infilename);
	reader->Update();
  }
  catch (itk::ExceptionObject &ex)
  {
    std::cout << ex << std::endl;
 	fprintf(fp,"Read error on input file\n");
	fclose(fp);
   return 2;
  }
  cout << infilename << " successfully read." << endl;

  // Define the thinning filter 
  typedef itk::BinaryThinningImageFilter3D< ImageType, ImageType > ThinningFilterType;
  ThinningFilterType::Pointer thinningFilter = ThinningFilterType::New();
  thinningFilter->SetInput( reader->GetOutput() );
  thinningFilter->Update();

  // output to file
  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
  writer->SetInput( thinningFilter->GetOutput() );
  writer->SetFileName( outfilename );

  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject &ex)
  {
    std::cout << ex << std::endl;
 	fprintf(fp,"Write error on output file\n");
	fclose(fp);
    return 3;
  }
  cout << outfilename << " successfully written." << endl;

  cout << "Program terminated normally." << endl;
  return 0;
}


