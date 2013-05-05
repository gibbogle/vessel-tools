/*=========================================================================
 To find all voxels on the boundary (surface) of a 3D binary image.
 The image is eroded, inverted, then anded with the original image
 The resulting image contains only voxels on the boundary of the original image
 These voxels are written to the boundary file.
 *=========================================================================*/
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include <itkInvertIntensityImageFilter.h>
#include <itkAndImageFilter.h>

unsigned char *p;
#define V(a,b,c)  p[(c)*imsize+(b)*width+(a)]

int main( int argc, char* argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " <inputImage> <outputImage> <radius>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }
  typedef unsigned char PixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image< PixelType, Dimension >    ImageType;
  typedef itk::ImageFileReader< ImageType >     ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  typedef itk::BinaryBallStructuringElement< PixelType, Dimension > StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius( atoi( argv[3] ) );
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryErodeImageFilter< ImageType, ImageType, StructuringElementType > BinaryErodeImageFilterType;
  BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
  erodeFilter->SetInput( reader->GetOutput() );
  erodeFilter->SetKernel( structuringElement );

  typedef itk::InvertIntensityImageFilter< ImageType, ImageType > InvertIntensityImageFilterType;
  InvertIntensityImageFilterType::Pointer invertFilter = InvertIntensityImageFilterType::New();
//  invertFilter->SetInput( reader->GetOutput() );
  invertFilter->SetInput(  erodeFilter->GetOutput() );

  typedef itk::AndImageFilter< ImageType, ImageType > AndImageFilterType;
  AndImageFilterType::Pointer andFilter = AndImageFilterType::New();
  andFilter->SetInput( 0, reader->GetOutput() );
  andFilter->SetInput( 1, invertFilter->GetOutput() );

	typedef itk::Image<unsigned char,3> ImageType;
	ImageType::Pointer im;
	im = andFilter->GetOutput();

	try
    {
		andFilter->Update();
    }
    catch( itk::ExceptionObject & e )
    {
		std::cerr << "Error: " << e << std::endl;
		return EXIT_FAILURE;
    }


	int width = im->GetLargestPossibleRegion().GetSize()[0];
	int height = im->GetLargestPossibleRegion().GetSize()[1];
	int depth = im->GetLargestPossibleRegion().GetSize()[2];
	int imsize = width*height;

	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	p = (unsigned char *)(im->GetBufferPointer());

	FILE *fpout = fopen("bdry.dat","w");

	int n=0;
	for (int x=0; x<width; x++) {
		for (int y=0; y<height; y++) {
			for (int z=0; z<depth; z++) {
				if (V(x,y,z) != 0) {
					n++;
					fprintf(fpout,"%6d %6d %6d\n",x,y,z);
				}
			}
		}
	}
	fprintf(fpout,"np: %d\n",n);
	fclose(fpout);
	printf("Number of bdry voxels: %d\n",n);


  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
//  writer->SetInput( erodeFilter->GetOutput() );
//  writer->SetInput( invertFilter->GetOutput() );
  writer->SetInput( andFilter->GetOutput() );
  writer->SetFileName( argv[2] );
  writer->UseCompressionOn();

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Error: " << e << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
