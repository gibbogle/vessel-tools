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
	char *imagefile, *bdryfile;
	char bdryimagefile[] = "bdry.tif";
	int radius;
	int width, height, depth, count;

  if( argc != 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " <inputImage> <outputBdryData> ";
    std::cerr << std::endl;
    return 1;
    }
  imagefile =  argv[1] ;
  bdryfile = argv[2];
//  radius =  atoi(argv[3]) ;
  radius = 1;
  printf("Lymphatic image file: %s\n",imagefile);
  printf("Boundary data file: %s\n",bdryfile);
  printf("Ball radius: %d\n",radius);

  typedef unsigned char PixelType;
  typedef itk::Image<unsigned char,3> ImageType;
  const unsigned int Dimension = 3;

//  typedef itk::Image< PixelType, Dimension >    ImageType;
	typedef itk::ImageFileReader< ImageType >     ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(imagefile);
	ImageType::Pointer im;
	im = reader->GetOutput();
	reader->Update();
	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	p = (unsigned char *)(im->GetBufferPointer());
	count = 0;
	for (int i=0; i<width*height*depth; i++) {
		if (p[i] != 0) count++;
	}
	printf("Lit voxel count: %d\n",count);
//	return 0;

  typedef itk::BinaryBallStructuringElement< PixelType, Dimension > StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();

  printf("BinaryErodeImageFilter\n");
  typedef itk::BinaryErodeImageFilter< ImageType, ImageType, StructuringElementType > BinaryErodeImageFilterType;
  BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
  erodeFilter->SetInput( reader->GetOutput() );
  erodeFilter->SetKernel( structuringElement );

  printf("InvertIntensityImageFilter\n");
  typedef itk::InvertIntensityImageFilter< ImageType, ImageType > InvertIntensityImageFilterType;
  InvertIntensityImageFilterType::Pointer invertFilter = InvertIntensityImageFilterType::New();
//  invertFilter->SetInput( reader->GetOutput() );
  invertFilter->SetInput(  erodeFilter->GetOutput() );

  printf("AndImageFilter\n");
  typedef itk::AndImageFilter< ImageType, ImageType > AndImageFilterType;
  AndImageFilterType::Pointer andFilter = AndImageFilterType::New();
  andFilter->SetInput( 0, reader->GetOutput() );
  andFilter->SetInput( 1, invertFilter->GetOutput() );

//	ImageType::Pointer im;
	im = andFilter->GetOutput();

	try
    {
		andFilter->Update();
    }
    catch( itk::ExceptionObject & e )
    {
		std::cerr << "Error: " << e << std::endl;
		return 2;
    }


	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	int imsize = width*height;

	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	p = (unsigned char *)(im->GetBufferPointer());

	printf("Writing boundary data file\n");
	FILE *fpout = fopen(bdryfile,"w");
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

	printf("Writing boundary image file\n");
  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
//  writer->SetInput( erodeFilter->GetOutput() );
//  writer->SetInput( invertFilter->GetOutput() );
  writer->SetInput( im );
  writer->SetFileName(bdryimagefile );
  writer->UseCompressionOn();

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Error: " << e << std::endl;
    return 3;
    }

  return 0;
}
