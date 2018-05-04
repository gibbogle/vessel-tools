//--------------------------------------------------------------------------------------
// To convert a 16-bit grayscale tiff file to a 8-bit tiff
//--------------------------------------------------------------------------------------
#include <cstdio>
#include <vector>
/*
#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>
*/
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSize.h"

typedef itk::Image<unsigned char,3> ImageType_u8;
typedef itk::Image<unsigned short,3> ImageType_u16;
//typedef itk::Image<short,3> ImageType_u16;
ImageType_u8::Pointer im_8;
ImageType_u16::Pointer im_16;
unsigned char *p_8;
unsigned short *p_16;
//short *p_16;

#define V_16(a,b,c)  p_16[(c)*xysize+(b)*width+(a)]
#define V_8(a,b,c)  p_8[(c)*xysize+(b)*width+(a)]
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int main(int argc, char**argv)
{
	int x, y, z, nbits, thresh;
	long long width, height, depth, xysize;
	bool use_compression = true;

	if (argc != 5) {
		printf("Usage: convert_16-8bit input_tiff_16 output_tiff_8 Nb thresh\n");
		printf("       16-bit grayscale input_tiff with data stored in Nb-bits is converted to 8-bit output_tiff\n");
		printf("       8-bit values < thresh are set to 0\n");
		return 1;
	}
	
	printf("Input image file: %s\n",argv[1]);
	sscanf(argv[3],"%d",&nbits);
	sscanf(argv[4],"%d",&thresh);
/*
  using ReaderType = itk::ImageFileReader< ImageType >;
  auto reader = ReaderType::New();
  reader->SetFileName( inputImage );
  // Defense against dodgy tiff images (wrong metadata)
  // Read incorrectly with SCIFIO module.
  auto tiffIO = itk::TIFFImageIO::New();
  if(tiffIO->CanReadFile(inputImage.c_str()))
    reader->SetImageIO( tiffIO );
*/
	typedef itk::ImageFileReader<ImageType_u16> FileReaderType_u16;
	FileReaderType_u16::Pointer reader = FileReaderType_u16::New();
	reader->SetFileName(argv[1]);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 2;
	}

	im_16 = reader->GetOutput();
	width = im_16->GetLargestPossibleRegion().GetSize()[0];
	height = im_16->GetLargestPossibleRegion().GetSize()[1];
	depth = im_16->GetLargestPossibleRegion().GetSize()[2];
	xysize = width*height;
	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);
	reader->GetImageIO()->Print( std::cout );

	p_16 = (unsigned short *)(im_16->GetBufferPointer());

	ImageType_u8::Pointer im_8 = ImageType_u8::New();
	ImageType_u8::SizeType imsize; 
	ImageType_u8::IndexType imstart; 
	ImageType_u8::RegionType imregion; 

	imsize[0] = width;
	imsize[1] = height;
	imsize[2] = depth;
	imstart[0] = 0;
	imstart[1] = 0;
	imstart[2] = 0;
	imregion.SetSize(imsize);
	imregion.SetIndex(imstart);
	im_8->SetRegions(imregion);
	im_8->Allocate();
	p_8 = (unsigned char *)(im_8->GetBufferPointer());

	printf("Converting...\n");
	int vmin = 9999;
	int vmax = 0;
	double vscale = 255./(pow(2.,nbits) - 1);
	for (x=0; x<width; x++) {
		for (y=0; y<height; y++) {
			for (z=0; z<depth; z++) {
				unsigned short val = V_16(x,y,z);
				vmin = MIN(val,vmin);
				vmax = MAX(val,vmax);
				val *= vscale;
				if (val < thresh) val = 0;
				V_8(x,y,z) = (unsigned char)(val);
			}
		}
	}
	printf("16-bit vmin, vmax: %d %d\n",vmin,vmax);

	printf("Writing 8-bit file: %s  dimensions: width, height, depth: %d %d %d\n",argv[2],width,height,depth);
	typedef itk::ImageFileWriter<ImageType_u8> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(argv[2]);
	writer->SetInput(im_8);
	if (use_compression) {
		writer->UseCompressionOn();
	}
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 3;
	}
	printf("Wrote 8-bit file\n");
	
	return 0;
}