//--------------------------------------------------------------------------------------
// To convert a binary (0,1) tiff file to a (0,255) tiff
//--------------------------------------------------------------------------------------
#include <cstdio>
#include <vector>

#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSize.h"

typedef itk::Image<unsigned char,3> ImageType_u8;
ImageType_u8::Pointer im_1, im_255;
unsigned char *p_1, *p_255;

#define V_1(a,b,c)  p_1[(c)*xysize+(b)*width+(a)]
#define V_255(a,b,c)  p_255[(c)*xysize+(b)*width+(a)]
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int main(int argc, char**argv)
{
	int x, y, z, xx, yy, zz, pos;
	int width, height, depth, xysize;
	char csigned;
//	char LEfile[256], GTfile[256], tag[13];
//	char drive[32], dir[256],filename[256], ext[32];
	bool use_compression = true;
	bool issigned;

	if (argc != 3) {
		printf("Usage: convert_1-255 input_tiff_1 output_tiff_255\n");
		printf("       binary (0-1) input_tiff is converted to binary (0-255) output_tiff\n");
		return 1;
	}
	
	printf("Input image file: %s\n",argv[1]);

	typedef itk::ImageFileReader<ImageType_u8> FileReaderType_u8;
	FileReaderType_u8::Pointer reader = FileReaderType_u8::New();
	reader->SetFileName(argv[1]);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 3;
	}

	im_1 = reader->GetOutput();
	width = im_1->GetLargestPossibleRegion().GetSize()[0];
	height = im_1->GetLargestPossibleRegion().GetSize()[1];
	depth = im_1->GetLargestPossibleRegion().GetSize()[2];
	p_1 = (unsigned char *)(im_1->GetBufferPointer());
	xysize = width*height;
	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	ImageType_u8::Pointer im_255 = ImageType_u8::New();
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
	im_255->SetRegions(imregion);
	im_255->Allocate();
	p_255 = (unsigned char *)(im_255->GetBufferPointer());

	printf("Converting...\n");
	for (x=0; x<width; x++) {
		for (y=0; y<height; y++) {
			for (z=0; z<depth; z++) {
				if (V_1(x,y,z) == 1)
					V_255(x,y,z) = 255;
				else
					V_255(x,y,z) = 0;
			}
		}
	}
	printf("Writing binary 0-255 file: %s  dimensions: width, height, depth: %d %d %d\n",argv[2],width,height,depth);
	typedef itk::ImageFileWriter<ImageType_u8> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(argv[2]);
	writer->SetInput(im_255);
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
		return 5;
	}
	printf("Wrote binary 0-255 file\n");
	
	return 0;
}