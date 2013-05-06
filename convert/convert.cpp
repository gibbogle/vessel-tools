//--------------------------------------------------------------------------------------
// To convert a tiff file to unsigned 8-bit format
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
typedef itk::Image<signed short,3> ImageType_s16;
typedef itk::Image<unsigned short,3> ImageType_u16;
ImageType_u8::Pointer im_u8;
ImageType_s16::Pointer im_s16;
ImageType_u16::Pointer im_u16;
unsigned char *p_u8;
signed short *p_s16;
unsigned short *p_u16;

#define V_u8(a,b,c)  p_u8[(c)*xysize+(b)*width+(a)]
#define V_s16(a,b,c)  p_s16[(c)*xysize+(b)*width+(a)]
#define V_u16(a,b,c)  p_u16[(c)*xysize+(b)*width+(a)]
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

	if (argc != 4) {
		printf("Usage: convert input_tiff output_tiff signed\n");
		printf("       16-bit input_tiff is converted to 8-bit output_tiff\n");
		printf("       signed = 'S' if input_tiff is signed 16-bit, 'U' otherwise\n");
//		printf("       comp = 'C' to compress image, 'U' otherwise\n");
		return 1;
	}
	
	printf("Input image file: %s\n",argv[1]);
	sscanf(argv[3],"%c",&csigned);
//	if (comp == 'C') {
//		use_compression = true;
//	} else {
//		use_compression = false;
//	}

//	_splitpath(argv[1],drive,dir,filename,ext);

	if (csigned =='S' || csigned == 's')
		issigned = true;
	else
		issigned = false;

	if (issigned) {
		typedef itk::ImageFileReader<ImageType_s16> FileReaderType_s16;
		FileReaderType_s16::Pointer reader = FileReaderType_s16::New();
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

		im_s16 = reader->GetOutput();
		width = im_s16->GetLargestPossibleRegion().GetSize()[0];
		height = im_s16->GetLargestPossibleRegion().GetSize()[1];
		depth = im_s16->GetLargestPossibleRegion().GetSize()[2];
		p_s16 = (signed short *)(im_s16->GetBufferPointer());
	} else {
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
			return 3;
		}

		im_u16 = reader->GetOutput();
		width = im_u16->GetLargestPossibleRegion().GetSize()[0];
		height = im_u16->GetLargestPossibleRegion().GetSize()[1];
		depth = im_u16->GetLargestPossibleRegion().GetSize()[2];
		p_u16 = (unsigned short *)(im_u16->GetBufferPointer());
	}
	xysize = width*height;
	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	ImageType_u8::Pointer im_u8 = ImageType_u8::New();
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
	im_u8->SetRegions(imregion);
	im_u8->Allocate();
	p_u8 = (unsigned char *)(im_u8->GetBufferPointer());

	for (x=0; x<width; x++) {
		for (y=0; y<height; y++) {
			for (z=0; z<depth; z++) {
				if (issigned)
					V_u8(x,y,z) = V_s16(x,y,z)/256;
				else
					V_u8(x,y,z) = V_u16(x,y,z)/256;
			}
		}
	}
	printf("Writing 8-bit file: %s  dimensions: width, height, depth: %d %d %d\n",argv[2],width,height,depth);
	typedef itk::ImageFileWriter<ImageType_u8> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(argv[2]);
	writer->SetInput(im_u8);
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
	printf("Wrote 8-bit file\n");
	
	return 0;
}