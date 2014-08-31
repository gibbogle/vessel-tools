//--------------------------------------------------------------------------------------
// To cut a tiff file into four pieces, named: xxx_x0y0, xxx_x0y1, xxx_x1y0, xxx_x1y1
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
#include "itkTIFFImageIO.h"

typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im, im1, im2;
unsigned char *p, *p1, *p2;

#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)]
#define V1(a,b,c)  p1[(c)*xysize1+(b)*width1+(a)]
//#define V2(a,b,c)  p2[(c)*xysize2+(b)*width2+(a)]
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#ifdef __linux__
void _splitpath(const char* Path, char* Drive, char* Directory, char* Filename, char* Extension);
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int main(int argc, char**argv)
{
	long long int x, y, z, xx, yy, zz, pos;
	long long int width, height, depth, xysize;
	long long int width1, height1, depth1, xysize1;
	long long xstart, ystart;
//	long long int width2, height2, depth2, xysize2;
//	long long int xstart2, ystart2, zstart2;
	char axis, comp;
	char Qfile[256], tag[13];
	char drive[32], dir[256],filename[256], ext[32];
	bool use_compression=true;

	if (argc != 2) {
		printf("Usage: cut input_tiff\n");
		return 1;
	}
	
	printf("Input image file: %s\n",argv[1]);
	
	_splitpath(argv[1],drive,dir,filename,ext);
	//strcpy(LEfile,drive);
	//strcat(LEfile,dir);
	//strcat(LEfile,filename);
	//strcpy(GTfile,LEfile);
	//sprintf(tag,"_%cLE%04d.tif",axis,pos);
	//strcat(LEfile,tag);
	//sprintf(tag,"_%cGT%04d.tif",axis,pos);
	//strcat(GTfile,tag);
	//printf("LEfile: %s\n",LEfile);
	//printf("GTfile: %s\n",GTfile);

	typedef itk::ImageFileReader<ImageType> FileReaderType;
	FileReaderType::Pointer reader = FileReaderType::New();

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

	im = reader->GetOutput();

	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	xysize = width*height;
	printf("Image dimensions: width, height, depth: %ld %ld %ld\n",width,height,depth);
	p = (unsigned char *)(im->GetBufferPointer());

	//if (axis == 'x' && pos >= width-1) {
	//	printf("Bad cut position: %c = %ld\n",axis,pos);
	//	return 4;
	//}
	//if (axis == 'y' && pos >= height-1) {
	//	printf("Bad cut position: %c = %ld\n",axis,pos);
	//	return 4;
	//}
	//if (axis == 'z' && pos >= depth-1) {
	//	printf("Bad cut position: %c = %ld\n",axis,pos);
	//	return 4;
	//}
	/*
	width1 = width2 = width;
	height1 = height2 = height;
	depth1 = depth2 = depth;
	xstart2 = 0;
	ystart2 = 0;
	zstart2 = 0;

	if (axis == 'x') {
		width1 = pos+1;
		width2 = width - width1;
		xstart2 = width1;
	}
	if (axis == 'y') {
		height1 = pos+1;
		height2 = height - height1;
		ystart2 = height1;
	}
	if (axis == 'z') {
		depth1 = pos+1;
		depth2 = depth - depth1;
		zstart2 = depth1;
	}
*/	
	ImageType::Pointer im1 = ImageType::New();
	ImageType::SizeType imsize; 
	ImageType::IndexType imstart; 
	ImageType::RegionType imregion; 

	depth1 = depth;
	for (int ix=0; ix<=1; ix++) {
		if (ix == 0) {
			xstart = 0;
			width1 = width/2;
		} else {
			xstart = width/2;
			width1 = width - xstart;
		}
		for (int iy=0; iy<=1; iy++) {
			if (iy == 0) {
				ystart = 0;
				height1 = height/2;
			} else {
				ystart = height/2;
				height1 = height - ystart;
			}
			if (ix==0 && iy==0) {
				strcpy(Qfile,"x0y0.tif");
			} else if (ix==0 && iy==1) {
				strcpy(Qfile,"x0y1.tif");
			} else if (ix==1 && iy==0) {
				strcpy(Qfile,"x1y0.tif");
			} else if (ix==1 && iy==1) {
				strcpy(Qfile,"x1y1.tif");
			}
			imsize[0] = width1;
			imsize[1] = height1;
			imsize[2] = depth1;
			imstart[0] = 0;
			imstart[1] = 0;
			imstart[2] = 0;
			imregion.SetSize(imsize);
			imregion.SetIndex(imstart);
			im1->SetRegions(imregion);
			im1->Allocate();
			p1 = (unsigned char *)(im1->GetBufferPointer());
			xysize1 = width1*height1;

			for (x=0; x<width1; x++) {
				for (y=0; y<height1; y++) {
					for (z=0; z<depth1; z++) {
						V1(x,y,z) = V(x+xstart,y+ystart,z);
					}
				}
			}
			printf("Writing Q file: %s  dimensions: width, height, depth: %ld %ld %ld\n",Qfile,width1,height1,depth1);
			typedef itk::ImageFileWriter<ImageType> FileWriterType;
			FileWriterType::Pointer writer = FileWriterType::New();
			writer->SetInput(im1);
			writer->SetFileName(Qfile);
			typedef  itk::TIFFImageIO TIFFIOType;
			TIFFIOType::Pointer tiffIO = TIFFIOType::New();
			writer->SetImageIO(tiffIO);
			if (use_compression) {
				writer->UseCompressionOn();
				tiffIO->SetCompressionToDeflate();
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
			printf("Wrote Q file: %d %d\n",ix,iy);
		}
	}
	return 0;
}