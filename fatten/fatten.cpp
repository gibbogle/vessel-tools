/*
 * To fatten a skeleton .tif 
 */

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

typedef itk::Image<unsigned char,3> ImageType_u8; 
ImageType_u8::Pointer im_in, im_out;
long long width, height, depth, xysize;
unsigned char *p_in, *p_out;

#define V_in(a,b,c)   p_in[(c)*xysize+(b)*width+(a)]
#define V_out(a,b,c)  p_out[(c)*xysize+(b)*width+(a)]

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

int main(int argc, char**argv)
{
	int linewidth, radius;
	float r2, d2;
	FILE *fp;
	char errfile[] = "error.log";
	fp = fopen(errfile,"w");

	if (argc != 4) {
		printf("Usage: fatten input_tiff output_tiff linewidth(odd)\n");
		fprintf(fp,"Usage: compress input_tiff output_tiff linewidth(odd)\n");
		fprintf(fp,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fp,"argv: %d: %s\n",i,argv[i]);
		}		fclose(fp);
		return 1;	// Wrong command line
	}

	printf("Input image file: %s\n",argv[1]);
	printf("Output image file: %s\n",argv[2]);
	sscanf(argv[3],"%d",&linewidth);
	radius = (linewidth-1)/2;
	r2 = radius*radius;
	printf("Line width: %d  radius: %d\n",linewidth,radius);


	typedef itk::ImageFileReader<ImageType_u8> FileReaderType;
	FileReaderType::Pointer reader = FileReaderType::New();

	reader->SetFileName(argv[1]);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		fprintf(fp,"Read error on input file: %s\n",argv[1]);
		fclose(fp);
		return 2;	// Read error on input file
	}

	im_in = reader->GetOutput();

	width = im_in->GetLargestPossibleRegion().GetSize()[0];
	height = im_in->GetLargestPossibleRegion().GetSize()[1];
	depth = im_in->GetLargestPossibleRegion().GetSize()[2];
	xysize = width*height;
	p_in = (unsigned char *)(im_in->GetBufferPointer());

	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	ImageType_u8::Pointer im_out = ImageType_u8::New();
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
	im_out->SetRegions(imregion);
	im_out->Allocate();
	p_out = (unsigned char *)(im_out->GetBufferPointer());

	printf("Created an output file buffer\n");

	// Now fatten the lines
	for (int x=0; x<width; x++) {
//		printf("x: %d\n",x);
		for (int y=0; y<height; y++) {
			for (int z=0; z<depth; z++) {
//				printf("x,y,z: %d %d %d\n",x,y,z);
				int indx = (z)*xysize+(y)*width+(x);
//				printf("index: %d  p_in[0]: %d  V_in: %d\n",indx,p_in[0],V_in(x,y,z));
				if (V_in(x,y,z) != 0) {
					int xfrom = MAX(0,x-radius);
					int xto = MIN(width-1,x+radius);
					int yfrom = MAX(0,y-radius);
					int yto = MIN(height-1,y+radius);
					int zfrom = MAX(0,z-radius);
					int zto = MIN(depth-1,z+radius);
//					printf("V_in: %d  %d %d %d %d %d %d\n",V_in(x,y,z),xfrom,xto,yfrom,yto,zfrom,zto);
					for (int xx = xfrom; xx <= xto; xx++) {
						for (int yy = yfrom; yy <= yto; yy++) {
							for (int zz = zfrom; zz <= zto; zz++) {
//								d2 = (x-xx)*(x-xx) + (y-yy)*(y-yy) + (z-zz)*(z-zz);
								V_out(xx,yy,zz) = 255;
							}
						}
					}
				}
			}
		}
	}

	typedef itk::ImageFileWriter<ImageType_u8> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(argv[2]);
	writer->SetInput(im_out);

	typedef  itk::TIFFImageIO TIFFIOType;
	TIFFIOType::Pointer tiffIO = TIFFIOType::New();
	writer->SetImageIO(tiffIO);
	writer->UseCompressionOn();
	//if (method == 'P') {
	//	tiffIO->SetCompressionToPackBits();
	//	printf("Compressing with Packbits\n");
	//} else if (method == 'L') {
	//	tiffIO->SetCompressionToLZW();
	//	printf("Compressing with LZW\n");
	//} else if (method == 'Z') {
	//	tiffIO->SetCompressionToDeflate();
	//	printf("Compressing with Zip/Deflate\n");
	//}
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		fprintf(fp,"Write error on output file\n");
		fclose(fp);
		return 3;	// Write error on output file
	}
	printf("Created compressed image file: %s\n",argv[2]);

	return 0;
}