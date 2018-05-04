/*
 * To chop out a piece of a .tif
 * Note: now range indices are 0-based to make them conform with the coordinates in a .am file.
 * This implies a divergence from Irfanview, which uses 1-based indices.
 */

#include <cstdio>
#include <vector>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <itkDanielssonDistanceMapImageFilter.h>
#include "itkSize.h"
#include "itkThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMeanImageFilter.h"

typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im;
unsigned char *p, *p2;

#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)]
#define V2(a,b,c)  p2[(c)*xysize2+(b)*width2+(a)]

int main(int argc, char**argv)
{
	int x1, x2, y1, y2, z1, z2;
	int x0, y0, z0, dx, dy, dz;
	int x, y, z, xx, yy, zz;
	long long width, height, depth, xysize;
	long long width2, height2, depth2, xysize2;
	int R, R2;
	bool use_compression = true;
	bool use_sphere;

	if (argc != 9 && argc != 7) {
		printf("Usage: crop input_tiff output_tiff x1 x2 y1 y2 z1 z2\n");
		printf("       where the ranges (x1,x2), (y1,y2), (z1,z2) define the selected region (0-based)\n");
		printf("       (any value < 0 implies use full range for this axis)\n");
		printf("   or:\n");
		printf("Usage: crop input_tiff output_tiff x0 y0 z0 R\n");
		printf("       where (x0,y0,z0) is the selected sphere centre and R is the radius\n");
		return 1;
	}

	printf("Input image file: %s\n",argv[1]);
	printf("Output image file: %s\n",argv[2]);
	if (argc == 9) {
		use_sphere = false;
		sscanf(argv[3],"%d",&x1);
		sscanf(argv[4],"%d",&x2);
		sscanf(argv[5],"%d",&y1);
		sscanf(argv[6],"%d",&y2);
		sscanf(argv[7],"%d",&z1);
		sscanf(argv[8],"%d",&z2);
	} else {
		use_sphere = true;
		sscanf(argv[3],"%d",&x0);
		sscanf(argv[4],"%d",&y0);
		sscanf(argv[5],"%d",&z0);
		sscanf(argv[6],"%d",&R);
		R2 = R*R;
	}
//	sscanf(argv[9],"%c",&comp);
// Change from 1-based to 0-based indices
//	x1--; x2--; y1--; y2--; z1--; z2--;

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
		return 2;
	}

	im = reader->GetOutput();

	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	xysize = width*height;
	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);
	p = (unsigned char *)(im->GetBufferPointer());

	if (use_sphere) {
		width2 = 2*(R+5);
		height2 = 2*(R+5);
		depth2 = 2*(R+5);
		x1 = x0 - R - 5;
		y1 = y0 - R - 5;
		z1 = z0 - R - 5;
		x2 = width - (x0 + R + 5);
		y2 = height - (y0 + R + 5);
		z2 = depth - (z0 + R + 5);
		if (x1 < 0 || y1 < 0 || z1 < 0) return 3;
		if (x2 < 0 || y2 < 0 || z2 < 0) return 3;
		printf("R: %d width,height,depth: %d %d %d\n",R,width2,height2,depth2);
	} else {
		if (x1 < 0 || x2 < 0) {
			x1 = 0;
			x2 = width-1;
			printf("Using full range of x\n");
		}
		if (y1 < 0 || y2 < 0) {
			y1 = 0;
			y2 = height-1;
			printf("Using full range of y\n");
		}
		if (z1 < 0 || z2 < 0) {
			z1 = 0;
			z2 = depth-1;
			printf("Using full range of z\n");
		}
		width2 = x2 - x1 + 1;
		height2 = y2 - y1 + 1;
		depth2 = z2 - z1 + 1;
	}
	xysize2 = width2*height2;
	printf("Desired image dimensions: width, height, depth: %d %d %d\n",width2,height2,depth2);

	ImageType::Pointer im2 = ImageType::New();
	ImageType::SizeType imsize; 
	imsize[0] = width2;
	imsize[1] = height2;
	imsize[2] = depth2;
	ImageType::IndexType imstart; 
	imstart[0] = 0;
	imstart[1] = 0;
	imstart[2] = 0;
	ImageType::RegionType imregion; 
	imregion.SetSize(imsize);
	imregion.SetIndex(imstart);
	im2->SetRegions(imregion);
	im2->Allocate();
	p2 = (unsigned char *)(im2->GetBufferPointer());

	width2 = im2->GetLargestPossibleRegion().GetSize()[0];
	height2 = im2->GetLargestPossibleRegion().GetSize()[1];
	depth2 = im2->GetLargestPossibleRegion().GetSize()[2];
	printf("Cropped image dimensions: width, height, depth: %d %d %d\n",width2,height2,depth2);

	for (xx=0; xx<width2; xx++) {
		x = xx + x1;
		for (yy=0; yy<height2; yy++) {
			y = yy + y1;
			for (zz=0; zz<depth2; zz++)	{
				z = zz + z1;
				if (use_sphere) {
					dx = x - x0;
					dy = y - y0;
					dz = z - z0;
					if (dx*dx+dy*dy+dz*dz <= R2) {
						V2(xx,yy,zz) = V(x,y,z);
					} else {
						V2(xx,yy,zz) = 0;
					}
				} else {
					V2(xx,yy,zz) = V(x,y,z);
				}
			}
		}
	}

	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(argv[2]);
	writer->SetInput(im2);
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
		return 4;
	}
	if (use_compression) {
		printf("Created compressed image file: %s\n",argv[2]);
	} else {
		printf("Created uncompressed image file: %s\n",argv[2]);
	}

	return 0;
}