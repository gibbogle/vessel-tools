//--------------------------------------------------------------------------------------
// To check and correct for intensity variations from frame to frame,
// created by the sample milling.
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

typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im;
unsigned char *p, *p2;

#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)]
#define V2(a,b,c)  p2[(c)*xysize2+(b)*width2+(a)]
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
int getdma(int z, int ndelta, int delta[], int end[], int depth)
{
	int i, dma;
	for (i=0; i<ndelta; i++) {
		if (z <= end[i]) break;
	}
	if (z == end[i])
		dma = delta[i];
	else {
		if (delta[i] > delta[i-1])
			dma = MIN(delta[i],delta[i-1] + z - end[i-1]);
		else
			dma = MAX(delta[i],delta[i-1] - z + end[i-1]);
	}
	dma = MIN(dma,z);
	dma = MIN(dma,depth-1-z);
	return dma;
}

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
int main(int argc, char**argv)
{
	int x1, x2, y1, y2, z1, z2;
	int x, y, z, xx, yy, zz;
	int width, height, depth, xysize;
	int width2, height2, depth2, xysize2;
	double threshold, sum, *ave, *movingave, *scale, *newave;
	double sum0, newave0;
	char mode;
	char comp = 'C';
	bool use_compression;
	bool use_average = false;
	double alpha;

	if (argc != 6) {
		printf("Usage: equalize input_tiff output_tiff threshold mode alpha\n");
		printf("       threshold = intensity level above which pixels are summed\n");
		printf("       mode = 'A' for average, 'S' for sum\n");
		printf("       alpha is fraction of scaling to use (<1)\n");
		return 0;
	}

	sscanf(argv[3],"%lf",&threshold);
	sscanf(argv[4],"%c",&mode);
	sscanf(argv[5],"%lf",&alpha);
	if (mode == 'A' || mode == 'a')
		use_average = true;
	else
		use_average = false;
	printf("Input image file: %s\n",argv[1]);
	printf("Output image file: %s\n",argv[2]);
	printf("Threshold intensity: %f\n",threshold);
	printf("Scaling alpha: %f\n",alpha);
	if (use_average)
		printf("Using the average intensity\n");
	else
		printf("Using the summed intensity\n");

	if (comp == 'C') {
		use_compression = true;
	} else {
		use_compression = false;
	}

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
		return 1;
	}

	im = reader->GetOutput();

	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	xysize = width*height;
	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);
	p = (unsigned char *)(im->GetBufferPointer());

	ave = (double *)malloc(depth*sizeof(double));
	movingave = (double *)malloc(depth*sizeof(double));
	scale = (double *)malloc(depth*sizeof(double));

	FILE *fpout = fopen("eq.out","w");
	for (z=0; z<depth; z++) {
		int n = 0;
		sum = 0;
		for (y=0; y<height; y++) {
			for (x=0; x<width; x++) {
				if (V(x,y,z) > threshold) {
					sum += V(x,y,z);
					n++;
				}
			}
		}
		if (sum == 0) {
			printf("Sum of pixel intensities > threshold = 0, slice: %d\n",z);
			return 1;
		}
		if (use_average)
			ave[z] = sum/n;
		else
			ave[z] = sum;
	}
	fprintf(fpout,"Frame sums:\n");
	for (z=0; z<depth; z++) {
		fprintf(fpout,"%4d %8.0f\n",z,ave[z]);
	}
	return 0;

	// Moving average range tailored to the data
	int ndelta = 4;
	int *delta = (int *)malloc(ndelta*sizeof(int));
	int *end = (int *)malloc(ndelta*sizeof(int));
	delta[0] = 0;
	delta[1] = 6;
	delta[2] = 4;
	delta[3] = 0;
	end[0] = 0;
	end[1] = 185;
	end[2] = depth-1;
	for (z=0; z<depth; z++) {
		int dma = getdma(z,ndelta,delta,end,depth);
		sum = 0;
		for (int k=-dma; k<=dma; k++)
			sum += ave[z+k];
		movingave[z] = sum/(2*dma+1);
		scale[z] = alpha*movingave[z]/ave[z] + (1-alpha);

//		printf("z: %3d dma: %2d movingave: %6.2f\n",z,dma,movingave[z]);
//		fprintf(fpout,"z: %3d dma: %2d movingave: %6.2f\n",z,dma,movingave[z]);
//		movingave[z] = ave[z]
//		printf("%4d %8.4f %8.4f %8.4f\n",z,ave[z],movingave[z],scale[z]);
//		fprintf(fpout,"%4d %8.4f %8.4f %8.4f\n",z,ave[z],movingave[z],scale[z]);
	}
//	fclose(fpout);
//	return 0;

	ImageType::Pointer im2 = ImageType::New();
	ImageType::SizeType imsize; 
	imsize[0] = width;
	imsize[1] = height;
	imsize[2] = depth;
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
	width2 = width;
	xysize2 = xysize;

	bool use_original_val = true;

	newave = (double *)malloc(depth*sizeof(double));
	fprintf(fpout,"   z    scale      ave  mov_ave   newave\n");
	for (z=0; z<depth; z++) {
		int n = 0;
		int n0 = 0;
		sum = 0;
		sum0 = 0;
		for (x=0; x<width; x++)
		{
			for (y=0; y<height; y++)
			{
				int val = V(x,y,z);
				V2(x,y,z) = MIN(255,int(scale[z]*val + 0.5));
				// alpha-adjusted scaling
				bool over = false;
				if (use_original_val) {
					if (val > threshold)
						over = true;
				} else {
					if (V2(x,y,z) > threshold)
						over = true;
				}
				if (over) {
					sum += V2(x,y,z);	
					n++;
				}
				/*
				// non-adjusted scaling
				int val0 = int(val*movingave[z]/ave[z] + 0.5);
//				if (val0 > threshold) {
				if (val > threshold) {
					sum0 += MIN(255,val0);
					n0++;
				}
				*/
			}
		}
		if (use_average) {
			newave[z] = sum/n;
//			newave0 = sum0/n0;
		} else {
			newave[z] = sum;
//			newave0 = sum0;
		}
		if (use_average)
			fprintf(fpout,"%4d %8.4f %8.2f %8.2f %8.2f\n",z,scale[z],ave[z],movingave[z],newave[z]);
		else
			fprintf(fpout,"%4d %8.4f %8.0f %8.0f %8.0f\n",z,scale[z],ave[z],movingave[z],newave[z]);
	}
	return 0;

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
		return 1;
	}
	fclose(fpout);

	return 0;
}