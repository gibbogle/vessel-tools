/*
 * Apply thresholding (Local, Binary or Continuous) to a tiff.
 * The output file is an ITK-compressed version (MACpackbits)
 * Smoothing (mean) is optionally applied.
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
#include <itkDanielssonDistanceMapImageFilter.h>
#include "itkSize.h"
#include "itkThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMeanImageFilter.h"

typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im, background;
long long width, height, depth, imsize;

#define V(a,b,c)  p[(c)*imsize+(b)*width+(a)]
#define B(a,b,c)  back[(c)*imsize+(b)*width+(a)]

unsigned char *p, *back;
FILE *fp;

/*
int x1=17;
int x2=21;
int yy1=1007;
int y2=1011;
int z1=202;
int z2=203;

void show_V()
{
	int x, y, z, val;
	for (z=z1; z<=z2; z++) {
		for (y=yy1; y<=y2; y++) {
			for (x=x1; x<=x2; x++) {
				val = V(x,y,z);
				fprintf(fp,"%d   ",val);
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
}

void show_B()
{
	int x, y, z, val;
	for (z=z1; z<=z2; z++) {
		for (y=yy1; y<=y2; y++) {
			for (x=x1; x<=x2; x++) {
				val = B(x,y,z);
				fprintf(fp,"%d   ",val);
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
}
*/


//------------------------------------------------------------------------------------------------------
// There are two threshold ranges, depending on the background level b.
// For b < T, the threshold = b + delta
// For b > T, the threshold = T + delta
//------------------------------------------------------------------------------------------------------
void LocalThreshold(int T, int delta)
{
	int x, y, z;
	double thresh;
	unsigned char b, v;
	p = (unsigned char *)(im->GetBufferPointer());
	back = (unsigned char *)(background->GetBufferPointer());

//	show_V();
//	show_B();

	for (z=0; z<depth; z++)
//	for (z=214; z<depth; z++)
	{
		printf(".");
		for (y=0; y<height; y++)
		{
			for (x=0; x<width; x++)
			{
				b = B(x,y,z);
				v = V(x,y,z);
				/*
				if (b >= T)
					thresh = T + delta;
				else 
					thresh = b + delta;
				*/
				if (b >= T-delta)
					thresh = T;
				else 
					thresh = b + delta;
				if (v >= thresh)
					V(x,y,z) = 255;
				else
					V(x,y,z) = 0;
			}
		}
	}
	printf("\n");
//	show_V();
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
int main(int argc, char**argv)
{
	int T, delta;
	char *infile, *backfile, *outfile;
	char errfile[] = "error.log";
	fp = fopen(errfile,"w");

	if (argc != 6) {
		printf("Usage: threshold input_tiff back_tiff output_tiff T delta\n");
		printf("       where: input_tiff is the image to be thresholded\n");
		printf("              back_tiff is the background level image (smoothed)\n");
		printf("              output_tiff is the thresholded image\n");
		printf("              T is the basic threshold value\n");
		printf("              delta is the excess-over-background criterion at low-intensities ( < T )\n");
		printf("Note: with v = voxel value, b = background value:\n");
		printf("       the threshold is b+delta when b < T-delta, T when b >= T-delta\n");
		printf("Note: if T = delta the threshold is T and the background level is ignored\n");
		fprintf(fp,"Usage: threshold input_tiff back_tiff output_tiff T delta\n");
		fprintf(fp,"       where: input_tiff is the image to be thresholded\n");
		fprintf(fp,"              back_tiff is the background level image (smoothed)\n");
		fprintf(fp,"              output_tiff is the thresholded image\n");
		fprintf(fp,"              T is the basic threshold value\n");
		fprintf(fp,"              delta is the excess-over-background criterion at low-intensities ( < T )\n");
		fprintf(fp,"Note: with v = voxel value, b = background value:\n");
		fprintf(fp,"       the threshold is b+delta when b < T-delta, T when b >= T-delta\n");
		fprintf(fp,"Note: if T = delta the threshold is T and the background level is ignored\n\n");

		fprintf(fp,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fp,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fp);
		return 1;
	}

	infile = argv[1];
	backfile = argv[2];
	outfile = argv[3];
	printf("Input image file: %s\n",infile);
	printf("Background (smoothed) image file: %s\n",backfile);
	printf("Output image file: %s\n",outfile);
	sscanf(argv[4],"%d",&T);
	printf("Applying basic threshold level T: %d\n",T);
	sscanf(argv[5],"%d",&delta);
	printf("Applying low-intensity criterion delta: %d\n",delta);

	typedef itk::ImageFileReader<ImageType> FileReaderType;
	FileReaderType::Pointer reader1 = FileReaderType::New();

	reader1->SetFileName(infile);
	try
	{
		reader1->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		fprintf(fp,"Read error on input file\n");
		fclose(fp);
		return 2;
	}

	im = reader1->GetOutput();

	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	imsize = width*height;

	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	FileReaderType::Pointer reader2 = FileReaderType::New();
	reader2->SetFileName(backfile);
	try
	{
		reader2->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		fprintf(fp,"Read error on background file\n");
		fclose(fp);
		return 3;
	}

	background = reader2->GetOutput();
	LocalThreshold(T, delta);

	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(outfile);
	writer->SetInput(im);
	writer->UseCompressionOn();
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		fprintf(fp,"Write error on output file\n");
		fclose(fp);
		return 4;
	}
	printf("Created compressed image file: %s\n",outfile);

	return 0;
}