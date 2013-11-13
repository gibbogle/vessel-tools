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
#include <itkBinaryCloseParaImageFilter.h>
#include <itkShiftScaleImageFilter.h>

typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer inputImage;


unsigned char *p, *p2;

#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)]
#define V2(a,b,c)  p2[(c)*xysize2+(b)*width2+(a)]

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int ImWriter(ImageType::Pointer im, char *filename)
{
	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(im);
	writer->UseCompressionOn();
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 1;
	}
	return 0;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int BinImWriter(char *dataFile, unsigned char *p, int nx, int ny, int nz, bool compressdata)
{
	int nbytes, nx8, ny8, nz8, k;
	unsigned char *pc;
	FILE *fpdata;
	fpdata = fopen(dataFile,"wb");
	if (fpdata == NULL) {
		return 2;
	}

	if (compressdata) {
		// create compressed buffer pc from p using (0,1)
		// First round up nx, ny, nz to a multiple of 8
		k = nx%8;
		if (k == 0)
			nx8 = nx;
		else
			nx8 = nx + 8-k;
		k = ny%8;
		if (k == 0)
			ny8 = ny;
		else
			ny8 = ny + 8-k;
		k = nx%8;
		if (k == 0)
			nz8 = nz;
		else
			nz8 = nz + 8-k;
		nbytes = nx8*ny8*nz8/8;

		pc = (unsigned char *)malloc(nbytes*sizeof(unsigned char));
		int kbit = 0;
		int kbyte = 0;
		k = 0;
		unsigned char pcbyte = 0;
		for (int iz=0; iz<nz8; iz++) {
			bool zok = (iz<nz);
			for (int iy=0; iy<ny8; iy++) {
				bool yok = (iy<ny);
				for (int ix=0; ix<nx8; ix++) {
					bool xok = (ix<nx);
					if (xok && yok && zok) {
						// set bit kbit of pcbyte to (0,1) based on p[k]
						if (p[k] == 1 || p[k] == 255) {
							pcbyte |= 1 << kbit;
						}
						if (kbyte < 4) printf("ix,iy,iz k p[k] pcbyte: %d %d %d  %d  %d  %d\n", ix,iy,iz,k,p[k],pcbyte);
						k++;
					}
					kbit++;
					if (kbit == 8) {
						pc[kbyte] = pcbyte;
						if (kbyte == 0) {
							printf("byte:%d  %d\n",kbyte,pc[kbyte]);
						}
						kbyte++;
						pcbyte = 0;
						kbit = 0;
					}
				}
			}
		}
	}
	fwrite(&nx,4,1,fpdata);
	fwrite(&ny,4,1,fpdata);
	fwrite(&nz,4,1,fpdata);
	fwrite(&nx8,4,1,fpdata);
	fwrite(&ny8,4,1,fpdata);
	fwrite(&nz8,4,1,fpdata);
	fwrite(&nbytes,4,1,fpdata);
	printf("nx,ny,nz  nx8,ny8,nz8: %d %d %d  %d %d %d\n",nx,ny,nz,nx8,ny8,nz8);

	if (compressdata) {
		fwrite(pc,1,nbytes,fpdata);
	} else {
		fwrite(p,1,nx*ny*nz,fpdata);
	}
	fclose(fpdata);
	return 0;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int main(int argc, char**argv)
{
	int width, height, depth, xysize, compress, err;
	float threshold, closingsize;
	char *inputFile, *closeFile, *binFile;
	bool compressdata = true;

	if (argc != 7) {
		printf("Usage: close input_tiff close_tiff data_file threshold closingsize compress_data\n");
		return 0;
	}
	inputFile = argv[1];
	closeFile = argv[2];
	binFile = argv[3];
	sscanf(argv[4],"%f",&threshold);
	sscanf(argv[5],"%f",&closingsize);
	sscanf(argv[6],"%d",&compress);
	compressdata = (compress == 1);
	printf("Input image file: %s\n",inputFile);
	printf("Close image file: %s\n",closeFile);
	printf("Close binary data file: %s\n",binFile);
	printf("threshold: %f\n",threshold);
	printf("closingsize: %f\n",closingsize);
	printf("compressdata: %d\n",compressdata);

	typedef itk::ImageFileReader<ImageType> FileReaderType;
	FileReaderType::Pointer reader = FileReaderType::New();

	reader->SetFileName(inputFile);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 1;
	}

	inputImage = reader->GetOutput();

	width = inputImage->GetLargestPossibleRegion().GetSize()[0];
	height = inputImage->GetLargestPossibleRegion().GetSize()[1];
	depth = inputImage->GetLargestPossibleRegion().GetSize()[2];
	xysize = width*height;
	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholdFilterType;
	ThresholdFilterType::Pointer Thresh = ThresholdFilterType::New();
	if (threshold > 0) {
		printf("Thresholding\n");
		Thresh->SetInput(inputImage);
		Thresh->SetUpperThreshold(threshold);
		Thresh->SetOutsideValue(1);	
		Thresh->SetInsideValue(0);
		Thresh->ReleaseDataFlagOn();
		//err = ImWriter(Thresh->GetOutput(),"thresh.tif");
		//if (err != 0) {
		//	printf("ImWriter error on threshold file\n");
		//	return 1;
		//}
	}

	printf("Closing\n");
	typedef itk::BinaryCloseParaImageFilter<ImageType, ImageType>  BinCloseFilterType;
	BinCloseFilterType::Pointer BinClose = BinCloseFilterType::New();
	if (threshold > 0) {
		BinClose->SetInput(Thresh->GetOutput());
	} else {
		BinClose->SetInput(inputImage);
	}
	BinClose->SetUseImageSpacing(true);
	BinClose->SetRadius(closingsize);
	BinClose->ReleaseDataFlagOn();
	err = ImWriter(BinClose->GetOutput(),"close1.tif");
	//if (err != 0) {
	//	printf("ImWriter error on close1 file\n");
	//	return 1;
	//}

	printf("Scaling\n");
    typedef itk::ShiftScaleImageFilter<ImageType, ImageType>  ScaleFilterType;
	ScaleFilterType::Pointer Scale = ScaleFilterType::New();
	Scale->SetInput(BinClose->GetOutput());
	Scale->SetScale(255.0);

	printf("Writing close tiff\n");
	err = ImWriter(Scale->GetOutput(),closeFile);
	if (err != 0) {
		printf("ImWriter error on close file\n");
		return 1;
	}

	printf("Writing close binary data file\n");
	p = (unsigned char *)(Scale->GetOutput()->GetBufferPointer());
	err = BinImWriter(binFile,p,width,height,depth,compressdata);
	return 0;
}