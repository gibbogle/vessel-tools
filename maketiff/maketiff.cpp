/*
 * Make a tiff file from data (direction = 1)
 * This direction does not convert voxel values - if the data is (0,1) the tiff will be (0,1)
 * or
 * Make a data file from a tiff file (direction = 0)
 * This direction is currently hard-wired to generate a compressed (0,1) data file from a binary (0,255) tiff.
 * This is the format of a close (hull) file.
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

typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im;
unsigned char *p, *pc;

#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)]

int main(int argc, char**argv)
{
	int nx, ny, nz, nbytes, direction, nx8, ny8, nz8, k;
//	int n, x, y, z, dx, dy;
	long long width, height, depth, xysize;
	char *dataFile, *tiffFile;
	bool use_compression = true, data2tiff, compressdata = true;
	FILE *fpdata;

	if (argc != 4) {
		printf("Usage: maketiff dataFile tiffFile direction\n");
		printf("       direction = 0 (tiff to data) = 1 (data to tiff)\n");
		return 1;
	}

	dataFile = argv[1];
	tiffFile = argv[2];
	sscanf(argv[3],"%d",&direction);
	if (direction == 0) {
		data2tiff = false;
		printf("Output data file: %s\n",dataFile);
	} else {
		data2tiff = true;
		printf("Output image file: %s\n",tiffFile);
	}

	if (data2tiff) {		// create tiff file
		fpdata = fopen(dataFile,"rb");
		if (fpdata == NULL) {
			return 2;
		}
		fread(&nx,4,1,fpdata);
		fread(&ny,4,1,fpdata);
		fread(&nz,4,1,fpdata);
//		fread(&nbytes,4,1,fpdata);

		//if (nbytes != nx*ny*nz) {
		//	printf("Error: this is a compressed data file\n");
		//	return 10;
		//}
		width = nx;
		height = ny;
		depth = nz;
		xysize = width*height;
		printf("Desired image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

		ImageType::Pointer im = ImageType::New();
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
		im->SetRegions(imregion);
		im->Allocate();
		p = (unsigned char *)(im->GetBufferPointer());

		nbytes = nx*ny*nz;
		fread(p,nbytes,1,fpdata);

		typedef itk::ImageFileWriter<ImageType> FileWriterType;
		FileWriterType::Pointer writer = FileWriterType::New();
		writer->SetFileName(tiffFile);
		writer->SetInput(im);
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
		if (use_compression) {
			printf("Created compressed image file: %s\n",tiffFile);
		} else {
			printf("Created uncompressed image file: %s\n",tiffFile);
		}
	} else {		// create data file
		typedef itk::ImageFileReader<ImageType> FileReaderType;
		FileReaderType::Pointer reader = FileReaderType::New();

		reader->SetFileName(tiffFile);
		try
		{
			reader->Update();
		}
		catch (itk::ExceptionObject &e)
		{
			std::cout << e << std::endl;
			printf("Read error on input file\n");
			return 2;	// Read error on input tiff file
		}

		im = reader->GetOutput();
		width = im->GetLargestPossibleRegion().GetSize()[0];
		height = im->GetLargestPossibleRegion().GetSize()[1];
		depth = im->GetLargestPossibleRegion().GetSize()[2];
		printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

		p = (unsigned char *)(im->GetBufferPointer());

		fpdata = fopen(dataFile,"wb");
		if (fpdata == NULL) {
			return 2;
		}
		nx = width;
		ny = height;
		nz = depth;

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

		width = nx;
		height = ny;
		depth = nz;
		xysize = width*height;
// write(nfdata) nx,ny,nz,(((imagedata(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
		// try this...
		if (compressdata) {
			fwrite(pc,1,nbytes,fpdata);
		} else {
			fwrite(p,1,nx*ny*nz,fpdata);
		}
		fclose(fpdata);
	}

	return 0;
}
