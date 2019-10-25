//------------------------------------------------------------------------------------------------
// Approximate flood-filling of vessels is achieved by filling each 2D slice.
// This entails identifying all connected objects and removing those less than a specified size.
//------------------------------------------------------------------------------------------------
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkSize.h"
#include "itkInvertIntensityImageFilter.h"

#define V3D(a,b,c)  p3D[(c)*imsize_xy+(b)*width+(a)]
#define V_xy(a,b)  p[(b)*width+(a)]
#define V_xz(a,b)  p[(b)*width+(a)]
#define V_yz(a,b)  p[(b)*height+(a)]
#define label_xy(a,b)  plabel[(b)*width+(a)]
#define label_xz(a,b)  plabel[(b)*width+(a)]
#define label_yz(a,b)  plabel[(b)*height+(a)]

short *p, *plabel;
unsigned char *p3D;
unsigned int *npixels;
bool *fill;
long long width, height, depth, imsize_xy;
int nobjects;
int niter;
double gapwidth;
double dx_voxel, dy_voxel, dz_voxel;

//typedef itk::Image<unsigned char,2> ImageType2D;
typedef itk::Image<short,2> ImageType2D;
typedef itk::Image<unsigned char,3> ImageType3D;

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
int CountPixels(char axis)
{
	int x, y, z;
	short v, lab;

	for (int k=1; k <= nobjects; k++)
		npixels[k] = 0;
	if (axis == 'Z') {
		for (x=0; x<width; x++) {
			for (y=0; y<height; y++) {
				v = V_xy(x,y);
				lab = label_xy(x,y);
				if (lab > 0) {
					npixels[lab]++;
				}
			}
		}
	} else if (axis == 'Y') {
		for (x=0; x<width; x++) {
			for (z=0; z<depth; z++) {
				v = V_xz(x,z);
				lab = label_xz(x,z);
				if (lab > 0) {
					npixels[lab]++;
				}
			}
		}
	} else if (axis == 'X') {
		for (y=0; y<height; y++) {
			for (z=0; z<depth; z++) {
				v = V_yz(y,z);
				lab = label_yz(y,z);
				if (lab > 0) {
					npixels[lab]++;
				}
			}
		}
	}
	return 0;
}

//------------------------------------------------------------------------------------------------
// Need to exclude apparent objects that touch the boundary!
//------------------------------------------------------------------------------------------------
int FillPixels(char axis)
{
	int x, y, z;
//	unsigned char v, lab;
	short v, lab;

	if (axis == 'Z') {
		for (x=0; x<width; x++) {
			lab = label_xy(x,0);
			if (lab > 0) fill[lab] = false;
			lab = label_xy(x,height-1);
			if (lab > 0) fill[lab] = false;
		}
		for (y=0; y<height; y++) {
			lab = label_xy(0,y);
			if (lab > 0) fill[lab] = false;
			lab = label_xy(width-1,y);
			if (lab > 0) fill[lab] = false;
		}
		for (x=0; x<width; x++) {
			for (y=0; y<height; y++) {
				v = V_xy(x,y);
				lab = label_xy(x,y);
				if (lab > 0) {
					if (fill[lab]) {
						V_xy(x,y) = 0;
					}
				}
			}
		}
	} else if (axis == 'Y') {
		for (x=0; x<width; x++) {
			lab = label_xz(x,0);
			if (lab > 0) fill[lab] = false;
			lab = label_xz(x,depth-1);
			if (lab > 0) fill[lab] = false;
		}
		for (z=0; z<depth; z++) {
			lab = label_xz(0,z);
			if (lab > 0) fill[lab] = false;
			lab = label_xz(width-1,z);
			if (lab > 0) fill[lab] = false;
		}
		for (x=0; x<width; x++) {
			for (z=0; z<depth; z++) {
				v = V_xz(x,z);
				lab = label_xz(x,z);
				if (lab > 0) {
					if (fill[lab]) {
						V_xz(x,z) = 0;
					}
				}
			}
		}
	} else if (axis == 'X') {
		// First check for boundary touch
		for (y=0; y<height; y++) {
			lab = label_yz(y,0);
			if (lab > 0) fill[lab] = false;
			lab = label_yz(y,depth-1);
			if (lab > 0) fill[lab] = false;
		}
		for (z=0; z<depth; z++) {
			lab = label_yz(0,z);
			if (lab > 0) fill[lab] = false;
			lab = label_yz(height-1,z);
			if (lab > 0) fill[lab] = false;
		}
		for (y=0; y<height; y++) {
			for (z=0; z<depth; z++) {
				v = V_yz(y,z);
				lab = label_yz(y,z);
				if (lab > 0) {
					if (fill[lab]) {
						V_yz(y,z) = 0;
					}
				}
			}
		}
	}
	return 0;
}

//------------------------------------------------------------------------------------------------
// Look for pixels that are in gaps in an object.  A pixel P is in a gap if a line through the pixel
// in at least one of the 4 principal directions intersects with the same object at two points P1
// and P2 on opposite sides of P and the distance P1-P2 <= gapwidth.
//------------------------------------------------------------------------------------------------
int CloseGaps_xy(void)
{
	int x, y, x0, y0;
	int kdir, j, sgn, kstep, P1[2], P2[2], kdirmin, labmin;
	double  dx, dy, d2, w2min;
	int dir[4][2];
	short v, lab, hitlab[2];

	dir[0][0] = 1; dir[0][1] = 0;
	dir[1][0] = 1; dir[1][1] = 1;
	dir[2][0] = 0; dir[2][1] = 1;
	dir[3][0] = 1; dir[3][1] = -1;
	for (x0=0; x0<width; x0++) {
		for (y0=0; y0<height; y0++) {
			v = V_xy(x0,y0);
			if (v != 0) continue;
			w2min = 999;
			for (kdir=0; kdir<4; kdir++) {
				hitlab[0] = hitlab[1] = 0;
				for (j=0; j<2; j++) {
					sgn = 2*j-1;
					for (kstep=1; kstep<gapwidth; kstep++) {
						x = x0 + sgn*kstep*dir[kdir][0];
						if (x < 0 || x > width-1) continue;
						y = y0 + sgn*kstep*dir[kdir][1];
						if (y < 0 || y > height-1) continue;
						lab = label_xy(x,y);
						if (lab > 0) {
							hitlab[j] = lab;
							if (sgn > 0) {
								P1[0] = x; P1[1] = y;
							} else {
								P2[0] = x; P2[1] = y;
							}
							break;
						}
					}
				}
				if (hitlab[0] > 0 && hitlab[0] == hitlab[1]) {
					dx = P1[0] - P2[0];
					dy = P1[1] - P2[1];
					dx *= dx_voxel;
					dy *= dy_voxel;
					d2 = dx*dx + dy*dy;
					if (d2 < w2min) {
						w2min = d2;
						kdirmin = kdir;
						labmin = hitlab[0];
					}
				}
			}
			if (w2min <= gapwidth*gapwidth) {
				V_xy(x0,y0) = 255;
				label_xy(x0,y0) = labmin;
				for (j=0; j<2; j++) {
					sgn = 2*j-1;
					for (kstep=1; kstep<gapwidth; kstep++) {
						x = x0 + sgn*kstep*dir[kdirmin][0];
						y = y0 + sgn*kstep*dir[kdirmin][1];
						lab = label_xy(x,y);
						if (lab > 0) break;
						V_xy(x,y) = 255;
						label_xy(x,y) = labmin;
					}
				}
			}
		}
	}
	return 0;
}

int CloseGaps_xz(void)
{
	int x, z, x0, z0;
	int kdir, j, sgn, kstep, P1[2], P2[2], kdirmin, labmin;
	double  dx, dz, d2, w2min;
	int dir[4][2];
	short v, lab, hitlab[2];

	dir[0][0] = 1; dir[0][1] = 0;
	dir[1][0] = 1; dir[1][1] = 1;
	dir[2][0] = 0; dir[2][1] = 1;
	dir[3][0] = 1; dir[3][1] = -1;
	for (x0=0; x0<width; x0++) {
		for (z0=0; z0<depth; z0++) {
			v = V_xz(x0,z0);
			if (v != 0) continue;
			w2min = 999;
			for (kdir=0; kdir<4; kdir++) {
				hitlab[0] = hitlab[1] = 0;
				for (j=0; j<2; j++) {
					sgn = 2*j-1;
					for (kstep=1; kstep<gapwidth; kstep++) {
						x = x0 + sgn*kstep*dir[kdir][0];
						if (x < 0 || x > width-1) continue;
						z = z0 + sgn*kstep*dir[kdir][1];
						if (z < 0 || z > depth-1) continue;
						lab = label_xz(x,z);
						if (lab > 0) {
							hitlab[j] = lab;
							if (sgn > 0) {
								P1[0] = x; P1[1] = z;
							} else {
								P2[0] = x; P2[1] = z;
							}
							break;
						}
					}
				}
				if (hitlab[0] > 0 && hitlab[0] == hitlab[1]) {
					dx = P1[0] - P2[0];
					dz = P1[1] - P2[1];
					dx *= dx_voxel;
					dz *= dz_voxel;
					d2 = dx*dx + dz*dz;
					if (d2 < w2min) {
						w2min = d2;
						kdirmin = kdir;
						labmin = hitlab[0];
					}
				}
			}
			if (w2min <= gapwidth*gapwidth) {
				V_xz(x0,z0) = 255;
				label_xz(x0,z0) = labmin;
				for (j=0; j<2; j++) {
					sgn = 2*j-1;
					for (kstep=1; kstep<gapwidth; kstep++) {
						x = x0 + sgn*kstep*dir[kdirmin][0];
						z = z0 + sgn*kstep*dir[kdirmin][1];
						lab = label_xz(x,z);
						if (lab > 0) break;
						V_xz(x,z) = 255;
						label_xz(x,z) = labmin;
					}
				}
			}
		}
	}
	return 0;
}

int CloseGaps_yz(void)
{
	int y, z, y0, z0;
	int kdir, j, sgn, kstep, P1[2], P2[2], kdirmin, labmin;
	double  dz, dy, d2, w2min;
	int dir[4][2];
	short v, lab, hitlab[2];

//	printf("Closing gaps\n");
	dir[0][0] = 1; dir[0][1] = 0;
	dir[1][0] = 1; dir[1][1] = 1;
	dir[2][0] = 0; dir[2][1] = 1;
	dir[3][0] = 1; dir[3][1] = -1;
	for (y0=0; y0<height; y0++) {
		for (z0=0; z0<depth; z0++) {
			v = V_yz(y0,z0);
			if (v != 0) continue;
			w2min = 999;
			for (kdir=0; kdir<4; kdir++) {
				hitlab[0] = hitlab[1] = 0;
				for (j=0; j<2; j++) {
					sgn = 2*j-1;
					for (kstep=1; kstep<gapwidth; kstep++) {
						y = y0 + sgn*kstep*dir[kdir][0];
						if (y < 0 || y > height-1) continue;
						z = z0 + sgn*kstep*dir[kdir][1];
						if (z < 0 || z > depth-1) continue;
						lab = label_yz(y,z);
						if (lab > 0) {
							hitlab[j] = lab;
							if (sgn > 0) {
								P1[0] = y; P1[1] = z;
							} else {
								P2[0] = y; P2[1] = z;
							}
							break;
						}
					}
				}
				if (hitlab[0] > 0 && hitlab[0] == hitlab[1]) {
					dy = P1[0] - P2[0];
					dz = P1[1] - P2[1];
					dy *= dy_voxel;
					dz *= dz_voxel;
					d2 = dy*dy + dz*dz;
					if (d2 < w2min) {
						w2min = d2;
						kdirmin = kdir;
						labmin = hitlab[0];
					}
				}
			}
			if (w2min <= gapwidth*gapwidth) {
				V_yz(y0,z0) = 255;
				label_yz(y0,z0) = labmin;
				for (j=0; j<2; j++) {
					sgn = 2*j-1;
					for (kstep=1; kstep<gapwidth; kstep++) {
						y = y0 + sgn*kstep*dir[kdirmin][0];
						z = z0 + sgn*kstep*dir[kdirmin][1];
						lab = label_yz(y,z);
						if (lab > 0) break;
						V_yz(y,z) = 255;
						label_yz(y,z) = labmin;
					}
				}
			}
		}
	}
	return 0;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
int WriteImage(ImageType2D::Pointer image, const char filename[])
{
	typedef itk::ImageFileWriter<ImageType2D> FileWriterType;

	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(image);
//	writer->UseCompressionOn();
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
//		fprintf(fp,"Write error on output file\n");
//		fclose(fp);
		return 3;	// Write error on output file
	}
	return 0;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
int GetSlice(char axis, int xyz)
{
	int x, y, z;

	if (axis == 'Z') {
		z = xyz;
		for (x=0; x<width; x++) {
			for (y=0; y<height; y++) {
				V_xy(x,y) = V3D(x,y,z);
			}
		}
	} else if (axis == 'Y') {
		y = xyz;
		for (x=0; x<width; x++) {
			for (z=0; z<depth; z++) {
				V_xz(x,z) = V3D(x,y,z);
			}
		}
	} else if (axis == 'X') {
		x = xyz;
		for (y=0; y<height; y++) {
			for (z=0; z<depth; z++) {
				V_yz(y,z) = V3D(x,y,z);
			}
		}
	}
	return 0;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
int PutSlice(char axis, int xyz)
{
	int x, y, z;

	if (axis == 'Z') {
		z = xyz;
		for (x=0; x<width; x++) {
			for (y=0; y<height; y++) {
				V3D(x,y,z) = V_xy(x,y);
			}
		}
	} else if (axis == 'Y') {
		y = xyz;
		for (x=0; x<width; x++) {
			for (z=0; z<depth; z++) {
				V3D(x,y,z) = V_xz(x,z);
			}
		}
	} else if (axis == 'X') {
		x = xyz;
		for (y=0; y<height; y++) {
			for (z=0; z<depth; z++) {
				V3D(x,y,z) = V_yz(y,z);
			}
		}
	}
	return 0;
}


//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
int main( int argc, char *argv[])
{
	char *infile3D, *outfile3D;
	unsigned int limit = 4000;
	int x, y, z, nfilled;
	double pixel_area;
	char axis;
	ImageType3D::Pointer image3D;
	typedef itk::ImageFileReader<ImageType3D> FileReaderType;
	FileReaderType::Pointer reader = FileReaderType::New();
	typedef itk::InvertIntensityImageFilter <ImageType2D> InvertIntensityImageFilterType;
	typedef itk::ConnectedComponentImageFilter <ImageType2D, ImageType2D> ConnectedComponentImageFilterType;
	FILE *fp;
	char errfile[] = "error.log";

	fp = fopen(errfile,"w");

	if (argc != 9) {
		printf("Usage: flood2D input_tiff output_tiff arealimit gapwidth niter dx dy dz\n");
		printf("       where: arealimit is the maximum area of a cavity to be filled (um^2)\n");
		printf("              gapwidth is the maximum length of a gap to be closed (um)\n");
		printf("              iterations is the number of iterations\n");
		printf("              dx, dy, dz are the voxel dimensions (um)\n");
		fprintf(fp,"Usage: flood2D input_tiff output_tiff arealimit gapwidth niter dx dy dz\n");
		fprintf(fp,"       where: arealimit is the maximum area of a cavity to be filled (um^2)\n");
		fprintf(fp,"              gapwidth is the maximum length of a gap to be closed\n");
		fprintf(fp,"              iterations is the number of iterations\n");
		fprintf(fp,"              dx, dy, dz are the voxel dimensions (um)\n");
		fprintf(fp,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fp,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fp);
		return 1;
	}
	infile3D = argv[1];
	outfile3D = argv[2];
	printf("Input image file: %s\n",infile3D);
	printf("Output image file: %s\n",outfile3D);
	sscanf(argv[3],"%d",&limit);
	printf("limit: %d\n",limit);
	sscanf(argv[4],"%d",&gapwidth);
	printf("gapwidth: %d\n",gapwidth);
	sscanf(argv[5],"%d",&niter);
	printf("niter: %d\n",niter);
	sscanf(argv[6],"%lf",&dx_voxel);
	sscanf(argv[7],"%lf",&dy_voxel);
	sscanf(argv[8],"%lf",&dz_voxel);
	printf("voxel dimensions: %lf %lf %lf\n",dx_voxel,dy_voxel,dz_voxel);

	printf("Reading 3D image file\n");
	reader->SetFileName(infile3D);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
//		fprintf(fp,"Read error on input file\n");
//		fclose(fp);
		return 2;	// Read error on input file
	}

	image3D = reader->GetOutput();
	p3D = (unsigned char *)(image3D->GetBufferPointer());

	width = image3D->GetLargestPossibleRegion().GetSize()[0];
	height = image3D->GetLargestPossibleRegion().GetSize()[1];
	depth = image3D->GetLargestPossibleRegion().GetSize()[2];
	imsize_xy = width*height;
	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	/*
	// Test looping time - about 5 sec
	double sum=0;
	for (z=0; z<depth; z++) {
		printf(".");
		for (y=0; y<height; y++) {
			for (x=0; x<width; x++) {
				short v = V3D(x,y,z);
				sum += v;
			}
		}
	}
	return 0;
	*/
	// Create 2D image of desired size to hold an xy slice

	for (int iter=0;iter<niter;iter++) {

	// Do xz slices
	{
	printf("Processing xz slices\n");
	axis = 'Y';
	pixel_area = dx_voxel*dz_voxel;
	ImageType2D::Pointer image;
	ImageType2D::Pointer labelimage;
	image = ImageType2D::New();
	ImageType2D::SizeType imsize; 
	imsize[0] = width;
	imsize[1] = depth;
	ImageType2D::IndexType imstart; 
	imstart[0] = 0;
	imstart[1] = 0;
	ImageType2D::RegionType imregion; 
	imregion.SetSize(imsize);
	imregion.SetIndex(imstart);
	image->SetRegions(imregion);
	image->Allocate();

	nfilled = 0;
	for (y=0; y<height; y++)
	{
//		p = (unsigned char *)(image->GetBufferPointer());
		p = (short *)(image->GetBufferPointer());
		GetSlice(axis,y);
		printf("xz slice: %d\n",y);

		// Find connected objects
		ConnectedComponentImageFilterType::Pointer labelFilter = ConnectedComponentImageFilterType::New ();
		if (gapwidth > 0) {
			labelFilter->SetInput(image);
			labelFilter->Update();
			nobjects = labelFilter->GetObjectCount();
			labelimage = labelFilter->GetOutput();
			plabel = (short *)(labelimage->GetBufferPointer());

//		if (gapwidth > 0) {
			// Close gaps in connected objects
			CloseGaps_xz();
		}

		// Invert the image
		InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();
		invertIntensityFilter->SetInput(image);
		invertIntensityFilter->SetMaximum(255);
		invertIntensityFilter->Update();
		image = invertIntensityFilter->GetOutput();
//		p = (unsigned char *)(image->GetBufferPointer());
		p = (short *)(image->GetBufferPointer());

//		WriteImage(image,"inverted.tif");
//		printf("Created inverted image file: %s\n","inverted.tif");

		labelFilter->SetInput(image);
		labelFilter->Update();
		nobjects = labelFilter->GetObjectCount();
		labelimage = labelFilter->GetOutput();

		npixels = (unsigned int *)malloc((nobjects+1)*sizeof(unsigned int));
		fill = (bool *)malloc((nobjects+1)*sizeof(bool));

//		plabel = (unsigned char *)(labelimage->GetBufferPointer());
		plabel = (short *)(labelimage->GetBufferPointer());
		CountPixels(axis);

		for (int k=1; k<=nobjects; k++) {
			if (npixels[k]*pixel_area <= limit) {
				fill[k] = true;
				nfilled++;
			} else {
				fill[k] = false;
			}
		}

		FillPixels(axis);

		// Invert back again
		invertIntensityFilter->SetInput(image);
		invertIntensityFilter->SetMaximum(255);
		invertIntensityFilter->Update();
		image = invertIntensityFilter->GetOutput();
//		p = (unsigned char *)(image->GetBufferPointer());
		p = (short *)(image->GetBufferPointer());

		PutSlice(axis,y);
		free (npixels);
	}
	printf("Number of holes filled: %d\n\n",nfilled);
	}


	// Do yz slices
	{
	printf("Processing yz slices\n");
	axis = 'X';
	pixel_area = dy_voxel*dz_voxel;
	ImageType2D::Pointer image;
	ImageType2D::Pointer labelimage;
	image = ImageType2D::New();
	ImageType2D::SizeType imsize; 
	imsize[0] = height;
	imsize[1] = depth;
	ImageType2D::IndexType imstart; 
	imstart[0] = 0;
	imstart[1] = 0;
	ImageType2D::RegionType imregion; 
	imregion.SetSize(imsize);
	imregion.SetIndex(imstart);
	image->SetRegions(imregion);
	image->Allocate();

	nfilled = 0;
	for (x=0; x<width; x++)
	{
//		p = (unsigned char *)(image->GetBufferPointer());
		p = (short *)(image->GetBufferPointer());
		GetSlice(axis,x);
		printf("yz slice: %d\n",x);

		// Find connected objects
		ConnectedComponentImageFilterType::Pointer labelFilter = ConnectedComponentImageFilterType::New ();
		if (gapwidth > 0) {
			labelFilter->SetInput(image);
			labelFilter->Update();
			nobjects = labelFilter->GetObjectCount();
			labelimage = labelFilter->GetOutput();
			plabel = (short *)(labelimage->GetBufferPointer());

//		if (gapwidth > 0) {
			// Close gaps in connected objects
			CloseGaps_yz();
		}

		// Invert the image
		InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();
		invertIntensityFilter->SetInput(image);
		invertIntensityFilter->SetMaximum(255);
		invertIntensityFilter->Update();
		image = invertIntensityFilter->GetOutput();
//		p = (unsigned char *)(image->GetBufferPointer());
		p = (short *)(image->GetBufferPointer());

		labelFilter->SetInput(image);
		labelFilter->Update();
		nobjects = labelFilter->GetObjectCount();
		labelimage = labelFilter->GetOutput();

		npixels = (unsigned int *)malloc((nobjects+1)*sizeof(unsigned int));
		fill = (bool *)malloc((nobjects+1)*sizeof(bool));

//		plabel = (unsigned char *)(labelimage->GetBufferPointer());
		plabel = (short *)(labelimage->GetBufferPointer());
		CountPixels(axis);

		for (int k=1; k<=nobjects; k++) {
			if (npixels[k]*pixel_area <= limit) {
				fill[k] = true;
				nfilled++;
			} else {
				fill[k] = false;
			}
		}

		FillPixels(axis);

		// Invert back again
		invertIntensityFilter->SetInput(image);
		invertIntensityFilter->SetMaximum(255);
		invertIntensityFilter->Update();
		image = invertIntensityFilter->GetOutput();
//		p = (unsigned char *)(image->GetBufferPointer());
		p = (short *)(image->GetBufferPointer());

		PutSlice(axis,x);
		free (npixels);
	}
	printf("Number of holes filled: %d\n\n",nfilled);
	}

	// Do xy slices
	{
	printf("Processing xy slices\n");
	axis = 'Z';
	pixel_area = dx_voxel*dy_voxel;
	ImageType2D::Pointer image;
	ImageType2D::Pointer labelimage;
	image = ImageType2D::New();
	ImageType2D::SizeType imsize; 
	imsize[0] = width;
	imsize[1] = height;
	ImageType2D::IndexType imstart; 
	imstart[0] = 0;
	imstart[1] = 0;
	ImageType2D::RegionType imregion; 
	imregion.SetSize(imsize);
	imregion.SetIndex(imstart);
	image->SetRegions(imregion);
	image->Allocate();

	nfilled = 0;
	for (z=0; z<depth; z++)
	{
//		p = (unsigned char *)(image->GetBufferPointer());
		p = (short *)(image->GetBufferPointer());
		GetSlice(axis,z);
		printf("xy slice: %d\n",z);

		// Find connected objects
		ConnectedComponentImageFilterType::Pointer labelFilter = ConnectedComponentImageFilterType::New ();
		if (gapwidth > 0) {
			labelFilter->SetInput(image);
			labelFilter->Update();
			nobjects = labelFilter->GetObjectCount();
			labelimage = labelFilter->GetOutput();
			plabel = (short *)(labelimage->GetBufferPointer());

//		if (gapwidth > 0) {
			// Close gaps in connected objects
			CloseGaps_xy();
		}

		// Invert the image
		InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();
		invertIntensityFilter->SetInput(image);
		invertIntensityFilter->SetMaximum(255);
		invertIntensityFilter->Update();
		image = invertIntensityFilter->GetOutput();
//		p = (unsigned char *)(image->GetBufferPointer());
		p = (short *)(image->GetBufferPointer());

		labelFilter->SetInput(image);
		labelFilter->Update();
		nobjects = labelFilter->GetObjectCount();
		labelimage = labelFilter->GetOutput();

		npixels = (unsigned int *)malloc((nobjects+1)*sizeof(unsigned int));
		fill = (bool *)malloc((nobjects+1)*sizeof(bool));

//		plabel = (unsigned char *)(labelimage->GetBufferPointer());
		plabel = (short *)(labelimage->GetBufferPointer());
		CountPixels(axis);

		for (int k=1; k<=nobjects; k++) {
			if (npixels[k]*pixel_area <= limit) {
				fill[k] = true;
				nfilled++;
			} else {
				fill[k] = false;
			}
		}

		FillPixels(axis);

		// Invert back again
		invertIntensityFilter->SetInput(image);
		invertIntensityFilter->SetMaximum(255);
		invertIntensityFilter->Update();
		image = invertIntensityFilter->GetOutput();
//		p = (unsigned char *)(image->GetBufferPointer());
		p = (short *)(image->GetBufferPointer());

		PutSlice(axis,z);
		free (npixels);
	}
	printf("Number of holes filled: %d\n\n",nfilled);
	}
	}

	// As a final step, need to check that the generated image is indeed a single connected object

	typedef itk::ImageFileWriter<ImageType3D> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(outfile3D);
	writer->SetInput(image3D);
	writer->UseCompressionOn();
	printf("Writing 3D image file\n");
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
//		fprintf(fp,"Write error on output file\n");
//		fclose(fp);
		return 3;	// Write error on output file
	}
	printf("Created filled 3D image file: %s\n",outfile3D);
 
  return EXIT_SUCCESS;
}
