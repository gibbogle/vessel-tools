#ifndef QUANTIFY_H
#define QUANTIFY_H

#include <bitset>
#include <stdio.h>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSize.h"
#include "itkTIFFImageIO.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define PI 3.14159
#define ALPHAMAX PI/12

#define STR_LEN 128
#define BIG 1.0e6

FILE *fperr=NULL;
typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im;
long long imwidth, imheight, imdepth, imsize;
unsigned char *pbuffer;

#define V(a,b,c)  pbuffer[(c)*imsize+(b)*imwidth+(a)]

int nxc, nyc, nzc, nx8, ny8, nz8, nbytes;
double voxelsize[3];
unsigned char *closedata = NULL;

bool is_setup = false;

#endif
