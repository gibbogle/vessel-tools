#include "itkImage.h"
#include "itkImageFileReader.h"

unsigned char *p;
//#define V(a,b,c)  p[(c)*imsize+(b)*width+(a)]

int main( int argc, char* argv[] )
{
	char *imagefile;
	int width, height, depth, count;

	if( argc != 2 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <inputImage> ";
		std::cerr << std::endl;
		return 1;
	}
	imagefile =  argv[1] ;
	printf("Image file: %s\n",imagefile);

	typedef unsigned char PixelType;
	typedef itk::Image<unsigned char,3> ImageType;
	typedef itk::ImageFileReader< ImageType >     ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(imagefile);
	ImageType::Pointer im;
	im = reader->GetOutput();
	reader->Update();
	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	p = (unsigned char *)(im->GetBufferPointer());
	count = 0;
	for (int i=0; i<width*height*depth; i++) {
		if (p[i] != 0) count++;
	}
	printf("Lit voxel count: %d\n",count);
	return 0;
}
