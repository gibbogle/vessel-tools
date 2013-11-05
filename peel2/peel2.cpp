// Segmentation of a tissue section stack that is imaged using a stain
// that shows up blood vessels and faint tissue edges. The faint
// tissue edges aren't interesting and need to be got rid of.
//
// The idea is to create a tissue mask, erode it a bit and then mask
// the original.
//
// Problem from Gib Bogle, via the ITK list.

#include <iostream>
#include "tclap/CmdLine.h"
#include "ioutils.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkBinaryShapeOpeningImageFilter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include <itkMorphologicalWatershedFromMarkersImageFilter.h>
#include <itkBinaryCloseParaImageFilter.h>
#include <itkBinaryErodeParaImageFilter.h>
#include <itkBinaryDilateParaImageFilter.h>
#include <itkMaximumImageFilter.h>
#include <itkMaskImageFilter.h>
#include "itkBinaryShapeKeepNObjectsImageFilter.h"

#include <itkSmartPointer.h>
namespace itk
{
    template <typename T>
    class Instance : public T::Pointer {
    public:
        Instance() : SmartPointer<T>( T::New() ) {}
    };
}

//#include "morphutils.h"

typedef class CmdLineType
{
public:
  std::string InputIm, OutputImPrefix, suffix;
  float xspacing, yspacing, zspacing;
  float closingsize;
  float sigma;
  float SizeOpening;
  float FirstThresh;
  float peel;
} CmdLineType;

void ParseCmdLine(int argc, char* argv[],
                  CmdLineType &CmdLineObj
                  )
{
  using namespace TCLAP;
  try
    {
    // Define the command line object.
    CmdLine cmd("peelTissueMask ", ' ', "0.9");

    ValueArg<std::string> inArg("","input","input image",true,"result","string");
    cmd.add( inArg );

    ValueArg<std::string> outArg("","outputprefix","output image prefix",true,"result","string");
    cmd.add( outArg );

    ValueArg<std::string> sufArg("","outputsuffix","output image suffix",false,".tif","string");
    cmd.add( sufArg );

    ValueArg<float> xArg("","xspacing","input image xspacing",false, 1.0 ,"float");
    cmd.add( xArg );

    ValueArg<float> yArg("","yspacing","input image xspacing",false, 1.0 ,"float");
    cmd.add( yArg );

    ValueArg<float> zArg("","zspacing","input image xspacing",false, 1.0 ,"float");
    cmd.add( zArg );

    ValueArg<float> sigArg("","sigma","smoothing sigma - can be zero",false, 0.0 ,"float");
    cmd.add( sigArg );

    ValueArg<float> threshArg("","thresh","initial threshold",false, 1.0 ,"float");
    cmd.add( threshArg );

    ValueArg<float> clArg("","closingsize","size of closing filter for marker selection",false, 1.0 ,"string");
    cmd.add( clArg );

    ValueArg<float> pArg("","peelsize","size of erosion applied to mask",false, 1.0 ,"string");
    cmd.add( pArg );

    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputIm = inArg.getValue();
    CmdLineObj.OutputImPrefix = outArg.getValue();
    CmdLineObj.suffix = sufArg.getValue();
    CmdLineObj.xspacing = xArg.getValue();
    CmdLineObj.yspacing = yArg.getValue();
    CmdLineObj.zspacing = zArg.getValue();
    CmdLineObj.sigma = sigArg.getValue();
    CmdLineObj.closingsize = clArg.getValue();
    CmdLineObj.SizeOpening = 1000;
    CmdLineObj.FirstThresh = threshArg.getValue();
    CmdLineObj.peel = pArg.getValue();
    }
  catch (ArgException &e)  // catch any exceptions
    {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

template <class PixType, int dim>
void doSeg(const CmdLineType &CmdLineObj)
{
	typedef typename itk::Image<PixType, dim> RawImType;
	typedef typename RawImType::Pointer PRawImType;

	typedef typename itk::Image<unsigned char, dim> MaskImType;
	typedef typename MaskImType::Pointer PMaskImType;

	typedef typename itk::Image<unsigned int, dim> LabImType;
	typedef typename LabImType::Pointer PLabImType;

	// load input
	PRawImType input = readIm<RawImType>(CmdLineObj.InputIm);
	// stuff to do with spacing information should go here - deal with
	// this later.

	// begin filtering to produce a foreground and background marker -
	// this uses world units for filter parameters. There is nothing
	// special about this order of operations.

	// Binarise first as the data seems clean
	// filter, and it probably wouldn't make much difference. I've
	// chosen to apply the threshold estimation step to the stage which
	// will have the most non zero voxels.

	float delta = CmdLineObj.closingsize/2;
	//float delta = 10;

	itk::Instance<itk::BinaryThresholdImageFilter<MaskImType, MaskImType> > Thresh;
	Thresh->SetInput(input);
	Thresh->SetUpperThreshold(CmdLineObj.FirstThresh);
	Thresh->SetOutsideValue(1);	
	Thresh->SetInsideValue(0);
	Thresh->ReleaseDataFlagOn();	//Gib

	//	writeImComp<MaskImType>(Thresh->GetOutput(), CmdLineObj.OutputImPrefix + "_thresh" + CmdLineObj.suffix);

	itk::Instance<itk::BinaryCloseParaImageFilter<MaskImType, MaskImType> > BinClose;

	itk::Instance<itk::BinaryErodeParaImageFilter<MaskImType, MaskImType> > BinErode;

	itk::Instance<itk::BinaryDilateParaImageFilter<MaskImType, MaskImType> > BinDilate;

	BinClose->SetInput(Thresh->GetOutput());
	printf("Did Thresh\n");
	BinClose->SetUseImageSpacing(true);
	BinClose->SetRadius(CmdLineObj.closingsize);
	BinClose->ReleaseDataFlagOn();	//Gib

	// writeImComp<MaskImType>(BinClose->GetOutput(), CmdLineObj.OutputImPrefix + "_close" + CmdLineObj.suffix);	// (0,1)

    itk::Instance<itk::ShiftScaleImageFilter<MaskImType, MaskImType> > Scale255;
	Scale255->SetScale(255.0);
	Scale255->SetInput(BinClose->GetOutput());
    writeImComp<MaskImType>(Scale255->GetOutput(), CmdLineObj.OutputImPrefix + "_close255" + CmdLineObj.suffix);	// (0,255)

	BinErode->SetInput(BinClose->GetOutput());
	BinErode->SetUseImageSpacing(true);
	BinErode->SetRadius(delta);
	printf("Did BinClose\n");

	// only keep big objects
	itk::Instance<itk::BinaryShapeOpeningImageFilter<MaskImType> > KeepBig;
	KeepBig->SetInput(BinErode->GetOutput());
	printf("Did BinErode\n");
	KeepBig->SetForegroundValue(255);	//1
	KeepBig->SetLambda(CmdLineObj.SizeOpening);
	KeepBig->SetAttribute("PhysicalSize");
	// writeImComp<MaskImType>(KeepBig->GetOutput(), CmdLineObj.OutputImPrefix + "_big" + CmdLineObj.suffix);	// (0,1)

	BinDilate->SetInput(BinClose->GetOutput());
	BinDilate->SetUseImageSpacing(true);
	BinDilate->SetRadius(delta);

	itk::Instance<itk::BinaryErodeParaImageFilter<MaskImType, MaskImType> > Peeler;
	itk::Instance<itk::BinaryThresholdImageFilter<MaskImType, MaskImType> > Selector;

	// writeImComp<MaskImType>(BinDilate->GetOutput(), CmdLineObj.OutputImPrefix + "_dilate" + CmdLineObj.suffix);	// (0,1)

	// invert to create background marker
	itk::Instance<itk::BinaryThresholdImageFilter<MaskImType, MaskImType> > Invert;
	Invert->SetInput(BinDilate->GetOutput());
	printf("Did BinDilate\n");
	Invert->SetUpperThreshold(1);
	Invert->SetLowerThreshold(1);
	Invert->SetInsideValue(0);
	Invert->SetOutsideValue(2);	

	//      writeImComp<MaskImType>(Invert->GetOutput(), CmdLineObj.OutputImPrefix + "_invert" + CmdLineObj.suffix);	// (0,0)?

	// Modification (16/7/2013) to eliminate interior voids in the marker tiff
	// Note change to Comb->SetInput
	itk::Instance<itk::BinaryShapeKeepNObjectsImageFilter<MaskImType> > SizeFilter;
	itk::Instance<itk::MaximumImageFilter<MaskImType, MaskImType, MaskImType> > Comb;

	SizeFilter->SetInput(Invert->GetOutput());
	printf("Did Invert\n");
	SizeFilter->SetBackgroundValue(0);
	SizeFilter->SetForegroundValue(2);
	SizeFilter->SetNumberOfObjects(1);
	SizeFilter->SetAttribute("PhysicalSize");
	Comb->SetInput(SizeFilter->GetOutput());
//      writeImComp<MaskImType>(SizeFilter->GetOutput(), CmdLineObj.OutputImPrefix + "_size" + CmdLineObj.suffix);	// (0,0)?
	printf("Did SizeFilter\n");
	Comb->SetInput2(KeepBig->GetOutput());

	//	writeImComp<MaskImType>(Comb->GetOutput(), CmdLineObj.OutputImPrefix + "_marker" + CmdLineObj.suffix);	//(0,1,2) 1 = desired

	// now for watershed
	itk::Instance< itk::SmoothingRecursiveGaussianImageFilter <RawImType, RawImType> > Smoother;
	Smoother->SetInput(input);
	Smoother->SetSigma(CmdLineObj.sigma);

	itk::Instance< itk::MorphologicalWatershedFromMarkersImageFilter<RawImType, MaskImType> > WS;
	if (CmdLineObj.sigma > 0)
	{
		WS->SetInput(Smoother->GetOutput());
	}
	else
	{
		WS->SetInput(input);
	}

	WS->SetMarkerImage(Comb->GetOutput());	//_marker.tif
	// throw away the background label
	Selector->SetInput(WS->GetOutput());
	Selector->SetUpperThreshold(1);
	Selector->SetLowerThreshold(1);
	Selector->SetInsideValue(1);
	Selector->SetOutsideValue(0);

//	Scale255->SetInput(Selector->GetOutput());
//	writeImComp<MaskImType>(Scale255->GetOutput(), CmdLineObj.OutputImPrefix + "_selector255" + CmdLineObj.suffix);	// (0,255)

	// writeImComp<MaskImType>(Selector->GetOutput(), CmdLineObj.OutputImPrefix + "_wsseg" + CmdLineObj.suffix);

	//itk::Instance<itk::ShiftScaleImageFilter<MaskImType, MaskImType> > Scale127;
	//Scale127->SetInput(Selector->GetOutput());
	//printf("Did Selector\n");
	//Scale127->SetScale(127.0);

	writeImComp<MaskImType>(Scale255->GetOutput(), CmdLineObj.OutputImPrefix + "_wsseg255" + CmdLineObj.suffix);

	Peeler->SetInput(Selector->GetOutput());
	Peeler->SetUseImageSpacing(true);
	Peeler->SetRadius(CmdLineObj.peel);
	printf("Did BinDilate\n");

//	writeImComp<MaskImType>(Peeler->GetOutput(), CmdLineObj.OutputImPrefix + "_peeler" + CmdLineObj.suffix);

	itk::Instance<itk::MaskImageFilter<RawImType, MaskImType> > masker;
	masker->SetInput(input);
	masker->SetInput2(Peeler->GetOutput());
	printf("Do Peeler\n");

//	writeImComp<RawImType>(masker->GetOutput(), CmdLineObj.OutputImPrefix + "_masked" + CmdLineObj.suffix);

	masker->SetInput(input);
	masker->SetInput2(Selector->GetOutput());
	writeImComp<RawImType>(masker->GetOutput(), CmdLineObj.OutputImPrefix + "_nopeel" + CmdLineObj.suffix);
	printf("Did Peeler\n");
}

int main(int argc, char * argv[])
{
  //itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);

  CmdLineType CmdLineObj;
  ParseCmdLine(argc, argv, CmdLineObj);
  const int dimension = 3;

  doSeg<unsigned char, dimension>(CmdLineObj);

  return(EXIT_SUCCESS);
}
