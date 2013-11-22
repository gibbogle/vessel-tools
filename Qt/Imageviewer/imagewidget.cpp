/*
 * File:   imagewidget.cpp
 * Author: zian fanti
 * 
 * Created on 19 de septiembre de 2011, 19:34
 */
#include "QVBoxLayout"

#include <vtkImageReader2.h>
#include <vtkImageReader2Factory.h>
#include <vtkImageLuminance.h>

#include <vtkImageFlip.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCastImageFilter.h>
//#include <itkMedianImageFilter.h>
//#include <itkGradientAnisotropicDiffusionImageFilter.h>
//#include <itkVectorGradientAnisotropicDiffusionImageFilter.h>

#include <vtkInteractorStyleImage.h>
#include <vtkObjectFactory.h>

#include "imagewidget.h"
//#include "medianFilterDialog.h"
//#include "GADFilterDialog.h"

//ImageWidget::ImageWidget(QWidget *parent) : QWidget(parent)
ImageWidget::ImageWidget(QWidget *page)
{

	this->setAttribute(Qt::WA_DeleteOnClose);

	qvtkWidget = new QVTKWidget(this);
	//    qvtkWidget->resize(640, 480);

	QVBoxLayout *layout = new QVBoxLayout;
	layout->setContentsMargins(0, 0, 0, 0);
	layout->setSpacing(0);
	layout->addWidget(qvtkWidget);
//	this->setLayout(layout);
    // Associate the layout with page (copied from myvtk.cpp)
    page->setLayout(layout);

	// Create image actor
	actor = vtkSmartPointer<vtkImageActor>::New();

    // Create a camera 
    camera = vtkSmartPointer<vtkCamera>::New();
    
	// A renderer and render window
	renderer = vtkSmartPointer<vtkRenderer>::New();
//	renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow = qvtkWidget->GetRenderWindow();   // (as in myvtk.cpp)
    renderWindow->AddRenderer(renderer);

//    qvtkWidget->SetRenderWindow(renderWindow);

}

ImageWidget::~ImageWidget()
{
	renderWindow->Finalize();
	qvtkWidget = NULL;
	itkImage2 = NULL;
	itkImage3 = NULL;
	vtkImage = NULL;
}

// Need to find the way to remove itkImage2, itkImage3, vtkImage
void ImageWidget::clearImages() {
    renderer->RemoveActor(actor);
    printf("removed actor\n");
//    if (vtkImage) vtkImage->Delete();
//    printf("deleted vtkImage\n");
    if (itkImage2) itkImage2->Initialize();
    printf("initialized itkImage2\n");
    if (itkImage3) itkImage3->Initialize();
    printf("initialized itkImage3\n");
//    itkImage2 = NULL;
//    itkImage3 = NULL;
//    vtkImage = NULL;
    printf("did not set image pointers to NULL\n");
}

void ImageWidget::openWithVTK()
{
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::currentPath());
    
	if (!fileName.isEmpty()) {
		/* 		
		// This code is currently not used because the file is read with vtkImageReader2

		// Obtain image information
		this->setImageProperties(fileName.toAscii().data(), true);

		// set itk image depending on the image type 
		// if image type is grayscale
		if (imageType.compare("scalar") == 0) {
			// read the image
			typedef itk::ImageFileReader <ImageType> ReaderType;
			ReaderType::Pointer reader = ReaderType::New();
			reader->SetFileName(fileName.toAscii().data());
			reader->Update();

			// set the image data provided by the reader
			itkImage = reader->GetOutput();

		} else {
			// if the image is RGB
			typedef itk::ImageFileReader <RGBImageType> ReaderType;
			ReaderType::Pointer reader = ReaderType::New();
			reader->SetFileName(fileName.toAscii().data());
			reader->Update();

			// set the image data provided by the reader
			rgbItkImage = reader->GetOutput();
		}
		*/
		// reads a vtkImage for display purposes
		vtkSmartPointer <vtkImageReader2Factory> readerFactory =
			vtkSmartPointer <vtkImageReader2Factory>::New();
		vtkSmartPointer <vtkImageReader2> reader =
			readerFactory->CreateImageReader2(fileName.toAscii().data());

		reader->SetFileName(fileName.toAscii().data());
		reader->Update();

		vtkImage = reader->GetOutput();

        this->isFlipped = true;
		this->displayImage(vtkImage);


		readerFactory = NULL;
		reader = NULL;

	} else {
		QErrorMessage errorMessage;
		errorMessage.showMessage("No file specified for loading");
		errorMessage.exec();
		return;
	}
}

#define V2(a,b)  p2[(b)*imageWidth+(a)]
#define V3(a,b,c)  p3[(c)*xysize+(b)*imageWidth+(a)]

void ImageWidget::getFrame(int z)
{
	unsigned char *p2, *p3;
//	int width, height, depth;
    int xysize, x, y;

//	width = itkImage3->GetLargestPossibleRegion().GetSize()[0];
//	height = itkImage3->GetLargestPossibleRegion().GetSize()[1];
//	depth = itkImage3->GetLargestPossibleRegion().GetSize()[2];
    xysize = imageWidth*imageHeight;
    printf("Image dimensions: width, height, depth: %d %d %d\n",imageWidth,imageHeight,imageDepth);
	p3 = (unsigned char *)(itkImage3->GetBufferPointer());
	// create image itkImage2 that is the selected frame of itkImage3
    if (!itkImage2) itkImage2 = ImageType2::New();
	ImageType2::SizeType imsize; 
    imsize[0] = imageWidth;
    imsize[1] = imageHeight;
	ImageType2::IndexType imstart; 
	imstart[0] = 0;
	imstart[1] = 0;
	ImageType2::RegionType imregion; 
	imregion.SetSize(imsize);
	imregion.SetIndex(imstart);
	itkImage2->SetRegions(imregion);
	itkImage2->Allocate();
	p2 = (unsigned char *)(itkImage2->GetBufferPointer());
    for (x=0; x<imageWidth; x++) {
        for (y=0; y<imageHeight; y++) {
			V2(x,y) = V3(x,y,z);
		}
	}
//    imageDepth = depth;
}

void ImageWidget::openWithITK()
{

    fileName = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::currentPath());
    
	if (!fileName.isEmpty()) {
        printf("openWithITK\n");
        // Obtain image information
        this->setImageProperties(fileName.toAscii().data(), true);

		// read the image
		typedef itk::ImageFileReader <ImageType2> ReaderType2;
		typedef itk::ImageFileReader <ImageType3> ReaderType3;
		ReaderType2::Pointer reader2;
		ReaderType3::Pointer reader3;
		// set the image data provided by the reader

		// If this is a 3D image we need to create a 2D image from a selected slice
		if (numDimensions == 3) {
            printf("read itkImage3\n");
            reader3 = ReaderType3::New();
			reader3->SetFileName(fileName.toAscii().data());
			reader3->Update();
			itkImage3 = reader3->GetOutput();
            printf("got itkImage3\n");
        } else {
            printf("read itkImage2\n");
            reader2 = ReaderType2::New();
			reader2->SetFileName(fileName.toAscii().data());
			reader2->Update();
			itkImage2 = reader2->GetOutput();
            printf("got itkImage2\n");
        }
        if (numDimensions == 3) {
            imageWidth = itkImage3->GetLargestPossibleRegion().GetSize()[0];
            imageHeight = itkImage3->GetLargestPossibleRegion().GetSize()[1];
            imageDepth = itkImage3->GetLargestPossibleRegion().GetSize()[2];
            getFrame(0);
            printf("got frame\n");
        } else {
            imageWidth = itkImage2->GetLargestPossibleRegion().GetSize()[0];
            imageHeight = itkImage2->GetLargestPossibleRegion().GetSize()[1];
            imageDepth = 1;
        }

        printf("call show2DImage\n");
        show2DImage();

		reader2 = NULL;
		reader3 = NULL;
	}
}

QString ImageWidget::getFileName() {
    return fileName;
}

QString ImageWidget::getShortFileName() {
    QStringList pieces = fileName.split("/");
    return pieces.value( pieces.length() - 1 );
}

int ImageWidget::getDepth()
{
    return imageDepth;
}

void ImageWidget::getInfo(int *numdim, int *width, int *height, int *depth) {
    *numdim = numDimensions;
    *width = imageWidth;
    *height = imageHeight;
    *depth = imageDepth;
}

unsigned char * ImageWidget::getBuffer() {
    if (numDimensions == 3)
        return (unsigned char *)(itkImage3->GetBufferPointer());
    else
        return (unsigned char *)(itkImage2->GetBufferPointer());
}

void ImageWidget::subtractImage(unsigned char * psub) {
    unsigned char *p = getBuffer();
    for (int i=0; i<imageWidth*imageHeight*imageDepth; i++) {
        p[i] = MAX(0,p[i]-psub[i]);
    }
}

void ImageWidget::addImage(unsigned char * psub) {
    unsigned char *p = getBuffer();
    int sum;
    for (int i=0; i<imageWidth*imageHeight*imageDepth; i++) {
        sum = p[i]+psub[i];
        p[i] = MIN(255,sum);
    }
}

void ImageWidget::subtract(int val)
{
    int size, i;
    unsigned char *p;

    printf("subtract %d from image\n",val);
    if (numDimensions == 3) {
        size = imageWidth*imageHeight*imageDepth;
        p = (unsigned char *)(itkImage3->GetBufferPointer());
    } else {
        size = imageWidth*imageHeight;
        p = (unsigned char *)(itkImage2->GetBufferPointer());
    }
    for (i=0; i<size; i++) {
        if (p[i] > val)
            p[i] -= val;
        else
            p[i] = 0;
    }
}

bool ImageWidget::frameOK(int z) {
    if (z >= 0 && z < imageDepth)
        return true;
    else
        return false;
}

void ImageWidget::showFrame(int z)
{
    printf("showFrame: %d\n",z);
    if (imageDepth > 1) {
        itkImage2->Initialize();
        getFrame(z);
        printf("got frame\n");
    }
    show2DImage();
}

void ImageWidget::show2DImage()
{
    // setup and connect itk with vtk
    vtkConnectorType::Pointer connector = vtkConnectorType::New();
    connector->GetExporter()->SetInput(itkImage2);
    connector->GetImporter()->Update();
    printf("did connector\n");

    // flip image in Y axis
    vtkSmartPointer<vtkImageFlip> flipYFilter = vtkSmartPointer<vtkImageFlip>::New();
    flipYFilter->SetFilteredAxis(1); // flip Y axis
    flipYFilter->SetInput(connector->GetImporter()->GetOutput());
    flipYFilter->Update();
    printf("did flipYFilter\n");

    // create vtk image
    vtkImage = vtkSmartPointer<vtkImageData>::New();
    printf("new vtkImage\n");
    vtkImage->DeepCopy(flipYFilter->GetOutput());
    printf("DeepCopy vtkImage\n");
    vtkImage->SetScalarTypeToUnsignedChar();
    printf("update vtkImage\n");
    vtkImage->Update();

    printf("nullify connector\n");
    connector = NULL;
    printf("call displayImage\n");
    this->displayImage(vtkImage);
}

void ImageWidget::saveImage(QString saveFileName) {
    if (numDimensions == 2) {
        typedef itk::ImageFileWriter<ImageType2> FileWriterType2;
        FileWriterType2::Pointer writer2 = FileWriterType2::New();
        writer2->SetFileName(saveFileName.toAscii().constData());
        writer2->SetInput(itkImage2);
        writer2->UseCompressionOn();
        try
        {
            writer2->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cout << e << std::endl;
            return;	// Write error on output file
        }
    } else {
        typedef itk::ImageFileWriter<ImageType3> FileWriterType3;
        FileWriterType3::Pointer writer3 = FileWriterType3::New();
        writer3->SetFileName(saveFileName.toAscii().constData());
        writer3->SetInput(itkImage3);
        writer3->UseCompressionOn();
        try
        {
            writer3->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cout << e << std::endl;
            return;	// Write error on output file
        }
    }
    fileName = saveFileName;
}

/*
void ImageWidget::medianFilter()
{
	// create and show the median filter dialog 
	MedianFilterDialog filterDialog(this);
    
    // if the user don't cancel the action
	if (filterDialog.exec()) {
		// get selected value from dialog
		int intensity = filterDialog.spinBox->value();

		// if the itkImage is not loaded, then vtkImage is converted to itkImage 
		if (itkImage2.IsNull()) {
			this->setITKImageFromVTK();
		}

		// setup the itk median filter
		typedef itk::MedianImageFilter<ImageType2, ImageType2> FilterType;

		FilterType::Pointer filter = FilterType::New();
		FilterType::InputSizeType radius;
		radius.Fill(intensity);

		filter->SetRadius(radius);
		filter->SetInput(itkImage2);
        filter->Update();
        
        
		// setup and connect itk with vtk, to transform the itkImage to vtkImage
		vtkConnectorType::Pointer vtkConnector = vtkConnectorType::New();
		vtkConnector->GetExporter()->SetInput(filter->GetOutput());
		vtkConnector->GetImporter()->Update();

		itkImage2 = filter->GetOutput();
		// clear previous vtkImage
		vtkImage = NULL;

		// create new vtk image
		vtkImage = vtkSmartPointer <vtkImageData>::New();
		vtkImage->Initialize();
		vtkImage->DeepCopy(vtkConnector->GetImporter()->GetOutput());
		vtkImage->Update();

        isFlipped = false;
		this->displayImage(vtkImage);

		filter = NULL;
		vtkConnector = NULL;
	}
}

void ImageWidget::gradientAnisotropicDiffusionFilter()
{
	// create and show the median filter dialog 	
	GADFilterDialog filterDialog(this);
	if (filterDialog.exec()) {

		// if the image is in grayscale
		if (imageType.compare("scalar") == 0) {
			// set up gradient anisotropic diffusion filter
			typedef itk::GradientAnisotropicDiffusionImageFilter< ImageType2, FloatImageType > FilterType;
			FilterType::Pointer filter = FilterType::New();
			filter->SetInput(itkImage2);

			filter->SetNumberOfIterations(filterDialog.iterationsSpinBox->value());
			filter->SetTimeStep(filterDialog.timeStepSpinBox->value());
			filter->SetConductanceParameter(filterDialog.conductanceSpinBox->value());
			filter->Update();			
			
			// cast the float image to scalar image in order to display
			typedef itk::CastImageFilter< FloatImageType, ImageType2 > CastFilterType;
			CastFilterType::Pointer castFilter = CastFilterType::New();
			castFilter->SetInput(filter->GetOutput());
			
			itkImage2 = castFilter->GetOutput();

			// setup and connect itk with vtk, to transform the itkImage to vtkImage
			vtkConnectorType::Pointer vtkConnector = vtkConnectorType::New();
			vtkConnector->GetExporter()->SetInput(castFilter->GetOutput());
			vtkConnector->GetImporter()->Update();

			// clear previous vtkImage
			vtkImage = NULL;

			// create new vtk image
			vtkImage = vtkSmartPointer <vtkImageData>::New();
			vtkImage->Initialize();
			vtkImage->DeepCopy(vtkConnector->GetImporter()->GetOutput());
			vtkImage->Update();

            isFlipped = false;
			this->displayImage(vtkImage);

			filter = NULL;
			vtkConnector = NULL;
		} else {
			// if the image is RGB
			typedef itk::RGBPixel< float > FloatPixelType;
			typedef itk::Image< FloatPixelType, 2 > FloatRGBImageType;
			typedef itk::VectorGradientAnisotropicDiffusionImageFilter< RGBImageType, FloatRGBImageType > FilterType;


			FilterType::Pointer filter = FilterType::New();
			filter->SetInput(rgbItkImage);

			filter->SetNumberOfIterations(filterDialog.iterationsSpinBox->value());
			filter->SetTimeStep(filterDialog.timeStepSpinBox->value());
			filter->SetConductanceParameter(filterDialog.conductanceSpinBox->value());
			filter->Update();

			typedef itk::CastImageFilter< FloatRGBImageType, RGBImageType > CastFilterType;
			CastFilterType::Pointer castFilter = CastFilterType::New();
			castFilter->SetInput(filter->GetOutput());

			rgbItkImage = castFilter->GetOutput();

			// setup and connect itk with vtk, to transform the itkImage to vtkImage
			RGBVtkConnectorType::Pointer vtkConnector = RGBVtkConnectorType::New();
			vtkConnector->GetExporter()->SetInput(castFilter->GetOutput());
			vtkConnector->GetImporter()->Update();

			// clear previous vtkImage
			vtkImage = NULL;

			// create new vtk image
			vtkImage = vtkSmartPointer <vtkImageData>::New();
			vtkImage->Initialize();
			vtkImage->DeepCopy(vtkConnector->GetImporter()->GetOutput());
			vtkImage->Update();

            isFlipped = false;
			this->displayImage(vtkImage);

			filter = NULL;
			vtkConnector = NULL;
		}
	}
}
*/

// Define interaction style
class MouseInteractorStyle4 : public vtkInteractorStyleImage
{
  public:
    static MouseInteractorStyle4* New();
    vtkTypeMacro(MouseInteractorStyle4, vtkInteractorStyleImage);

    virtual void OnLeftButtonDown() 
    {
      //std::cout << "Pressed left mouse button." << std::endl;
      // Forward events
	  int shifted = this->Interactor->GetShiftKey();
	  if (shifted != 0) vtkInteractorStyleImage::OnLeftButtonDown();
    }
};
vtkStandardNewMacro(MouseInteractorStyle4);

void ImageWidget::displayImage(vtkImageData *image)
{

    int *dim = image->GetDimensions();
    double *spacing = image->GetSpacing(); 
    double *origin = image->GetOrigin();  

    float Cx = (dim[0] * spacing[0])/2. + origin[0];
    float Cy = (dim[1] * spacing[1])/2. + origin[1];
    camera->ParallelProjectionOn();
    camera->SetFocalPoint(Cx,Cy,0);    
    camera->SetPosition(Cx,Cy,1);
    //
    //    // to flip de image
    //    camera->SetViewUp (0, 1, 0);  
    //    
    // set actor properties
    actor->SetInput(image);
	actor->InterpolateOff();    

    renderer->AddActor(actor);
    renderer->SetActiveCamera(camera);
    renderer->ResetCamera();

//	qvtkWidget->SetRenderWindow(renderWindow); // moved to constructor
//    printf("returning from displayImage\n");    // apparently not necessary to set interactor style - default OK 
//    return;

    // window interactor style for display images 
//	vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
    vtkSmartPointer<MouseInteractorStyle4> style = vtkSmartPointer<MouseInteractorStyle4>::New();
	// set interactor style to the qvtkWidget Interactor
	qvtkWidget->GetInteractor()->SetInteractorStyle(style);

    qvtkWidget->GetInteractor()->Render();
    this->update();
}

void ImageWidget::setITKImageFromVTK()
{
	itkConnectorType::Pointer itkConnector = itkConnectorType::New();

	if (imageType.compare("rgb") == 0) {
		// Must convert image to grayscale because itkVTKImageToImageFilter only accepts single channel images
		vtkSmartPointer<vtkImageLuminance> luminanceFilter =
			vtkSmartPointer<vtkImageLuminance>::New();
		luminanceFilter->SetInput(vtkImage);
		luminanceFilter->Update();

		//get itkImage from vtkImage;  vtkImageData is unsigned char and single channel            
		itkConnector->SetInput(luminanceFilter->GetOutput());
		luminanceFilter = NULL;
	} else {
		itkConnector->SetInput(vtkImage);
	}

	itkConnector->Update();

	itkImage2 = ImageType2::New();
	itkImage2->Graft(itkConnector->GetOutput());

	itkConnector = NULL;
}

void ImageWidget::setImageProperties(std::string fileName, bool verbose)
{
	// Obtain image information
	typedef itk::ImageIOBase::IOComponentType ScalarPixelType;

	itk::ImageIOBase::Pointer imageIO =
		itk::ImageIOFactory::CreateImageIO(fileName.c_str(), itk::ImageIOFactory::ReadMode);

	imageIO->SetFileName(fileName);
	imageIO->ReadImageInformation();

	pixelType = imageIO->GetComponentTypeAsString(imageIO->GetComponentType());
	numDimensions = imageIO->GetNumberOfDimensions();
	imageType = imageIO->GetPixelTypeAsString(imageIO->GetPixelType());
    if (verbose) {
		std::cout << "Pixels type: " << pixelType << std::endl;
		std::cout << "Image type: " << imageType << std::endl;
		std::cout << "Num of Dimensions: " << numDimensions << std::endl;
	}
}
