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

#include <vtkInteractorStyleImage.h>
#include <vtkObjectFactory.h>

// For picking
#include <vtkAbstractPicker.h>
#include <vtkVectorText.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include "vtkPolyDataMapper.h"
#include <vtkCoordinate.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include "vtkDiskSource.h"

#include "imagewidget.h"

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
ImageWidget::ImageWidget(QWidget *parent, QWidget *page)
{
	this->setAttribute(Qt::WA_DeleteOnClose);

    m_parent = parent;

	qvtkWidget = new QVTKWidget(this);

	QVBoxLayout *layout = new QVBoxLayout;
	layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(0);
    layout->addWidget(qvtkWidget);
    // Associate the layout with page (copied from myvtk.cpp)
    page->setLayout(layout);

    connect(qvtkWidget,SIGNAL(mouseEvent(QMouseEvent*)),this, SLOT(tallyMark(QMouseEvent*)));

	// Create image actor
    imageActor = vtkSmartPointer<vtkImageActor>::New();

    // Create a camera 
    camera = vtkSmartPointer<vtkCamera>::New();
    
	// A renderer and render window
	renderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow = qvtkWidget->GetRenderWindow();   // (as in myvtk.cpp)
    renderWindow->AddRenderer(renderer);

    vtkSmartPointer<vtkDiskSource> diskSource = vtkSmartPointer<vtkDiskSource>::New();
    diskSource->SetOuterRadius(2);
    mapper = vtkSmartPointer<vtkPolyDataMapper2D>::New();
    mapper->SetInputConnection(diskSource->GetOutputPort());
    std::cout << "made mapper" << std::endl;

    connect(this,SIGNAL(SaveStatus(bool)),m_parent,SLOT(setLineSaveStatus(bool)));
    clearAll();

}

ImageWidget::~ImageWidget()
{
	renderWindow->Finalize();
	qvtkWidget = NULL;
	itkImage2 = NULL;
	itkImage3 = NULL;
	vtkImage = NULL;
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::clearAll()
{
    iframe = 0;
    LButtonDown = false;
    clearActors();
    if (imageActor) {
        printf("clearImages\n");
        fflush(stdout);
        clearImages();
    }
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::clearActors()
{
    vtkSmartPointer<vtkActor2D> actor;
    vtkSmartPointer<vtkTextActor> textActor;
    printf("clearActors\n");
    fflush(stdout);
    printf("Actors.size: %d\n",Actors.size());
    fflush(stdout);
    printf("vect_actor_mat.size: %d\n",vect_actor_mat.size());
    fflush(stdout);
    if (Actors.size() > 0) {
        for (int iactor=0; iactor < Actors.size(); iactor++) {
            actor = Actors.at(iactor);
            renderer->RemoveActor2D(actor);
            textActor = NumActors.at(iactor);
            renderer->RemoveActor2D(textActor);
        }
        Actors.clear();
        NumActors.clear();
    }
    if (vect_actor_mat.size() > 0) {
        vect_actor_mat.clear();
    }
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::clearImages() {
    renderer->RemoveActor(imageActor);
//    printf("removed imageActor\n");
    if (itkImage2) itkImage2->Initialize();
//    printf("initialized itkImage2\n");
    if (itkImage3) itkImage3->Initialize();
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::openWithVTK()
{
    QString imageFileName = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::currentPath());
    
    if (!imageFileName.isEmpty()) {
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
        vtkSmartPointer <vtkImageReader2Factory> readerFactory = vtkSmartPointer <vtkImageReader2Factory>::New();
        vtkSmartPointer <vtkImageReader2> reader = readerFactory->CreateImageReader2(imageFileName.toAscii().data());

        reader->SetFileName(imageFileName.toAscii().data());
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

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::getFrame(int z)
{
	unsigned char *p2, *p3;
    long long xysize, x, y;

    printf("getFrame: %d\n",z);
    fflush(stdout);
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

    vtkSmartPointer<vtkActor2D> actor;
    vtkSmartPointer<vtkTextActor> textActor;
    QMessageBox::StandardButton reply;
    vect_actor_t oldActors;
    std::vector<vtkSmartPointer<vtkTextActor>> oldNumActors;
    if (Actors.size() > 0) {
        if (!lineSaved) {
            reply = QMessageBox::question(this, "Existing line of marks has not been saved", "Save line?", QMessageBox::Yes|QMessageBox::No);
            if (reply == QMessageBox::Yes) {
                saveLine();
            }
        }
        oldActors = Actors;
        oldNumActors = NumActors;
    }
    showActors();
    iframe = z;
    while (vect_actor_mat.size() <= iframe) {
        std::vector<vtkSmartPointer<vtkActor2D>> A;
        A.clear();
        vect_actor_mat.push_back(A);
    }
    Actors = vect_actor_mat.at(iframe);
    showActors();
    NumActors.clear();
    if (Actors.size() > 0) {    // There are stored marks, which could be used
        reply = QMessageBox::question(this, "There are stored marks for this frame", "Use stored marks?", QMessageBox::Yes|QMessageBox::No);
        if (reply == QMessageBox::Yes) {
            // Delete current marks
            for (int i=0; i<oldActors.size(); i++) {
                printf("delete actor and textActor: %d\n",i);
                actor = oldActors.at(i);
                renderer->RemoveActor2D(actor);
                textActor = oldNumActors.at(i);
                renderer->RemoveActor2D(textActor);
            }
            printf("removed old actors from renderer\n");
            fflush(stdout);
            for (int i=0; i<Actors.size(); i++) {
                printf("add actor: %d\n",i);
                fflush(stdout);
                actor = Actors.at(i);
                renderer->AddActor2D(actor);
            }
            for (int i=0; i<Actors.size(); i++) {
                printf("add textActor: %d\n",i);
                fflush(stdout);
                textActor = vtkSmartPointer<vtkTextActor>::New();
                NumActors.push_back(textActor);
                setNumber(i);
                renderer->AddActor2D(textActor);
            }
            printf("added new actors to renderer\n");
            fflush(stdout);
        }
    }

    if (iframe > 0) {
        printf("iframe: %d actors list size: new: %d old: %d\n",iframe,Actors.size(),oldActors.size());
        fflush(stdout);
        if (Actors.size() == 0 && oldActors.size() > 0) {
            printf("Could use oldActors\n");
            fflush(stdout);
//            Actors.clear();
            vtkSmartPointer<vtkActor2D> oldactor;
            vtkSmartPointer<vtkTextActor> oldtextactor;
            double *pos;
            reply = QMessageBox::question(this, "Existing marks can be used", "Copy marks?", QMessageBox::Yes|QMessageBox::No);
            if (reply == QMessageBox::Yes) {
                for (int i=0; i<oldActors.size(); i++) {
                    printf("copy actor: %d\n",i);
                    oldactor = oldActors.at(i);
                    oldtextactor = oldNumActors.at(i);
                    pos = oldactor->GetPosition();
                    actor = makeMarkActor(pos);
                    Actors.push_back(actor);
                    textActor = vtkSmartPointer<vtkTextActor>::New();
                    NumActors.push_back(textActor);
                    setNumber(i);
                    renderer->RemoveActor2D(oldactor);
                    renderer->RemoveActor2D(oldtextactor);
                    renderer->AddActor2D(actor);
                    renderer->AddActor2D(textActor);
                }
                showActors();
            } else {
                // Delete marks
                for (int i=0; i<oldActors.size(); i++) {
                    printf("delete actor and textActor: %d\n",i);
                    actor = oldActors.at(i);
                    renderer->RemoveActor2D(actor);
                    textActor = oldNumActors.at(i);
                    renderer->RemoveActor2D(textActor);
                }
            }
        }
    }
    lineSaved = false;
    fflush(stdout);
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::saveLine()
{
    printf("saveLine\n");
    fflush(stdout);
    vect_actor_mat.at(iframe) = Actors;
    lineSaved = true;
    emit(SaveStatus(false));
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::showLines()
{
    vtkSmartPointer<vtkActor2D> actor;
    printf("showLines: %d\n",vect_actor_mat.size());
    fflush(stdout);
    for (int n = 0; n < vect_actor_mat.size(); n++) {
        for (int m = 0; m < vect_actor_mat[n].size(); m++) {
            actor = vect_actor_mat[n][m];
            double *pos = actor->GetPosition();
            printf("vector: %d actor: %d pos: %f %f\n",n,m,pos[0],pos[1]);
            fflush(stdout);
        }
    }
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::showActors()
{
    vtkSmartPointer<vtkActor2D> actor;
    printf("showActors: %d\n",Actors.size());
    fflush(stdout);
    for (int m = 0; m < Actors.size(); m++) {
        actor = Actors.at(m);
        int pixel[2];
        double *pos = actor->GetPosition();
        convertPosToPixel(pos,pixel);
        printf("actor: %d pos: %f %f  pixel: %d %d\n",m,pos[0],pos[1],pixel[0],pixel[1]);
        fflush(stdout);
    }
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::openWithITK()
{
    imageFileName = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::currentPath());
    
    if (!imageFileName.isEmpty()) {
//        printf("openWithITK\n");
        // Obtain image information
        this->setImageProperties(imageFileName.toAscii().data(), true);

		// read the image
		typedef itk::ImageFileReader <ImageType2> ReaderType2;
		typedef itk::ImageFileReader <ImageType3> ReaderType3;
		ReaderType2::Pointer reader2;
		ReaderType3::Pointer reader3;
		// set the image data provided by the reader

		// If this is a 3D image we need to create a 2D image from a selected slice
		if (numDimensions == 3) {
//            printf("read itkImage3\n");
            reader3 = ReaderType3::New();
            reader3->SetFileName(imageFileName.toAscii().data());
			reader3->Update();
			itkImage3 = reader3->GetOutput();
//            printf("got itkImage3\n");
        } else {
//            printf("read itkImage2\n");
            reader2 = ReaderType2::New();
            reader2->SetFileName(imageFileName.toAscii().data());
			reader2->Update();
			itkImage2 = reader2->GetOutput();
//            printf("got itkImage2\n");
        }
        if (numDimensions == 3) {
            imageWidth = itkImage3->GetLargestPossibleRegion().GetSize()[0];
            imageHeight = itkImage3->GetLargestPossibleRegion().GetSize()[1];
            imageDepth = itkImage3->GetLargestPossibleRegion().GetSize()[2];
            getFrame(0);
            printf("imageDepth: %d\n",imageDepth);
        } else {
            imageWidth = itkImage2->GetLargestPossibleRegion().GetSize()[0];
            imageHeight = itkImage2->GetLargestPossibleRegion().GetSize()[1];
            imageDepth = 1;
            iframe = 0;
            if (vect_actor_mat.size() <= iframe) {
                std::vector<vtkSmartPointer<vtkActor2D>> A;
                vect_actor_mat.push_back(A);
            }
        }
//        printf("call show2DImage\n");
        show2DImage();

		reader2 = NULL;
		reader3 = NULL;
    }
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
QString ImageWidget::getFileName() {
    return imageFileName;
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
QString ImageWidget::getShortFileName() {
    QStringList pieces = imageFileName.split("/");
    return pieces.value( pieces.length() - 1 );
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
int ImageWidget::getDepth()
{
    return imageDepth;
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::getInfo(int *numdim, int *width, int *height, int *depth) {
    *numdim = numDimensions;
    *width = imageWidth;
    *height = imageHeight;
    *depth = imageDepth;
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
unsigned char * ImageWidget::getBuffer() {
    if (numDimensions == 3)
        return (unsigned char *)(itkImage3->GetBufferPointer());
    else
        return (unsigned char *)(itkImage2->GetBufferPointer());
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
bool ImageWidget::frameOK(int z) {
    if (z >= 0 && z < imageDepth)
        return true;
    else
        return false;
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::showFrame(int z)
{
    printf("showFrame: %d\n",z);
    if (imageDepth > 1) {
        itkImage2->Initialize();
        getFrame(z);
        printf("got frame\n");
    }
    fflush(stdout);
    show2DImage();
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::show2DImage()
{
    // setup and connect itk with vtk
    vtkConnectorType::Pointer connector = vtkConnectorType::New();
    connector->GetExporter()->SetInput(itkImage2);
    connector->GetImporter()->Update();
//    printf("did connector\n");

    // flip image in Y axis
    vtkSmartPointer<vtkImageFlip> flipYFilter = vtkSmartPointer<vtkImageFlip>::New();
    flipYFilter->SetFilteredAxis(1); // flip Y axis
    flipYFilter->SetInput(connector->GetImporter()->GetOutput());
    flipYFilter->Update();
//    printf("did flipYFilter\n");

    // create vtk image
    vtkImage = vtkSmartPointer<vtkImageData>::New();
//    printf("new vtkImage\n");
    vtkImage->DeepCopy(flipYFilter->GetOutput());
//    printf("DeepCopy vtkImage\n");
    vtkImage->SetScalarTypeToUnsignedChar();
//    printf("update vtkImage\n");
    vtkImage->Update();

//    printf("nullify connector\n");
    connector = NULL;
//    printf("call displayImage\n");
    this->displayImage(vtkImage);
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::getImageXY(int pos[], int xy[])
{
    printf("pos: %4d %4d\n",pos[0],pos[1]);
    fflush(stdout);
}

//---------------------------------------------------------------------------------------------------
// Define interaction style
//---------------------------------------------------------------------------------------------------
class MouseInteractorStyle4 : public vtkInteractorStyleImage
{
  public:
    static MouseInteractorStyle4* New();
    vtkTypeMacro(MouseInteractorStyle4, vtkInteractorStyleImage);

    virtual void OnLeftButtonDown() // This prevents the interactor from responding to L-button events
    {
        std::cout << "Pressed left mouse button." << std::endl;
        // for ImageViewer, not ImageTrimmer
        // int shifted = this->Interactor->GetShiftKey();
        // Forward events
        // if (shifted != 0) vtkInteractorStyleImage::OnLeftButtonDown();
    }
    virtual void OnRightButtonDown()    // This prevents the interactor from responding to R-button events
    {
        std::cout << "Pressed right mouse button." << std::endl;
    }
    virtual void OnMiddleButtonDown()   // This prevents the interactor from responding to M-button events
    {
        std::cout << "Pressed middle mouse button." << std::endl;
    }
};
vtkStandardNewMacro(MouseInteractorStyle4);


//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::displayImage(vtkImageData *image)
{
    int *r = renderWindow->GetSize();
    rect[0] = r[0];
    rect[1] = r[1];
    printf("rect: %d %d\n",rect[0],rect[1]);
    image->SetSpacing(1.0,1.0,1.0);
    int *dim = image->GetDimensions();
    double *spacing = image->GetSpacing(); 
    double *origin = image->GetOrigin();  

    float Cx = (dim[0] * spacing[0])/2. + origin[0];
    float Cy = (dim[1] * spacing[1])/2. + origin[1];

    printf("dim: %4d %4d\n",dim[0],dim[1]);

    camera->ParallelProjectionOn();
    camera->SetFocalPoint(Cx,Cy,0);    
    camera->SetPosition(Cx,Cy,1);   // insensitive to the z-value

    double scale, ratio[2];
    ratio[0] = (double)dim[0]/rect[0];
    ratio[1] = (double)dim[1]/rect[1];
    printf("ratio: %8.3f %8.3f\n",ratio[0],ratio[1]);
    if (ratio[0] > ratio[1]) {      // OK for dim[1] > dim[0], converse for dim[0] > dim[1]
        scale = dim[0]/2;
        scale *= (double)rect[1]/rect[0];
        // scaled to fit in the horizontal: dim[0] <-> rect[0]
        horizontalFit = true;
    } else {
        scale = dim[1]/2;
        // scaled to fit in the vertical: dim[1] <-> rect[1]
        horizontalFit = false;
    }
    fflush(stdout);

    camera->SetParallelScale(scale);
    // note that in each case the image centre coincides with the window centre

    //
    //    // to flip de image
    //    camera->SetViewUp (0, 1, 0);  
    //

    // set imageActor properties
    imageActor->SetInput(image);
    imageActor->InterpolateOff();

    renderer->AddActor(imageActor);
    renderer->SetActiveCamera(camera);

    // window interactor style for display images 
    vtkSmartPointer<MouseInteractorStyle4> style = vtkSmartPointer<MouseInteractorStyle4>::New();
    // set interactor style to the qvtkWidget Interactor
    qvtkWidget->GetInteractor()->SetInteractorStyle(style);

    qvtkWidget->GetInteractor()->Render();
    this->update();
}

//---------------------------------------------------------------------------------------------------
// Actor position is double pos[2] <--> pixel position int pixel[2].
// If horizontalFit:
//      pixel width <--> rect[0] => scale = rect[0]/width
// Else:
//      pixel height <--> rect[1] <--> scale = rect[1]/height
// where scale takes pixel -> pos
// In each case the image centre coincides with the window centre, i.e.:
//      pos[0] =rect[0]/2 + scale*(pixel[0] - width/2)
//      pos[1] = rect[1]/2 + scale*(pixel[1] - height/2)
//---------------------------------------------------------------------------------------------------
void ImageWidget:: convertPosToPixel(double *pos, int *pixel)
{
    double scale;
    if (horizontalFit) {
        scale = (double)rect[0]/imageWidth;
    } else {
        scale = (double)rect[1]/imageHeight;
    }
    pixel[0] = (pos[0] - rect[0]/2.)/scale + imageWidth/2. + 0.5;
    pixel[1] = (pos[1] - rect[1]/2.)/scale + imageHeight/2. + 0.5;
// Commented out for now to avoid obscuring an error
//    pixel[0] = qMax(pixel[0],0);
//    pixel[0] = qMin(pixel[0],imageWidth-1);
//    pixel[1] = qMax(pixel[1],0);
//    pixel[1] = qMin(pixel[1],imageHeight-1);
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget:: convertPixelToPos(int *pixel, double *pos)
{
    double scale;
    if (horizontalFit) {
        scale = (double)rect[0]/imageWidth;
    } else {
        scale = (double)rect[1]/imageHeight;
    }
    pos[0] = rect[0]/2 + scale*(pixel[0] - imageWidth/2);
    pos[1] = rect[1]/2 + scale*(pixel[1] - imageHeight/2);
// Commented out for now to avoid obscuring an error
//    pos[0] = qMax(pos[0],0);
//    pos[0] = qMin(pos[0],rect[0]);
//    pos[1] = qMax(pos[1],0);
//    pos[1] = qMin(pos[1],rect[1]);
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::tallyMark(QMouseEvent *event)
{
    QPoint p;
    int pos[3];
    vtkSmartPointer<vtkActor2D> actor;

    if (event->button()== Qt::LeftButton) {
        if (event->type() == event->MouseButtonPress) {
            LButtonDown = true;
            p = event->pos();
            pos[0] = p.x();
            pos[1] = rect[1] - p.y();
            pos[2] = 0;
            if (event->modifiers() & Qt::ShiftModifier){
                printf("Detected a shift-L-button mouse press event at: %d %d\n",p.x(),p.y());
                fflush(stdout);
                // identify the selected mark from the list (if there is one nearby)
                selected = selectMark(pos, &iselect);
                if (selected) printf("Selected mark: %d\n",iselect);
            } else if (event->modifiers() & Qt::ControlModifier){
                printf("Detected a ctrl-L-button mouse press event at: %d %d\n",p.x(),p.y());
                int iactor;
                if (selectMark(pos, &iactor)) {
                    printf("Selected to delete: %d\n",iactor);
                    deleteMark(iactor);
                }
            } else {
                printf("Detected a L-button mouse press event at: %d %d\n",p.x(),p.y());
            }
        }
    }
    if (event->button()== Qt::LeftButton) {
        if (event->type() == event->MouseButtonRelease) {
            LButtonDown = false;
            p = event->pos();
            pos[0] = p.x();
            pos[1] = rect[1] - p.y();
            pos[2] = 0;
            if (event->modifiers() & Qt::ShiftModifier){
                printf("Detected a shift-L-button mouse release event at: %d %d\n",p.x(),p.y());
                selected = false;
            } else if (event->modifiers() & Qt::ControlModifier){
                    printf("Detected a shift-L-button mouse release event at: %d %d\n",p.x(),p.y());
            } else {
                printf("Detected a L-button mouse release event at: %d %d\n",p.x(),p.y());
                addMark(pos);
            }
        }
    }
    if (LButtonDown) {
        if ((event->type() == event->MouseMove)) {
            if (event->modifiers() & Qt::ShiftModifier){
                p = event->pos();
                // if a mark has been selected, change its position in the list and on screen
                if (selected) {
                    actor = Actors.at(iselect);
                    pos[0] = p.x();
                    pos[1] = rect[1] - p.y();
                    pos[2] = 0;
                    moveMark(iselect,pos);
                }
            }
        }
    }
    fflush(stdout);
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
vtkSmartPointer<vtkActor2D> ImageWidget::makeMarkActor(double *pos)
{
    vtkSmartPointer<vtkActor2D> actor = vtkSmartPointer<vtkActor2D>::New();
    actor->SetMapper( mapper );
    actor->GetProperty()->SetColor( 1.0, 0.7, 0.0 ); // orange
    actor->SetPosition(pos[0], pos[1]);
    return actor;
}

//---------------------------------------------------------------------------------------------------
// from http://www.vtk.org/Wiki/VTK/Examples/Cxx/WishList/Images/MarkKeypoints
//---------------------------------------------------------------------------------------------------
void ImageWidget::addMark(int p[3])
{
    double pos[3];
    pos[0] = (double)p[0];
    pos[1] = (double)p[1];
    vtkSmartPointer<vtkActor2D> actor;
    actor = makeMarkActor(pos);

    int i;
    int nlist = Actors.size();
    if (nlist > 0) {
        bool ok;
        i = QInputDialog::getInt(this, tr("Mark sequence position"),
                                   tr("Insert after number:"), nlist-1, 0, nlist-1, 1, &ok);
        if (!ok) return;
    } else {
        i = -1;
    }
    vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
    if (i == nlist-1) {
        printf("Add to the end of the actors list\n");
        fflush(stdout);
        Actors.push_back(actor);
        NumActors.push_back(textActor);
        setNumber(Actors.size()-1);
    } else {
        printf("Insert into the actors list at: %d\n",i+1);
        fflush(stdout);
        std::vector<vtkSmartPointer<vtkActor2D>>::iterator it;
        it = Actors.begin();
        it = Actors.insert ( it+i+1 , actor );
        std::vector<vtkSmartPointer<vtkTextActor>>::iterator nit;
        nit = NumActors.begin();
        nit = NumActors.insert ( nit+i+1 , textActor );
        // Change numbers from i+1 to the end
        for (int j=i+1; j<NumActors.size(); j++) {
            setNumber(j);
        }
    }
    std::cout << "Actors.size: " << Actors.size() << std::endl;
    renderer->AddActor2D(actor);
    renderer->AddActor2D ( textActor );
    qvtkWidget->GetInteractor()->Render();

    emit(SaveStatus(true));
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::moveMark(int iactor, int p[3])
{
    vtkSmartPointer<vtkActor2D> actor = Actors.at(iactor);
    actor->SetPosition((double)p[0], (double)p[1]);
    printf("moved actor: %d to: %d %d\n",iactor,p[0],p[1]);
    fflush(stdout);
    // move number
    setNumber(iactor);
    qvtkWidget->GetInteractor()->Render();
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::deleteMark(int iactor)
{
    vtkSmartPointer<vtkActor2D> actor;
    actor = Actors.at(iactor);
    renderer->RemoveActor2D(actor);
    Actors.erase(Actors.begin()+iactor);
    vtkSmartPointer<vtkTextActor> textActor;
    textActor = NumActors.at(iactor);
    renderer->RemoveActor2D(textActor);
    NumActors.erase(NumActors.begin()+iactor);
    // Need to renumber from iactor to the end
    char numstr[8];
    int nlist = NumActors.size();
    for (int inum = iactor; inum < nlist; inum++) {
        textActor = NumActors.at(inum);
        _itoa(inum,numstr,10);
        textActor->SetInput (numstr);
    }
    std::cout << "Actors.size: " << Actors.size() << std::endl;
    qvtkWidget->GetInteractor()->Render();
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
bool ImageWidget::selectMark(int pos[3], int *iactor)
{
    double *apos;
    double delta, dist;
    vtkSmartPointer<vtkActor2D> actor;

    printf("selectMark: Actors.size: %d\n",Actors.size());
    fflush(stdout);
    for (int i=0; i<Actors.size(); i++) {
        printf("get actor: %d\n",i);
        fflush(stdout);
        actor = Actors.at(i);
        apos = actor->GetPosition();
        dist = 0;
        delta = pos[0] - apos[0];
        dist += delta*delta;
        delta = pos[1] - apos[1];
        dist += delta*delta;
        dist = sqrt(dist);
        if (dist <= threshold) {
            *iactor = i;
            return true;
        }
    }
    *iactor = -1;
    return false;
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::setThreshold(double t)
{
    threshold = t;
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::setNumber(int iactor)
{
    double *pos;
    char numstr[8];
    vtkSmartPointer<vtkActor2D> actor;
    actor = Actors.at(iactor);
    pos = actor->GetPosition();
//    std::cout << "setNumber: actor at: " << pos[0] << " " << pos[1] << std::endl;
    vtkSmartPointer<vtkTextActor> textActor;
    textActor = NumActors.at(iactor);
    textActor->GetTextProperty()->SetFontSize ( 14 );
    textActor->GetTextProperty()->SetBold(1);
    textActor->SetPosition ( pos[0]+4, pos[1]+4 );
    _itoa(iactor,numstr,10);
    textActor->SetInput (numstr);
    textActor->GetTextProperty()->SetColor ( 1.0,0.7,.0 );
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
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

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
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

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::writeLineFile(QString lineFileName)
{
    char s[100];
    QString str;
    int pixel[2];
    vtkSmartPointer<vtkActor2D> actor;

    QFile file(lineFileName);
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    out << getShortFileName() << endl;   // or "\n"
    str = QString::number(imageWidth) + " " + QString::number(imageHeight) + " " + QString::number(imageDepth);
    out << str << endl;
    out << QString::number(vect_actor_mat.size()) << endl;
    for (int n = 0; n < vect_actor_mat.size(); n++) {
        out << QString::number(vect_actor_mat[n].size()) << endl;
        for (int m = 0; m < vect_actor_mat[n].size(); m++) {
            actor = vect_actor_mat[n][m];
            double *pos = actor->GetPosition();
            convertPosToPixel(pos,pixel);
            sprintf(s,"%4d %2d %6.1f %6.1f %4d %4d",n,m,pos[0],pos[1],pixel[0],pixel[1]);
            out << s << endl;
        }
    }
    file.close();
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
void ImageWidget::readLineFile(QString lineFileName)
{
    QString line;
    int k;
    int width, height,depth;
//    int pixel[2];
    double pos[2];
    vtkSmartPointer<vtkActor2D> actor;
    QStringList myStringList;

    QFile file(lineFileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        // Error message
        std::cout << "File open failed on: " << lineFileName.toStdString() << endl;
        return;
    } else {
        std::cout << "File opened: " << lineFileName.toStdString() << endl;
    }
    QTextStream in(&file);
    line = in.readLine();
    std::cout << "image file name: " << line.toStdString() << endl;
    // Check that the file name is consistent
    if (line.compare(getShortFileName()) != 0) {
        std::cout << " image file name in file: " << line.toStdString() << " != " << getShortFileName().toStdString() << endl;
        return;
    }

    line = in.readLine();
    line = line.simplified();
    std::cout << line.toStdString() <<endl;
    myStringList = line.split(" ");
//    for (int index =0; index < myStringList.length(); index++)
//    {
//       std::cout << myStringList.at(index).toStdString() << endl;
//    }
    width = myStringList.at(0).toInt();
    height = myStringList.at(1).toInt();
    depth = myStringList.at(2).toInt();
    // Check that the image dimensions are consistent
    if (width != imageWidth) {
        std::cout << "Bad width: " << width << " != " << imageWidth << endl;
        return;
    }
    if (height != imageHeight) {
        std::cout << "Bad height: " << height << " != " << imageHeight << endl;
        return;
    }
    if (depth != imageDepth) {
        std::cout << "Bad depth: " << depth << " != " << imageDepth << endl;
        return;
    }
    clearActors();
    std::cout << "Did clearActors\n";
    line = in.readLine();
    line = line.simplified();
    int nlines = line.toInt();
    std::cout << "nlines: " << nlines << endl;
    if (nlines == 0) {
        // error message
        return;
    }
    int newz;
    int lastz = -1;
    for (int kline = 0; kline<nlines; kline++) {
        line = in.readLine();
        line = line.simplified();
        int nactors = line.toInt();
        if (nactors == 0) newz++;
        for (int kactor=0; kactor<nactors; kactor++) {
            QString line = in.readLine();
            line = line.simplified();
            myStringList = line.split(" ");
            // parse line -> newz,iactor,pos[0],pos[1],pixel[0],pixel[1]
            newz = myStringList.at(0).toInt();
            int iactor = myStringList.at(1).toInt();
            pos[0] = myStringList.at(2).toDouble();
            pos[1] = myStringList.at(3).toDouble();
            printf("newz, lastz, kactor, iactor: %d %d %d %d\n",newz,lastz,kactor,iactor);
            fflush(stdout);
            if (iactor != kactor) {
                // error message
                return;
            }
            if (kactor == 0 && newz != lastz+1) {   // need to create blank line entries
                Actors.clear();
                for (k=lastz+1; k<newz; k++) {
                    vect_actor_mat.push_back(Actors);
                    std::cout << "stored blank line\n";
                }
            }
            actor = makeMarkActor(pos);
            Actors.push_back(actor);
        }
        vect_actor_mat.push_back(Actors);
        std::cout << "stored line\n";
        Actors.clear();
        lastz = newz;
    }
    file.close();
    Actors.clear();
    showLines();
}
