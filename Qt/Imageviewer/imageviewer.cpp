#include <QtGui>

#include "imageviewer.h"
#include "ui_mainwindow.h"

ImageViewer::ImageViewer(QWidget *parent) : QMainWindow(parent),
    ui(new Ui::ImageViewer)
{
    ui->setupUi(this);

//    setAttribute(Qt::WA_DeleteOnClose);

    createActions();
	createMenus();
//	createStatusBar();
	useITK = true;
    imageWidget1 = NULL;
    imageWidget2 = NULL;
}

ImageViewer::~ImageViewer()
{
    delete ui;
    imageWidget1 = NULL;
    imageWidget2 = NULL;

}

void ImageViewer::open1()
{
    opener(1);
}

void ImageViewer::open2()
{
    opener(2);
}

void ImageViewer::opener(int imageNum)
{
    if (imageNum ==1) {
        if (!imageWidget1) {
            printf("not imageWidget1\n");
    //        this->imageWidget1 = new ImageWidget();
            this->imageWidget1 = new ImageWidget(ui->mdiArea1);

            if (useITK)
                imageWidget1->openWithITK();
            else
                imageWidget1->openWithVTK();
            printf("file name: %s\n",imageWidget1->getShortFileName().toAscii().constData());
            ui->labelImage1->setText(imageWidget1->getShortFileName());
            ui->lineEditFrame1->setText("0");
    //      this->setCentralWidget(imageWidget1);
    //		this->setWindowTitle(tr("Image Viewer"));
    //		this->resize(640, 480);
        } else {
            printf("is imageWidget1\n");
            /*
            ImageViewer *viewer = new ImageViewer();
    //        viewer->imageWidget1 = new ImageWidget();
            viewer->imageWidget1 = new ImageWidget(ui->mdiArea);

            if (useITK)
                viewer->imageWidget1->openWithITK();
            else
                viewer->imageWidget1->openWithVTK();
            */
            imageWidget1->clearImages();
            printf("did clearImages\n");
            if (useITK)
                imageWidget1->openWithITK();
            else
                imageWidget1->openWithVTK();

    //		viewer->setCentralWidget(viewer->imageWidget1);
    //		viewer->setWindowTitle(tr("Image Viewer"));
    //		viewer->resize(640, 480);
    //		viewer->show();
        }
    } else if (imageNum == 2) {
        if (!imageWidget2) {
            printf("not imageWidget2\n");
            this->imageWidget2 = new ImageWidget(ui->mdiArea2);

            if (useITK)
                imageWidget2->openWithITK();
            else
                imageWidget2->openWithVTK();
            printf("file name: %s\n",imageWidget2->getShortFileName().toAscii().constData());
            ui->labelImage2->setText(imageWidget2->getShortFileName());
            ui->lineEditFrame2->setText("0");
        } else {
            printf("is imageWidget2\n");
            imageWidget2->clearImages();
            printf("did clearImages\n");
            if (useITK)
                imageWidget2->openWithITK();
            else
                imageWidget2->openWithVTK();
        }
    }
    if (imageWidget1 || imageWidget2) {
        if (imageWidget1 && imageWidget2) {
            if (imageWidget1->getDepth() != imageWidget2->getDepth()) {
                printf("Different image depths\n");
            }
        }
        iframe = 0;
        ui->groupBoxFrame->setEnabled(true);
        ui->spinBoxFrame->setValue(iframe);
        int maxframe = 0;
        if (imageWidget1)
            maxframe = MAX(maxframe,imageWidget1->getDepth()-1);
        if (imageWidget2)
            maxframe = MAX(maxframe,imageWidget2->getDepth()-1);
        ui->spinBoxFrame->setMaximum(maxframe);
    } else {
        ui->groupBoxFrame->setEnabled(false);
    }

}

void ImageViewer::prevFrame()
{
    int ifrm = ui->spinBoxFrame->value();
    ifrm--;
    if (ifrm < 0) ifrm = 0;
    ui->spinBoxFrame->setValue(ifrm);
    updateImages();
}

void ImageViewer::nextFrame()
{
    int ifrm = ui->spinBoxFrame->value();
    ifrm++;
    int max1 = 0;
    int max2 = 0;
    if (imageWidget1) max1 = imageWidget1->getDepth();
    if (imageWidget2) max2 = imageWidget2->getDepth();
    int maxdepth = MAX(max1,max2);
    if (ifrm >= maxdepth) ifrm = maxdepth-1;
    ui->spinBoxFrame->setValue(ifrm);
    updateImages();
}

void ImageViewer::updateImages()
{
    int current_iframe;

    current_iframe = iframe;
    printf("updateImages\n");
    if (!imageWidget1 && !imageWidget2 ) {
        printf("no imageWidgets\n");
        return;
    }
    iframe = ui->spinBoxFrame->value();
    if (iframe == current_iframe) {
        printf("no change to iframe: %d\n",iframe);
        return;
    }
    printf("iframe: %d\n",iframe);
    QString framestr = QString::number(iframe);
    if (imageWidget1) {
        if (imageWidget1->frameOK(iframe)) {
            imageWidget1->showFrame(iframe);
            ui->lineEditFrame1->setText(framestr);
        }
    }
    if (imageWidget2) {
        if (imageWidget2->frameOK(iframe)) {
            imageWidget2->showFrame(iframe);
            ui->lineEditFrame2->setText(framestr);
        }
    }
}

void ImageViewer::imageSelecter() {
    if (ui->radioButtonImage1->isChecked()) {
        ui->pushButtonSubtractImage->setText("Subtract Image 2");
        ui->pushButtonAddImage->setText("Add Image 2");
        ui->pushButtonSaveImage->setText("Save Image 1");
     } else {
        ui->pushButtonSubtractImage->setText("Subtract Image 1");
        ui->pushButtonAddImage->setText("Add Image 1");
        ui->pushButtonSaveImage->setText("Save Image 2");
    }
}

bool ImageViewer::imageEquivalence() {
    int n1, w1, h1, d1, n2, w2, h2, d2;

    if (!imageWidget1 || !imageWidget2) return false;
    imageWidget1->getInfo(&n1,&w1,&h1,&d1);
    imageWidget2->getInfo(&n2,&w2,&h2,&d2);
    if (n1 != n2 || w1 != w2 || h1 != h2 || d1 != d2) return false;
    return true;
}

void ImageViewer::subtracter(){
    int val = ui->lineEditNumber->text().toInt();
    imageWidget1->subtract(val);
    imageWidget1->showFrame(iframe);
}

void ImageViewer::subtractImage()
{
    unsigned char *p;

    if (!imageEquivalence()) return;
    if (ui->radioButtonImage1->isChecked()) {
        p = imageWidget2->getBuffer();
        imageWidget1->subtractImage(p);
        imageWidget1->showFrame(iframe);
    } else {
        p = imageWidget1->getBuffer();
        imageWidget2->subtractImage(p);
        imageWidget2->showFrame(iframe);
    }
}

void ImageViewer::on_pushButtonAddImage_clicked()
{
    unsigned char *p;

    if (!imageEquivalence()) return;
    if (ui->radioButtonImage1->isChecked()) {
        p = imageWidget2->getBuffer();
        imageWidget1->addImage(p);
        imageWidget1->showFrame(iframe);
    } else {
        p = imageWidget1->getBuffer();
        imageWidget2->addImage(p);
        imageWidget2->showFrame(iframe);
    }

}

void ImageViewer::saveImageFile() {
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"), QDir::currentPath());
    if (ui->radioButtonImage1->isChecked()) {
        imageWidget1->saveImage(fileName);
        ui->labelImage1->setText(imageWidget1->getShortFileName());
    } else {
        imageWidget2->saveImage(fileName);
        ui->labelImage2->setText(imageWidget2->getShortFileName());
    }
}

void ImageViewer::saveAs()
{

}

//void ImageViewer::medianFilter()
//{
//	this->imageWidget1->medianFilter();
//}

//void ImageViewer::gradientAnisotropicDiffusionFilter()
//{
//	this->imageWidget1->gradientAnisotropicDiffusionFilter();
//}

void ImageViewer::about()
{
	QMessageBox::about(this, tr("About Image Viewer"),
		tr("<p>The <b>Image Viewer</b> example shows how to combine QLabel "
		"and QScrollArea to display an image. QLabel is typically used "
		"for displaying a text, but it can also display an image. "
		"QScrollArea provides a scrolling view around another widget. "
		"If the child widget exceeds the size of the frame, QScrollArea "
		"automatically provides scroll bars. </p><p>The example "
		"demonstrates how QLabel's ability to scale its contents "
		"(QLabel::scaledContents), and QScrollArea's ability to "
		"automatically resize its contents "
		"(QScrollArea::widgetResizable), can be used to implement "
		"zooming and scaling features. </p><p>In addition the example "
		"shows how to use QPainter to print an image.</p>"));
}

void ImageViewer::createActions()
{
    openAct1 = new QAction(tr("Open Image 1..."), this);
    openAct1->setShortcut(tr("Ctrl+1"));
    connect(openAct1, SIGNAL(triggered()), this, SLOT(open1()));

    openAct2 = new QAction(tr("Open Image 2..."), this);
    openAct2->setShortcut(tr("Ctrl+2"));
    connect(openAct2, SIGNAL(triggered()), this, SLOT(open2()));

    saveAsAct = new QAction(tr("Save As..."), this);
    saveAsAct->setShortcut(tr("Ctrl+Shift+S"));
    connect(saveAsAct, SIGNAL(triggered()), this, SLOT(saveAs()));

	exitAct = new QAction(tr("E&xit"), this);
	exitAct->setShortcut(tr("Ctrl+Q"));
	exitAct->setStatusTip(tr("Exit the application"));
	connect(exitAct, SIGNAL(triggered()), qApp, SLOT(closeAllWindows()));

//    medianFilterAct = new QAction(tr("Median Filter"), this);
//    medianFilterAct->setStatusTip(tr("Apply a median filter to image"));
//    connect(medianFilterAct, SIGNAL(triggered()), this, SLOT(medianFilter()));

//    GADFilterAct = new QAction(tr("Gradient Anisotropic Filter"), this);
//    connect(GADFilterAct, SIGNAL(triggered()), this, SLOT(gradientAnisotropicDiffusionFilter()));

	aboutAct = new QAction(tr("&About"), this);
	connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));

	aboutQtAct = new QAction(tr("About &Qt"), this);
	connect(aboutQtAct, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void ImageViewer::createMenus()
{
	fileMenu = new QMenu(tr("&File"), this);
    fileMenu->addAction(openAct1);
    fileMenu->addAction(openAct2);
    fileMenu->addAction(saveAsAct);
	fileMenu->addSeparator();
	fileMenu->addAction(exitAct);

//    filterMenu = new QMenu(tr("&Filter"), this);
//    filterMenu->addAction(medianFilterAct);
//    filterMenu->addAction(GADFilterAct);

	helpMenu = new QMenu(tr("&Help"), this);
	helpMenu->addAction(aboutAct);
	helpMenu->addAction(aboutQtAct);

	menuBar()->addMenu(fileMenu);
//	menuBar()->addMenu(filterMenu);
	menuBar()->addMenu(helpMenu);
}

void ImageViewer::createStatusBar()
{
	statusLabel = new QLabel(" Basic Image Viewer     ");
	statusLabel->setAlignment(Qt::AlignHCenter);
    statusLabel->setMaximumSize(statusLabel->sizeHint());

	statusBar()->addWidget(statusLabel);
}
