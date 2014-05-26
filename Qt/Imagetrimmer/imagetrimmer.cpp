#include <QtGui>

#include "imagetrimmer.h"
#include "ui_mainwindow.h"

ImageTrimmer::ImageTrimmer(QWidget *parent) : QMainWindow(parent),
    ui(new Ui::ImageTrimmer)
{
    ui->setupUi(this);
    createActions();
	createMenus();
	useITK = true;
    imageWidget1 = NULL;
    connect(imageWidget1,SIGNAL(SaveStatus(bool)),this,SLOT(setLineSaveStatus(bool)));
}

ImageTrimmer::~ImageTrimmer()
{
    delete ui;
    imageWidget1 = NULL;
}

void ImageTrimmer::open1()
{
    opener(1);
}

void ImageTrimmer::opener(int imageNum)
{
    if (!imageWidget1) {
        printf("not imageWidget1\n");
        fflush(stdout);
        this->imageWidget1 = new ImageWidget(this, ui->mdiArea1);
        printf("created imageWidget1\n");
        fflush(stdout);
        if (useITK)
            imageWidget1->openWithITK();
        else
            imageWidget1->openWithVTK();
    } else {
        printf("is imageWidget1\n");
        fflush(stdout);
        imageWidget1->clearAll();
        printf("did clearAll\n");
        fflush(stdout);
        if (useITK)
            imageWidget1->openWithITK();
        else
            imageWidget1->openWithVTK();
    }
    printf("file name: %s\n",imageWidget1->getShortFileName().toAscii().constData());
    ui->labelImage1->setText(imageWidget1->getShortFileName());
    ui->lineEditFrame1->setText("0");
    if (imageWidget1) {
        iframe = 0;
        ui->groupBoxFrame->setEnabled(true);
        ui->spinBoxFrame->setValue(iframe);
        int maxframe = 0;
        maxframe = MAX(maxframe,imageWidget1->getDepth()-1);
        ui->spinBoxFrame->setMaximum(maxframe);
        ui->pushButtonLoadLineFile->setEnabled(true);
    } else {
        ui->groupBoxFrame->setEnabled(false);
    }
    on_lineEditThreshold_textChanged();
}

void ImageTrimmer::prevFrame()
{
    int ifrm = ui->spinBoxFrame->value();
    ifrm--;
    if (ifrm < 0) ifrm = 0;
    ui->spinBoxFrame->setValue(ifrm);
    updateImages();
}

void ImageTrimmer::nextFrame()
{
    int ifrm = ui->spinBoxFrame->value();
    ifrm++;
    int maxdepth = imageWidget1->getDepth();
    if (ifrm >= maxdepth) ifrm = maxdepth-1;
    ui->spinBoxFrame->setValue(ifrm);
    updateImages();
}

void ImageTrimmer::updateImages()
{
    int current_iframe;

    current_iframe = iframe;
    printf("updateImages\n");
    if (!imageWidget1 ) {
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
}

void ImageTrimmer::setLineSaveStatus(bool saveable)
{
    std::cout << "Received LineSaveStatus signal\n";
    ui->pushButtonSaveLine->setEnabled(saveable);
//        ui->pushButtonWriteLineFile->setEnabled(true);
}

void ImageTrimmer::saveLine()
{
    if (imageWidget1) {
        imageWidget1->saveLine();
    }
}

void ImageTrimmer::saveLineMatrix()
{
    QString fileName;
    fileName = QFileDialog::getSaveFileName(this,tr("Open Save Cutline File"), "", tr("Cutline Files (*.dat)"));
    imageWidget1->writeLineFile(fileName);
}

void ImageTrimmer::loadLineMatrix()
{
    QString fileName;
    fileName = QFileDialog::getOpenFileName(this,tr("Open Load Cutline File"), "", tr("Cutline Files (*.dat)"));
    imageWidget1->readLineFile(fileName);
    iframe = -1;
    updateImages();
}

void ImageTrimmer::on_lineEditThreshold_textChanged()
{
    double threshold = ui->lineEditThreshold->text().toDouble();
    imageWidget1->setThreshold(threshold);
}


void ImageTrimmer::about()
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

void ImageTrimmer::createActions()
{
    openAct1 = new QAction(tr("Open Image 1..."), this);
    openAct1->setShortcut(tr("Ctrl+1"));
    connect(openAct1, SIGNAL(triggered()), this, SLOT(open1()));
	exitAct = new QAction(tr("E&xit"), this);
	exitAct->setShortcut(tr("Ctrl+Q"));
	exitAct->setStatusTip(tr("Exit the application"));
	aboutAct = new QAction(tr("&About"), this);
	connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));

	aboutQtAct = new QAction(tr("About &Qt"), this);
	connect(aboutQtAct, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void ImageTrimmer::createMenus()
{
	fileMenu = new QMenu(tr("&File"), this);
    fileMenu->addAction(openAct1);
	fileMenu->addSeparator();
	fileMenu->addAction(exitAct);
	helpMenu = new QMenu(tr("&Help"), this);
	helpMenu->addAction(aboutAct);
	helpMenu->addAction(aboutQtAct);

	menuBar()->addMenu(fileMenu);
	menuBar()->addMenu(helpMenu);
}

void ImageTrimmer::createStatusBar()
{
	statusLabel = new QLabel(" Basic Image Viewer     ");
	statusLabel->setAlignment(Qt::AlignHCenter);
    statusLabel->setMaximumSize(statusLabel->sizeHint());

	statusBar()->addWidget(statusLabel);
}
