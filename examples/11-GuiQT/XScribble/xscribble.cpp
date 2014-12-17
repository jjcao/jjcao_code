#include "xscribble.h"
#include <QtGui/QFileDialog>

XScribble::XScribble(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	ui.setupUi(this);
	connect(ui.actionSave,SIGNAL(triggered()),this,SLOT(saveFile()));
	connect(ui.actionOpen,SIGNAL(triggered()),this,SLOT(openFile()));
	connect(ui.actionDrawLine,SIGNAL(triggered()),ui.centralWidget,SLOT(flagLine()));
   connect(ui.actionDrawCircle,SIGNAL(triggered()),ui.centralWidget,SLOT(flagCircle()));
   connect(ui.actionFreeHand,SIGNAL(triggered()),ui.centralWidget,SLOT(flagFree()));
   connect(ui.actionClear,SIGNAL(triggered()),ui.centralWidget,SLOT(flagclear()));
   connect(ui.actionClearall,SIGNAL(triggered()),ui.centralWidget,SLOT(flagclearall()));
}


XScribble::~XScribble()
{

}
void XScribble::saveFile()
{
	QString fileName = QFileDialog::getSaveFileName(this,tr("Save Image"),"",tr("Image (*.png)"));
	if (fileName.isEmpty())
		return;
	this->ui.centralWidget->getimage().save(fileName);	
}
void XScribble::openFile()
{
	QString fileName = QFileDialog::getOpenFileName(this,tr("Open Image"),"",tr("Image (*.png)"));
	if (fileName.isEmpty())
		return;
	this->ui.centralWidget->getimage().load(fileName);	
}