#ifndef XSCRIBBLE_H
#define XSCRIBBLE_H

#include <QtGui/QMainWindow>
#include "ui_xscribble.h"

class XScribble : public QMainWindow
{
	Q_OBJECT

public:
	XScribble(QWidget *parent = 0, Qt::WFlags flags = 0);
	~XScribble();
public slots:
   void saveFile();
   void openFile();
private:
	Ui::XScribbleClass ui;
};

#endif // XSCRIBBLE_H
