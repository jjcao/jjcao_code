#ifndef MESHVIEWER_H
#define MESHVIEWER_H

#include <QtGui/QMainWindow>
#include "ui_meshviewer.h"

class MeshViewer : public QMainWindow
{
	Q_OBJECT

public:
	MeshViewer(QWidget *parent = 0, Qt::WFlags flags = 0);
	~MeshViewer();
public slots:
	void open();
private slots:
	void logging(QString& str);
public:
	void open(QString& filename);
private:	
	void setupActions();
private:
	Ui::MeshViewerClass ui;
};

#endif // MESHVIEWER_H
