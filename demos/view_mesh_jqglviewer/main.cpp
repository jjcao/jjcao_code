#include "meshviewer.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	MeshViewer w;
	if ( argc > 1)
	{
		QString filename(argv[1]);
		w.open(filename);
	}
	else
	{
		w.open();
	}
	w.show();
	return a.exec();
}
