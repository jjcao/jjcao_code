#include "xscribble.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	XScribble w;
	w.show();
	return a.exec();
}
