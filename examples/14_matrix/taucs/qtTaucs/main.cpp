#include <QtGui/QApplication>
#include "qttaucs.h"

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	qtTaucs w;
	w.show();
	a.connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));
	return a.exec();
}
