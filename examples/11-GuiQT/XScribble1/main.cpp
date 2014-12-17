#include "xscribble.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	QPixmap pixmap(":/XScribble/Resources/splash.png");
	QSplashScreen *splash = new QSplashScreen(pixmap);
	splash->show();
	splash->showMessage("Loaded modules");
	qApp->processEvents();
	splash->showMessage("Established connections");
	qApp->processEvents();

	XScribble w;
	w.show();
	return a.exec();
}
