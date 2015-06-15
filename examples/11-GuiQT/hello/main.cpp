#include <QtWidgets/QApplication>
#include <QtWidgets/QLabel>
#include <QtWidgets/qmessagebox.h>
#include <qtimer.h>

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    QLabel l( "Hello World!" );
	l.show();

	{
		QMessageBox tbox;
		tbox.setText(QString("%1 ms - cost = %2").arg(11).arg(22));
		//tbox.setStandardButtons(QMessageBox::Ok);
		tbox.exec();
		//tbox.show(); //QTimer::singleShot(60, &tbox, SLOT(hide()));
	}

    return app.exec();
}
