#ifndef QTTAUCS_H
#define QTTAUCS_H

#include <QtGui/QMainWindow>
#include "ui_qttaucs.h"

class qtTaucs : public QMainWindow
{
	Q_OBJECT

public:
	qtTaucs(QWidget *parent = 0, Qt::WFlags flags = 0);
	~qtTaucs();

private:
	Ui::qtTaucsClass ui;

private slots:
	void on_actionTest_triggered();
};

#endif // QTTAUCS_H
