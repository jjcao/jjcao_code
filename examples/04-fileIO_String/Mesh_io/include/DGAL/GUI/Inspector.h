#ifndef INSPECTOR_H
#define INSPECTOR_H

#include "qt_vtk_global.h"

class QT_VTK_EXPORT Inspector : public QWidget
{
	Q_OBJECT

public:
	Inspector(QWidget *parent=0);
	~Inspector();
	void clear();
	void addWidget(QWidget* widget);
	void removeWidget(const char* name);
public slots:
	void setCurrentWidget(int index);

private:
	QStackedWidget *m_stackWidget;
	QScrollArea *m_scrollArea;
};

#endif // INSPECTOR_H
