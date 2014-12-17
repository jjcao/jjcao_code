#ifndef PROPERTYWIDGET_H
#define PROPERTYWIDGET_H

#include "qt_vtk_global.h"
#include <QWidget>
#include "ui_PropertyWidget.h"

class QT_VTK_EXPORT PropertyWidget : public QWidget
{
	Q_OBJECT

public:
	PropertyWidget(QWidget *parent = 0);
	~PropertyWidget();
signals:
	void representationChanged(int rep);
	void opacityChanged(double opacity);
public:
	void setRepresentation(int rep);
	void setOpacity(double opacity);
private:
	Ui::PropertyWidgetClass ui;
};

#endif // PROPERTYWIDGET_H
