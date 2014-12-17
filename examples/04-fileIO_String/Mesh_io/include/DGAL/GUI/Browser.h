#ifndef BROWSER_H
#define BROWSER_H

#include "qt_vtk_global.h"

class QT_VTK_EXPORT Browser : public QWidget
{
	Q_OBJECT

public:
	Browser(QWidget *parent=0);
	~Browser();
	QToolButton *addItem(const char* name, bool bVisible = true);
	void removeItem(const char* name);
	void clear();
	void setColumnWidth(int colId, int width);
signals:
	void currentItemChanged(int i);
private slots:
	void cellCliked(int i, int j);	
private:
	QTableWidget* m_tableWidget;
	QIcon m_icon;
	QVector<QToolButton *> m_buttons;
	QSize m_iconSize;
};

#endif // BROWSER_H
