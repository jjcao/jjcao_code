#ifndef SCRIBBLEAREA_H
#define SCRIBBLEAREA_H

#include <QWidget>
#include <QColor>
#include <QImage>
#include <QPoint>
#include <QtGui>
#include "ui_scribblearea.h"

class ScribbleArea : public QWidget
{
	Q_OBJECT
public:
	enum DrawFlag {
		FREE,
        LINE,
        CIRCLE
    };

public:
	ScribbleArea(QWidget *parent = 0);
	~ScribbleArea();
	bool isModified() const { return m_modified; }
	bool openImage(const QString &fileName);
	bool saveImage(const QString &fileName, const char *fileFormat);
public slots:
	void flagFree(){	m_drawFlag=FREE; }
	void flagLine(){	m_drawFlag=LINE; }
	void flagCircle(){ m_drawFlag=CIRCLE; }
	void setPenColor();
protected:
    void mousePressEvent(QMouseEvent *event);// overriding 림맨
    void mouseMoveEvent(QMouseEvent *event);// overriding 림맨
    void mouseReleaseEvent(QMouseEvent *event);// overriding 림맨
    void paintEvent(QPaintEvent *event);// overriding 림맨
	void resizeEvent(QResizeEvent *event);// overriding 림맨
private:
	void drawFree(const QPoint& pt);
	void resizeImage(QImage *image, const QSize &newSize);
private:
	QPoint m_lastPoint;
	int m_drawFlag;
    int m_penWidth;
	QPen m_pen;
	QPainter m_painter;
    QImage m_image;
	bool m_modified;
private:
	Ui::ScribbleArea ui;
};

#endif // SCRIBBLEAREA_H
