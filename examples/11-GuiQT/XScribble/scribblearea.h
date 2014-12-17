#ifndef SCRIBBLEAREA_H
#define SCRIBBLEAREA_H

#include <QWidget>
#include "ui_scribblearea.h"
#include <QtGui/QPainter>
#include "shape.h"

class ScribbleArea : public QWidget
{
	Q_OBJECT

public:
	ScribbleArea(QWidget *parent = 0);
	~ScribbleArea();
	void initial();

	enum DrawFlag{
		FREE,
		LINE,
		CIRCLE,
		CLEAR
	};


protected:
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	void paintEvent(QPaintEvent *event);
	void resizeEvent(QResizeEvent *event);
	void resizeImage(QImage *image, const QSize &newSize);
public:
	const QImage& getimage() const;
	QImage& getimage();
	void clearall();

	public slots:
		void flagFree(){ m_drawFlag = FREE; }
		void flagLine(){ m_drawFlag = LINE; }
		void flagCircle(){ m_drawFlag = CIRCLE; }
		void flagclear(){ m_drawFlag = CLEAR; }
		void flagclearall();

private:
	Ui::ScribbleAreaClass ui;
	QPoint m_lastPoint;
	DrawFlag m_drawFlag;
	QImage m_image;
	int m_penWidth;
	QColor m_penColor;
	shape drawshape;//声明一个类的变量
};

#endif // SCRIBBLEAREA_H
