#ifndef SHAPE_H
#define SHAPE_H
#include <QtGui/QPainter>

class shape
{
public:
	shape(){
	}
	~shape(){
	}
	void PaintFree(const QPoint& pt,QImage& m_image,QColor m_penColor,int m_penWidth,QPoint& m_lastPoint);//参数中QImage 后的&不能少
	void PaintCircle(const QPoint& pt,QImage& m_image,QColor m_penColor,int m_penWidth,QPoint& m_lastPoint);
	void PaintLine(const QPoint& pt,QImage& m_image,QColor m_penColor,int m_penWidth,QPoint& m_lastPoint);
	void clear(const QPoint& pt,QImage& m_image,QColor m_penColor,int m_penWidth,QPoint& m_lastPoint);
private:
	float r;
	float mid_r;
	QPoint m_midpoint;

};
#endif //SHAPE_H
