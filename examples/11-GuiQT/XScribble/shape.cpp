#include "shape.h"
#include "scribblearea.h"
#include <QtGui/QPainter>
#include <math.h>

void shape::PaintFree(const QPoint& pt,QImage& m_image,QColor m_penColor,int m_penWidth,QPoint& m_lastPoint)
{
    QPainter painter(&m_image);
		painter.setPen(QPen(m_penColor, m_penWidth, Qt::SolidLine, Qt::RoundCap,
			Qt::RoundJoin));
	painter.drawLine(m_lastPoint, pt);
			m_lastPoint = pt;
}
void shape:: PaintCircle(const QPoint& pt,QImage& m_image,QColor m_penColor,int m_penWidth,QPoint& m_lastPoint)
{
	QPainter painter(&m_image);
		painter.setPen(QPen(m_penColor, m_penWidth, Qt::SolidLine, Qt::RoundCap,
			Qt::RoundJoin));

		float r;
			//repaint the last circle
			if (m_midpoint != pt)
		 {
		 m_penColor = QColor(255,255,255);
		 painter.setPen(QPen(m_penColor, m_penWidth, Qt::SolidLine, Qt::RoundCap,
		 Qt::RoundJoin));
		 painter.drawEllipse(m_lastPoint , (int) mid_r, (int) mid_r);
		 m_penColor = Qt::green;
		 painter.setPen(QPen(m_penColor, m_penWidth, Qt::SolidLine, Qt::RoundCap,
		 Qt::RoundJoin));
		 }
			r = sqrt(float((m_lastPoint.x() - pt.x()) *(m_lastPoint.x() - pt.x()) + (m_lastPoint.y() - pt.y()) *(m_lastPoint.y() - pt.y())));
			painter.drawEllipse(m_lastPoint , (int) r, (int) r);
			m_midpoint = pt;
			mid_r = r;
}
void shape::PaintLine(const QPoint& pt,QImage& m_image,QColor m_penColor,int m_penWidth,QPoint& m_lastPoint)
{
	QPainter painter(&m_image);
		        m_penColor = Qt::green;
		        m_penWidth = 3;
		        painter.setPen(QPen(m_penColor, m_penWidth, Qt::SolidLine, Qt::RoundCap,
		        	Qt::RoundJoin));
				painter.drawLine(m_lastPoint, pt);
}
void shape::clear(const QPoint& pt,QImage& m_image,QColor m_penColor,int m_penWidth,QPoint& m_lastPoint)
{
	QPainter painter(&m_image);
	m_penColor = QColor(255,255,255);
			painter.setPen(QPen(m_penColor, m_penWidth, Qt::SolidLine, Qt::RoundCap,
				Qt::RoundJoin));
			painter.drawLine(m_lastPoint, pt);
			m_lastPoint = pt;
}