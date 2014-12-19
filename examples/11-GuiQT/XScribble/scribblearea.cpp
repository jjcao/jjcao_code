#include "scribblearea.h"
#include <QtGui/QMouseEvent>
#include <QtGui/QPainter>
#include <math.h>

ScribbleArea::ScribbleArea(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);
	initial();
}

ScribbleArea::~ScribbleArea()
{

}
void ScribbleArea::initial()
{
	m_penColor = Qt::green;
		m_penWidth = 3;
}
void ScribbleArea::mousePressEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton) {
		m_lastPoint = event->pos();
	}
}
void ScribbleArea::mouseMoveEvent(QMouseEvent *event)
{
	if ((event->buttons() & Qt::LeftButton) )
	{
		
		/*QPainter painter(&m_image);
		
		painter.setPen(QPen(m_penColor, m_penWidth, Qt::SolidLine, Qt::RoundCap,
			Qt::RoundJoin));*///为什么加了这些就不行了,已解决--这里占用了m_image,在这上面绑定了一个painter 却没用着，而他没有释放的话，别的函数上就绑定不了了

		switch (m_drawFlag)
		{
		case FREE:
			drawshape.PaintFree(event->pos(),m_image,m_penColor,m_penWidth,m_lastPoint);
			update();
			break;
		case LINE:
			break;
		case CIRCLE:
			drawshape.PaintCircle(event->pos(),m_image,m_penColor,m_penWidth,m_lastPoint);
			update();
			break;
		case CLEAR:
			drawshape.clear(event->pos(),m_image,m_penColor,m_penWidth,m_lastPoint);
			update();
			break;
		}
	}
}
void ScribbleArea::mouseReleaseEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton) 
	{
			switch (m_drawFlag)
			{
			case FREE:
			case CIRCLE:
			case CLEAR:
				break;
		    case LINE:
				drawshape.PaintLine(event->pos(),m_image,m_penColor,m_penWidth,m_lastPoint);
				update();
				break;
			}
	}
}
void ScribbleArea::paintEvent(QPaintEvent *event)
{
	QPainter painter(this);
	QRect dirtyRect = event->rect();
	painter.drawImage(dirtyRect, m_image, dirtyRect);
}
void ScribbleArea::resizeEvent(QResizeEvent *event){
	if (width() > m_image.width() || height() > m_image.height()) {
		int newWidth = qMax(width() + 128, m_image.width());
		int newHeight = qMax(height() + 128, m_image.height());
		resizeImage(&m_image, QSize(newWidth, newHeight));
		update();
	}
	QWidget::resizeEvent(event);
}
void ScribbleArea::resizeImage(QImage *image, const QSize &newSize){
	if (image->size() == newSize)
		return;
	QImage newImage(newSize, QImage::Format_RGB32);
	newImage.fill(qRgb(255, 255, 255));
	QPainter painter(&newImage);
	painter.drawImage(QPoint(0, 0), *image);
	*image = newImage;
}
void ScribbleArea::flagclearall()
{
	QPainter painter(&m_image);
	m_penColor = Qt::white;
	m_penWidth = 3;
	painter.setPen(QPen(m_penColor, m_penWidth, Qt::SolidLine, Qt::RoundCap,
		Qt::RoundJoin));
	painter.setBrush(QBrush(Qt::white));
	painter.drawRect(m_image.rect());
	update();
}
const QImage& ScribbleArea::getimage() const
{
	return this->m_image;
}
QImage& ScribbleArea::getimage()
{
	return this->m_image;
}
