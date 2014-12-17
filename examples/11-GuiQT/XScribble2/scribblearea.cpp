#include "scribblearea.h"

ScribbleArea::ScribbleArea(QWidget *parent)
	: QWidget(parent), m_drawFlag(FREE),m_modified(false)
{
	ui.setupUi(this);
	
	setAttribute(Qt::WA_StaticContents);//Indicates that the widget contents are north-west aligned and static. On resize, such a widget will receive paint events only for parts of itself that are newly visible.

	m_pen = QPen(Qt::blue, 1, Qt::SolidLine, Qt::RoundCap,Qt::RoundJoin);
}

ScribbleArea::~ScribbleArea()
{
}

bool ScribbleArea::openImage(const QString &fileName)
{
    QImage loadedImage;
    if (!loadedImage.load(fileName))
        return false;

    QSize newSize = loadedImage.size().expandedTo(size());// union of loadedImage.size() & ScribbleArea's size.
    resizeImage(&loadedImage, newSize);
    m_image = loadedImage;
    m_modified = false;
    update();
    return true;
}
bool ScribbleArea::saveImage(const QString &fileName, const char *fileFormat)
{
    if (m_image.save(fileName, fileFormat)) {
        m_modified = false;
        return true;
    } else {
        return false;
    }
}
void ScribbleArea::setPenColor()
{
	QColor oldColor = m_pen.color();
	QColor color = QColorDialog::getColor(oldColor, this);
	m_pen.setColor(color);
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
		switch (m_drawFlag)
		{
		case FREE:
			drawFree(event->pos());
			m_modified = true;
			break;
		}
	}
}
void ScribbleArea::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) 
	{
    }
}
void ScribbleArea::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    QRect dirtyRect = event->rect();
    painter.drawImage(dirtyRect, m_image, dirtyRect);
}
void ScribbleArea::resizeEvent(QResizeEvent *event)
{
    if (width() > m_image.width() || height() > m_image.height()) {
        int newWidth = qMax(width() + 128, m_image.width());
        int newHeight = qMax(height() + 128, m_image.height());
        resizeImage(&m_image, QSize(newWidth, newHeight));
        update();
    }
    QWidget::resizeEvent(event);
}
void ScribbleArea::resizeImage(QImage *image, const QSize &newSize)
{
    if (image->size() == newSize)
        return;

    QImage newImage(newSize, QImage::Format_RGB32);
    newImage.fill(qRgb(255, 255, 255));
    QPainter painter(&newImage);
    painter.drawImage(QPoint(0, 0), *image);
    *image = newImage;
}
void ScribbleArea::drawFree(const QPoint& pt)
{
    QPainter painter(&m_image);
    painter.setPen(m_pen);
    painter.drawLine(m_lastPoint, pt);

	// update what you just drew
	int rad = (m_pen.width() / 2) + 2;
    update(QRect(m_lastPoint, pt).normalized().adjusted(-rad, -rad, +rad, +rad));

    m_lastPoint = pt;
}

