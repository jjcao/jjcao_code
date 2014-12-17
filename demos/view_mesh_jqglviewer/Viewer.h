#ifndef VIEWER_H__
#define VIEWER_H__

#include "MyItem.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <iostream>
#include <fstream>

#include <QGLViewer/qglviewer.h>
using namespace qglviewer;

class Viewer : public QGLViewer
{
	Q_OBJECT

signals:
	void vertsPicked(QString &str);
public slots:
	void invertFace();
	void invertNormal();
	void showWhole();
	void showScalar();
	void computeShortestDistance();

	void clearSelectedPoints();
	void invertSelectedPoints();
	void saveSelectedPoints();

public :
	Viewer(QWidget *parent);
	bool openMesh(const QString &fileName);
protected :
	virtual void init();
	virtual void draw();
	virtual void drawWithNames();
	//virtual void postSelection(const QPoint& point);
	virtual void endSelection(const QPoint& point);

	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseMoveEvent (QMouseEvent *);
	virtual void mouseReleaseEvent(QMouseEvent *e);

	virtual void keyPressEvent(QKeyEvent *e);
	virtual QString helpString() const;
private:
	void addIdToSelection(int id);
	void drawSelectionRectangle() const;
private:
	Polyhedron mesh_;
	std::list<Polyhedron::Vertex_iterator> pickedVertices_;
	bool wireframe_;
	bool frontFace_;
	bool flatShading_;
	double pointSize_;
	bool showScalar_;
	double scalarRange_[2];


	enum SelectionMode { NONE, ADD, REMOVE };
	SelectionMode selectionMode_;
	QRect rectangle_;
	std::vector<bool> vertPickedStatus_;
	bool beDraging_;
	Polyhedron::Vertex_iterator current_picked_vertex;
	Polyhedron::Point_3 move_p;
};

#endif