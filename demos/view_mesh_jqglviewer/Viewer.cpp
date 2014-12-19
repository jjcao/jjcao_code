#include "Viewer.h"
#include <QMenu>
#include <QKeyEvent>
#include <algorithm>
#include "graphAlgorithms.h"
#include "glut/glut.h"
#include <sstream>

typedef Polyhedron::Point_3 Point3;

void compute_bounding_box(Polyhedron &mesh, qglviewer::Vec &minv, qglviewer::Vec &maxv)
{
	if(mesh.empty())
		return;

	Polyhedron::Vertex_iterator pv = mesh.vertices_begin();
	minv.x = maxv.x = pv->point().x();
	minv.y = maxv.y = pv->point().y();
	minv.z = maxv.z = pv->point().z();
	for(;	pv !=  mesh.vertices_end(); pv++)
	{
		const Kernel::Point_3& p = pv->point();

		minv.x =  min(minv.x,p.x());
		minv.y =  min(minv.y,p.y());
		minv.z =  min(minv.z,p.z());

		maxv.x =  max(maxv.x,p.x());
		maxv.y =  max(maxv.y,p.y());
		maxv.z =  max(maxv.z,p.z());
	}
}
void setupIndex(Polyhedron &mesh, double* scalarRange, std::vector<bool> &vertPickedStatus)
{
	if(mesh.empty())
		return;

	int i(0);
	Polyhedron::Vertex_iterator pv = mesh.vertices_begin();
	scalarRange[0] = scalarRange[1] = pv->point().x();
	for(; pv !=  mesh.vertices_end(); ++pv, ++i)
	{
		pv->index_ = i;
		pv->scalar_ = pv->point().x();
		if ( pv->scalar_ < scalarRange[0])
			scalarRange[0] = pv->scalar_;
		if ( pv->scalar_ > scalarRange[1])
			scalarRange[1] = pv->scalar_;

		vertPickedStatus.push_back(false);
	}
}
double minEdgeLen(Polyhedron &mesh, double threshold)
{
	if(mesh.empty())
		return 0;

	double ml(threshold), tmp;
	for(Polyhedron::Halfedge_iterator pe = mesh.edges_begin();	 pe !=  mesh.edges_end(); ++pe)
	{
		Kernel::Vector_3 vec = pe->vertex()->point()-	pe->opposite()->vertex()->point();
		tmp = std::sqrt(vec*vec);
		if ( tmp < ml)
			ml = tmp;
	}
	return ml;
}
number_type average_edge_length_around(Polyhedron::Vertex_iterator pv)
{
	number_type sum = 0.0;
	Polyhedron::Halfedge_around_vertex_circulator pHalfEdge = pv->vertex_begin();
	Polyhedron::Halfedge_around_vertex_circulator end = pHalfEdge;
	Kernel::Vector_3 vec(0.0,0.0,0.0);
	int degree = 0;
	CGAL_For_all(pHalfEdge,end)
	{
		Kernel::Vector_3 vec = pHalfEdge->vertex()->point()-	pHalfEdge->opposite()->vertex()->point();
		sum += std::sqrt(vec*vec);
		++degree;
	}
	return sum / (number_type) degree;
}
bool Viewer::openMesh(const QString &fileName)
{
	mesh_.clear();
	pickedVertices_.clear();
	vertPickedStatus_.clear();

	beDraging_=false;
	current_picked_vertex = 0;

	/////////////////////
    std::ifstream stream(fileName.toUtf8() );
    stream >> mesh_;
    if(!stream || !mesh_.is_valid() || mesh_.empty())
    {
        std::cerr << "Error: cannot read OFF file " << fileName.toStdString() << std::endl;
        return EXIT_FAILURE;
    }
	std::for_each(mesh_.facets_begin(), mesh_.facets_end(), Face_normal());
	std::for_each(mesh_.vertices_begin(), mesh_.vertices_end(), Vertex_normal());
	setupIndex(mesh_, scalarRange_, vertPickedStatus_);

	// setting camera so that the loaded mesh can be viewed entirely
	qglviewer::Vec minv(0,0,0), maxv(1,1,1);
	compute_bounding_box(mesh_, minv, maxv);	
	camera()->setSceneBoundingBox(minv, maxv);
	pointSize_ = minEdgeLen(mesh_, maxv.x-minv.x) * 0.3;

	updateGL(); // update draw
	//update(); // update for qt
	return EXIT_SUCCESS;
}
void Viewer::computeShortestDistance()
{
	if (pickedVertices_.empty())
		return;

	Polyhedron::Vertex_iterator vi = *pickedVertices_.begin();
	dijkstra(mesh_, vi->index_, scalarRange_);
	showScalar_ = true;
	updateGL(); // update draw
}
Viewer::Viewer(QWidget *parent)  : QGLViewer(parent), flatShading_(false), frontFace_(true), pointSize_(0.1), showScalar_(false)
{
	scalarRange_[0] = 0; scalarRange_[1] = 1;

	// The STEREO action is disabled
	setShortcut(STEREO, 0);
	// Add custom key description (see keyPressEvent).	
	setKeyDescription(Qt::CTRL+Qt::Key_F, "Toggles front face display");
	setKeyDescription(Qt::Key_S, "Toggles flat shading display");
	setKeyDescription(Qt::Key_W, "Toggles wire frame display");
	setKeyDescription(Qt::Key_Equal, "Increase picked point size");
	setKeyDescription(Qt::Key_Minus, "Decrease picked point size");
	
	restoreStateFromFile();  // Restore previous Viewer state.  
	help();// Opens help window
}
// init opengl
void Viewer::init()
{
	setBackgroundColor(QColor(255,255,255));
	//setManipulatedFrame(new ManipulatedFrame());
}
void drawSpiral()
{
	const float nbSteps = 200.0;

	glBegin(GL_QUAD_STRIP);
	for (int i=0; i<nbSteps; ++i)
	{
		const float ratio = i/nbSteps;
		const float angle = 21.0*ratio;
		const float c = cos(angle);
		const float s = sin(angle);
		const float r1 = 1.0 - 0.8f*ratio;
		const float r2 = 0.8f - 0.8f*ratio;
		const float alt = ratio - 0.5f;
		const float nor = 0.5f;
		const float up = sqrt(1.0-nor*nor);
		glColor3f(1.0-ratio, 0.2f , ratio);
		glNormal3f(nor*c, up, nor*s);
		glVertex3f(r1*c, alt, r1*s);
		glVertex3f(r2*c, alt+0.05f, r2*s);
	}
	glEnd();
}
// from hsv(240/360, 1, 1) to hsv(0/360, 1, 1)
QColor getColor(double scalar, double* scalarRange)
{
	double val = 1- (scalar-scalarRange[0])/(scalarRange[1] - scalarRange[0]);
	return QColor::fromHsvF(val*2.0/3.0, 1, 1);
}
void drawFaces(Polyhedron &mesh, bool showScalar, double* scalarRange)
{
	glBegin(GL_TRIANGLES);//GL_POLYGON
		for(Polyhedron::Facet_iterator pFacet =  mesh.facets_begin(); 
			  pFacet != mesh.facets_end(); ++pFacet)
		{
			//::glNormal3d(pFacet->normal_.x(), pFacet->normal_.y(), pFacet->normal_.z()); // normal per face

			// revolve around current face to get vertices
			Polyhedron::Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
			do
			{
				Polyhedron::Vertex_handle vh =  pHalfedge->vertex();
				if (showScalar){
					QColor c = getColor(vh->scalar_, scalarRange);
					::glColor3d( c.redF(), c.greenF(), c.blueF() );
				}
				// one normal per vertex
				Kernel::Vector_3& n = vh->normal_;
				::glNormal3f( n[0], n[1], n[2]);
				const Kernel::Point_3& point  = vh->point();		  
				::glVertex3d(point[0],point[1],point[2]);
			}
			while(++pHalfedge != pFacet->facet_begin());
		}
	glEnd(); // end polygon assembly
	//glFlush();
}
void drawVertices(std::list<Polyhedron::Vertex_iterator> &pickedVertices, double pointSize)
{
	// save previous state
	GLint polygonMode;
	glGetIntegerv(GL_POLYGON_MODE, &polygonMode);

	//
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);	
	GLUquadricObj* pQuadric = gluNewQuadric();	
	Polyhedron::Vertex_iterator pv;	
	for (std::list<Polyhedron::Vertex_iterator>::iterator it = pickedVertices.begin(); it != pickedVertices.end(); ++it)
	{
		::glPushMatrix();
		pv = *it;
		::glTranslatef( (float) pv->point().x(),	(float) pv->point().y(), (float) pv->point().z());
		::gluSphere(pQuadric, pointSize, 24, 24); 
		::glPopMatrix();
	}
	gluDeleteQuadric(pQuadric);


	// restore state
	glPolygonMode(GL_FRONT_AND_BACK, polygonMode);
}
void Viewer::drawSelectionRectangle() const
{
  startScreenCoordinatesSystem();
  glDisable(GL_LIGHTING);
  glEnable(GL_BLEND);

  glColor4f(0.0, 0.0, 0.3f, 0.3f);
  glBegin(GL_QUADS);
  glVertex2i(rectangle_.left(),  rectangle_.top());
  glVertex2i(rectangle_.right(), rectangle_.top());
  glVertex2i(rectangle_.right(), rectangle_.bottom());
  glVertex2i(rectangle_.left(),  rectangle_.bottom());
  glEnd();

  glLineWidth(2.0);
  glColor4f(0.4f, 0.4f, 0.5f, 0.5f);
  glBegin(GL_LINE_LOOP);
  glVertex2i(rectangle_.left(),  rectangle_.top());
  glVertex2i(rectangle_.right(), rectangle_.top());
  glVertex2i(rectangle_.right(), rectangle_.bottom());
  glVertex2i(rectangle_.left(),  rectangle_.bottom());
  glEnd();

  glDisable(GL_BLEND);
  glEnable(GL_LIGHTING);
  stopScreenCoordinatesSystem();
}
void Viewer::draw()
{
	if (mesh_.empty())	
	{
		drawSpiral();
		return;
	}

	if (selectionMode_ == ADD)
		drawSelectionRectangle();
	
	glColor3d(1, 0, 0);
	drawVertices(pickedVertices_, pointSize_);

	/////////////////////////////////////////////////////////////////////////
	//draw faces
	glColor3d(0.2, 0.4, 0.4);//set mesh color
    //glPolygonMode(GL_BACK, GL_FILL);  // 设置反面为线形模式
	drawFaces(mesh_, showScalar_, scalarRange_);
}
void Viewer::drawWithNames()
{
	GLUquadricObj* pQuadric = gluNewQuadric();	

	for(Polyhedron::Vertex_iterator pv = mesh_.vertices_begin(); pv !=  mesh_.vertices_end(); ++pv)
	{
		::glPushName(pv->index_);

		::glPushMatrix();
		double radius = average_edge_length_around(pv);
		::glTranslatef( (float) pv->point().x(),	(float) pv->point().y(), (float) pv->point().z());
		::gluSphere(pQuadric, 0.2*radius, 24, 24); 
		::glPopMatrix();

		::glPopName();
	}
	gluDeleteQuadric(pQuadric);
}
void Viewer::mousePressEvent(QMouseEvent* e)
{
    if(e->buttons() & Qt::RightButton &&  e->modifiers() & Qt::AltModifier )
    {
	    beDraging_ = true; 

		rectangle_ = rectangle_.normalized();
		setSelectRegionWidth(rectangle_.width());
		setSelectRegionHeight(rectangle_.height());
		select(rectangle_.center());

		if ( current_picked_vertex != 0)
		{
			move_p = current_picked_vertex->point();
			std::stringstream ss;
			ss<<"move_p"<<move_p.x()<<" "<<move_p.y()<<" "<<move_p.z()<<std::endl;
			QString str(ss.str().c_str());
			emit vertsPicked(str);
		}
    }
	else
		beDraging_ = false;

	rectangle_ = QRect(e->pos(), e->pos());
	if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::ShiftModifier))
		selectionMode_ = ADD;
	else
		QGLViewer::mousePressEvent(e);
}

void Viewer::mouseMoveEvent(QMouseEvent *e)
{	
	if (beDraging_ && 
		e->buttons() & Qt::RightButton &&
		e->modifiers() & Qt::AltModifier)
	{		
		if (current_picked_vertex!=0)
		{
			Point3 p3 = current_picked_vertex->point();
			Vec s3(p3.x(), p3.y(), p3.z());
			Vec s2 = camera()->projectedCoordinatesOf(s3);
			s2.x = e->pos().x();
			s2.y = e->pos().y();

			Vec ns3 = camera()->unprojectedCoordinatesOf(s2);
			current_picked_vertex->point() = Point3(ns3.x, ns3.y, ns3.z);
			updateGL();
		}
	}
	else if( selectionMode_ == ADD && e->buttons() & Qt::LeftButton &&
		e->modifiers() & Qt::ShiftModifier)
	{
		rectangle_.setBottomRight(e->pos());
		updateGL();
	}
	else
		QGLViewer::mouseMoveEvent(e);
}

void Viewer::mouseReleaseEvent(QMouseEvent* e)
{
	if (selectionMode_ == ADD)
	{
		rectangle_ = rectangle_.normalized();
		setSelectRegionWidth(rectangle_.width());
		setSelectRegionHeight(rectangle_.height());
		select(rectangle_.center());
		updateGL();
    }
	else
		QGLViewer::mouseReleaseEvent(e);
}

void Viewer::endSelection(const QPoint& point)
{
	glFlush();
	GLint nbHits = glRenderMode(GL_RENDER);
	
	// Interpret results : each object created 4 values in the selectBuffer().
	// (selectBuffer())[4*i+3] is the id pushed on the stack.
	for (int i=0; i<nbHits; ++i)
	{
		addIdToSelection((selectBuffer())[4*i+3]);
	}

    selectionMode_ = NONE;
//	QGLViewer::endSelection(point);

}

void Viewer::clearSelectedPoints()
{
	this->pickedVertices_.clear();
	vertPickedStatus_.assign(vertPickedStatus_.size(), false);
}
void Viewer::invertSelectedPoints()
{
	if (mesh_.empty())	
		return;


	this->pickedVertices_.clear();

	int i = 0;
	for(Polyhedron::Vertex_iterator vi =  mesh_.vertices_begin(); vi != mesh_.vertices_end(); ++vi, ++i)
	{
		if (vertPickedStatus_[i])
		{
			vertPickedStatus_[i] = false;
		}
		else
		{
			vertPickedStatus_[i] = true;
			pickedVertices_.push_back(vi);
		}
	}
}
void Viewer::saveSelectedPoints()
{
	std::vector<float> arr;
	arr.push_back(1.0);	arr.push_back(2.0);	arr.push_back(3.0);	arr.push_back(4.0);
	std::ofstream ofs("output.txt");
	for (std::list<Polyhedron::Vertex_iterator>::iterator it = this->pickedVertices_.begin(); it != pickedVertices_.end(); ++it)
	{
		Polyhedron::Vertex_iterator vi = *it;
		ofs << vi->index_ << std::endl;
	}
	ofs << std::endl;
	ofs.close();
	
}
void Viewer::addIdToSelection(int id)
{
	for(Polyhedron::Vertex_iterator vi = mesh_.vertices_begin(); vi != mesh_.vertices_end(); ++vi)
	{		
		if ( vi->index_ == id)
		{		
			vertPickedStatus_[id] = true;
			//std::cout << "Picked vertex name is " << id << std::endl;
					
			std::list<Polyhedron::Vertex_iterator>::iterator it = std::find(pickedVertices_.begin(), pickedVertices_.end(), vi);
			if ( it == pickedVertices_.end())
			{
				//std::cout << "vertex: " << vi->point() << " added to picked vertices" << std::endl;
				pickedVertices_.push_back(vi);	
				current_picked_vertex = vi;

				//ManipulatedFrame* mf = manipulatedFrame();
				//mf->setPosition(vi->point().x(),vi->point().y(),vi->point().z());
			}
			else
			{
				if (beDraging_)
				{
					current_picked_vertex = vi;
				}
				else
				{
					//std::cout << "vertex: " << vi->point() << " removed frompicked vertices" << std::endl;
					pickedVertices_.erase(it);
					vertPickedStatus_[id] = false;
				}
			}
			/*std::stringstream ss;	
			ss << vi->index_ << "(" << vi->point().x() << ", " << vi->point().y() << ", "<< vi->point().z() << ") ";
			QString str(ss.str().c_str());
			emit vertsPicked(str);*/
			break;
		}
	}
}
//void Viewer::postSelection(const QPoint& point)
//{
//	int choose = selectedName();
//	if ( choose == -1)
//	{
//		picked_ = false;
//		std::cout << "No vertex selected under pixel " << point.x()  << ","  << point.y() << std::endl;
//		return;
//	}
//
//	for(Polyhedron::Vertex_iterator vi = mesh_.vertices_begin(); vi != mesh_.vertices_end(); ++vi)
//	{		
//		if ( vi->index_ == choose)
//		{				
//			std::cout << "Picked vertex name is " << choose << std::endl;
//					
//			std::list<Polyhedron::Vertex_iterator>::iterator it = std::find(pickedVertices_.begin(), pickedVertices_.end(), vi);
//			if ( it == pickedVertices_.end())
//			{
//				std::cout << "vertex: " << vi->point() << " added to picked vertices" << std::endl;
//				pickedVertices_.push_back(vi);		
//
//				//ManipulatedFrame* mf = manipulatedFrame();
//				//mf->setPosition(vi->point().x(),vi->point().y(),vi->point().z());
//				picked_ = true;
//			}
//			else
//			{
//				std::cout << "vertex: " << vi->point() << " removed frompicked vertices" << std::endl;
//				pickedVertices_.erase(it);
//				picked_ = false;
//			}
//			std::stringstream ss;	
//			ss << vi->index_ << "(" << vi->point().x() << ", " << vi->point().y() << ", "<< vi->point().z() << ") ";
//			QString str(ss.str().c_str());
//			emit vertsPicked(str);
//			break;
//		}
//	}
//}
void Viewer::keyPressEvent(QKeyEvent *e)
{
	// Get event modifiers key
	const Qt::KeyboardModifiers modifiers = e->modifiers();

	bool handled = false;
	if ((e->key()==Qt::Key_W) && (modifiers==Qt::NoModifier)) // w
	{
		wireframe_ = !wireframe_;
		if (wireframe_)
			glPolygonMode(GL_FRONT, GL_LINE);
		else
			glPolygonMode(GL_FRONT, GL_FILL);
		handled = true;
		updateGL();
	}
	else
	if ((e->key()==Qt::Key_F) && (modifiers==Qt::CTRL)) // ctrl+f
	{		
		invertFace();
		handled = true;		
	}
	else
	if ((e->key()==Qt::Key_S) && (modifiers==Qt::NoModifier)) // s
	{
		flatShading_ =  !flatShading_;
		if (flatShading_)
			glShadeModel(GL_FLAT);
		else
			glShadeModel(GL_SMOOTH);
		handled = true;
		updateGL();
	}
	else
	if((e->key()==Qt::Key_Equal) && (modifiers==Qt::NoModifier))
	{
		pointSize_ = pointSize_ * 1.2;
		handled = true;
		updateGL();
	}
	else
	if((e->key()==Qt::Key_Minus) && (modifiers==Qt::NoModifier))
	{
		pointSize_ = pointSize_ * 0.8;
		handled = true;
		updateGL();
	}
  // ... and so on with other else/if blocks.

  if (!handled)
    QGLViewer::keyPressEvent(e);
}
void Viewer::invertFace()
{
	frontFace_ = !frontFace_;
	if (frontFace_)
		glFrontFace(GL_CCW);
	else
		glFrontFace(GL_CW);
	updateGL();
}
void Viewer::invertNormal()
{
	for ( Polyhedron::Facet_iterator fi = mesh_.facets_begin(); fi != mesh_.facets_end(); ++fi)
	{
		fi->normal_ = - fi->normal_;
	}
	for ( Polyhedron::Vertex_iterator vi = mesh_.vertices_begin(); vi != mesh_.vertices_end(); ++vi)
	{
		vi->normal_ = - vi->normal_;
	}
	updateGL();
}

void Viewer::showWhole()
{
	qglviewer::Vec minv(0,0,0), maxv(1,1,1);
	compute_bounding_box(mesh_, minv, maxv);	
	camera()->fitBoundingBox(minv, maxv);
	 	 
	updateGL();
}
void Viewer::showScalar()
{
	showScalar_ = !showScalar_;
	updateGL();
}
QString Viewer::helpString() const
{
  QString text("<h2>S i m p l e V i e w e r</h2>");
  text += "Use the mouse to move the camera around the object. ";
  text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
  text += "Left and middleb uttons pressed together rotate around the camera view direction axis<br><br>";
  text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
  text += "Simply press the function key again to restore it. Several keyFrames define a ";
  text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
  text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
  text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
  text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
  text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
  text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
  text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
  text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
  text += "Press <b>Escape</b> to exit the Viewer.";
  return text;
}
