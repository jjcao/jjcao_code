#pragma warning(disable : 4244 4290 4305 4800 4996)

#include <mex.h>
#pragma comment(lib, "libmex.lib")
#pragma comment(lib, "libmx.lib")

#include <vtkSmartPointer.h>
#define VTK_CREATE(type, name) \
	vtkSmartPointer<type> name = vtkSmartPointer<type>::New()
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCommand.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkLight.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkLODActor.h>
#include <vtkProperty.h>
#include <vtkAxesActor.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkSphereSource.h>
#include <vtkDiskSource.h>
#include <vtkConeSource.h>
#include <vtkGlyph3D.h>
#include <vtkActor.h>
#include <vtkPropAssembly.h>
#include <vtkPointPicker.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

#pragma comment(lib,"vtkCommon.lib ") // for vtkSmartPointer
#pragma comment(lib,"vtkFiltering.lib ") // for *->GetOutputPort()
#pragma comment(lib,"vtkGraphics.lib ") // for vtkConeSource
#pragma comment(lib,"vtkRendering.lib ")// for render
#pragma comment(lib,"vtkHybrid.lib ")// fir vtkAxesActor, widgets

#include <sstream>
#include <vector>
using namespace std;

double *pts(0), *radius(0), *pickPtsIdx(0);
int npts(0), radius_dim(0);
vector<vtkSmartPointer<vtkActor> > pickActors;

VTK_CREATE(vtkPolyData, polyData);
VTK_CREATE(vtkLODActor, actor);
VTK_CREATE(vtkPolyDataMapper, mapper);
VTK_CREATE(vtkSphereSource, pickSource);
VTK_CREATE(vtkPolyDataMapper, pickMapper);
VTK_CREATE(vtkSphereSource, sphereSource);
VTK_CREATE(vtkDiskSource, diskSource);
VTK_CREATE(vtkTransformFilter, diskTF);
vtkSmartPointer<vtkActor> splatAssembly(0);
VTK_CREATE(vtkGlyph3D, splatGlyph);
VTK_CREATE(vtkPolyDataMapper, glyphMapper);
VTK_CREATE(vtkLight, light);
VTK_CREATE(vtkTextActor, textActor);

// display Mode
enum SplatType {POINTS=0, SPHERES, SPLATS, UNKNOWN};
int splatType(POINTS);
double opacity(0.8);
bool tranparencyEnabled(false);

// GUI Variables
double sphereRadius=0.01;
double pickRadiusRadio = 1.6;
double clear_color[3] = { 1, 1, 1 }; //white
int window_height = 800; //1100
int window_width  = 800; //1100

// Colors
//double surfelcolor[4] = {1, .73, .0, 0.5f}; // gold
//double surfelcolor[4] = {0.275, .337, .60, 0.5f}; // blue
double surfelcolor[3] = {1.00, .65, .35}; // orange
double backcolor[3] = {.0, .0, .0};//black
double pointPickColor[3] = {0, 0, 1}; // added by jjcao

void print(const char * msg){
	textActor->SetInput(msg);
}

void matlabDataToVtkPolyData(double*pts, int npts, double* scalars, int scalars_dim, double* normals, double* radius, int radius_dim)
{
	VTK_CREATE(vtkPoints, points);
	points->Allocate(npts);
	VTK_CREATE(vtkCellArray, verts);
	verts->Allocate(verts->EstimateSize(1,npts));
	verts->InsertNextCell(npts);

	for(int i=0; i<npts; ++i){
		int tmp = points->InsertNextPoint(pts[i], pts[npts + i], pts[2*npts + i]);
		verts->InsertCellPoint(tmp);
	}
	double * bb = points->GetBounds();
	double xl = sqrt( (bb[1]-bb[0])*(bb[1]-bb[0])+(bb[3]-bb[2])*(bb[3]-bb[2])+(bb[5]-bb[4])*(bb[5]-bb[4]) );
	sphereRadius = xl*0.01;

	polyData->SetPoints(points);
	polyData->SetVerts(verts);

	if (scalars){
		VTK_CREATE(vtkFloatArray, ss);
		ss->SetNumberOfComponents(3);		ss->SetNumberOfTuples(npts);		ss->SetName("Scalars");
		float* scalar = new float[scalars_dim];
		for (int i=0; i < npts; i++)	{
			for(int j = 0; j < scalars_dim; ++j){
				scalar[j] = scalars[j*npts + i];
			}
			ss->SetTuple(i,scalar);
		}		
		delete[] scalar;
		polyData->GetPointData()->SetScalars(ss);
	}

	if (normals)	{
		VTK_CREATE(vtkFloatArray, ns);
		ns->SetNumberOfComponents(3);		ns->SetNumberOfTuples(npts);		ns->SetName("Normals");
		float normal[3];
		for (int i=0; i < npts; i++)	{
			normal[0] = normals[i]; normal[1] = normals[npts + i]; normal[2] = normals[2*npts + i];
			ns->SetTuple(i,normal);
		}
		polyData->GetPointData()->SetNormals(ns);
	}
}

// if want each splat be an ellipse, which long and short axis vary according to two scalars, we've to use vtkProgrammableGlyphFilter.
// Maybe refer to vtkTensorGlyph, 
// No time. So let it go.
void matlabDataToSplatActor(vtkSmartPointer<vtkPolyData> polyData, double* radius, int radius_dim, vtkRenderer* ren1){		
	splatAssembly = vtkSmartPointer<vtkActor>::New();
	
	sphereSource->SetRadius(sphereRadius);	
	diskSource->SetInnerRadius(0); 	diskSource->SetOuterRadius(sphereRadius);	
	diskSource->SetRadialResolution(1);	diskSource->SetCircumferentialResolution(20);	
	VTK_CREATE(vtkTransform, trans);
	trans->RotateY(90);	
	diskTF->SetInputConnection(diskSource->GetOutputPort());
	diskTF->SetTransform(trans);	

	splatGlyph->SetInput(polyData);
	splatGlyph->SetVectorMode(VTK_USE_NORMAL);	   
	splatGlyph->SetSourceConnection(sphereSource->GetOutputPort());

	glyphMapper->SetInputConnection(splatGlyph->GetOutputPort());
	splatAssembly->SetMapper(glyphMapper);
	splatAssembly->SetPickable(0);
	vtkProperty* prop = splatAssembly->GetProperty();
	prop->SetColor(surfelcolor);
	VTK_CREATE(vtkProperty, bprop);	
	bprop->SetColor(backcolor);
	splatAssembly->SetBackfaceProperty(bprop);

	ren1->AddActor(splatAssembly);	
}
static void keyCallback( vtkObject* object,
					   unsigned long event,
					   void* clientdata,
					   void* vtkNotUsed(calldata) )
{
	vtkInteractorStyleTrackballCamera * style = (vtkInteractorStyleTrackballCamera*)clientdata;
	vtkRenderWindowInteractor* rwi=style->GetInteractor();
	vtkRenderer* ren1 = rwi->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

	vtkProperty* prop = actor->GetProperty();
	stringstream ss;
	switch (rwi->GetKeyCode()){
		case 'v':	// view type switch
			splatType = splatType + 1;
			splatType = splatType%UNKNOWN;
			if (splatType==POINTS){
				actor->SetVisibility(true);	
				if (splatAssembly) splatAssembly->SetVisibility(false);		
			}else{
				actor->SetVisibility(false);	
				if (!splatAssembly) matlabDataToSplatActor(polyData, radius, radius_dim, ren1);
				splatAssembly->SetVisibility(true);
				switch(splatType){
					case SPHERES:
						splatGlyph->SetSourceConnection(sphereSource->GetOutputPort());
						break;
					case SPLATS:
						if (polyData->GetPointData()->GetNormals())	{
							splatGlyph->SetSourceConnection(diskTF->GetOutputPort());
						}						
						break;
				}
			}
			break;
		case 't':// transparency 
			if (actor->GetProperty()->GetOpacity() < 1){
				actor->GetProperty()->SetOpacity(opacity);
			}else{
				actor->GetProperty()->SetOpacity(1.0);
			}			
			break;	
		case '-' : //transparency
			opacity = opacity * 0.9;
			if (opacity < 0.0) opacity = 0.0;
			actor->GetProperty()->SetOpacity(opacity);
			if (splatAssembly) splatAssembly->GetProperty()->SetOpacity(opacity);
			ss << "- transparency " << opacity << endl;
			print(ss.str().c_str());
			break;
		case '=' : 
			opacity = opacity * 1.1;
			if (opacity > 1.0) opacity = 1.0;
			actor->GetProperty()->SetOpacity(opacity);
			if (splatAssembly) splatAssembly->GetProperty()->SetOpacity(opacity);
			ss << "+ transparency " << opacity << endl;
			print(ss.str().c_str());
			break;
		case 's': // scalar switch
			mapper->SetScalarVisibility(!mapper->GetScalarVisibility());
			glyphMapper->SetScalarVisibility(!glyphMapper->GetScalarVisibility());
		case 'a': // scalar splat radius
			if ( splatGlyph->GetScaleMode()==VTK_DATA_SCALING_OFF){
				splatGlyph->SetScaleMode(VTK_SCALE_BY_SCALAR);
				splatGlyph->SetScaleFactor(0.25);	
			}else{
				splatGlyph->SetScaleMode(VTK_DATA_SCALING_OFF);
			}
			mapper->SetScalarVisibility(!mapper->GetScalarVisibility());
			glyphMapper->SetScalarVisibility(!glyphMapper->GetScalarVisibility());
		case 'c':// light follow camera
			rwi->SetLightFollowCamera( !rwi->GetLightFollowCamera() );
			ss << "light follow camera: " << rwi->GetLightFollowCamera() << endl;	
			print(ss.str().c_str());
			break;
		case 'w':
			prop->SetEdgeVisibility(!prop->GetEdgeVisibility());
			break;
		case 'l' : // light switch
		case 'L' :					
			light->SetSwitch(!light->GetSwitch());
			ss << "Lighting switch to " << light->GetSwitch() << endl;
			print(ss.str().c_str());
			break;
		case 'f' : // shading
			prop->SetInterpolation(VTK_FLAT);
			ss << "shading VTK_FLAT" << endl;
			print(ss.str().c_str());
			break;
		case 'g' : 
			prop->SetInterpolation(VTK_GOURAUD);
			ss << "shading VTK_GOURAUD" << endl;
			print(ss.str().c_str());
			break;
		case 'p' : 
			prop->SetInterpolation(VTK_PHONG);
			ss << "shading VTK_PHONG" << endl;
			print(ss.str().c_str());
			break;
		//case '-' : //Intensity of light
		//	light->SetIntensity( light->GetIntensity()*0.9 );
		//	ss << "- Intensity of light " << light->GetIntensity() << endl;
		//	print(ss.str().c_str());
		//	break;
		//case '=' : 
		//	light->SetIntensity( light->GetIntensity()*1.1 );
		//	ss << "+ Intensity of light " << light->GetIntensity() << endl;
		//	print(ss.str().c_str());
		//	break;
		case ',' : //point size
			prop->SetPointSize( prop->GetPointSize()*0.9 );
			ss << "- point size" << prop->GetPointSize() << endl;		
			print(ss.str().c_str());
			break;
		case '.' : 
			prop->SetPointSize( prop->GetPointSize()*1.1 );
			ss << "+ point size " << prop->GetPointSize() << endl;
			print(ss.str().c_str());
			break;
		case '[' : //picked point size
			sphereRadius = sphereRadius*0.9;
			pickSource->SetRadius(sphereRadius*pickRadiusRadio);
			sphereSource->SetRadius(sphereRadius);
			diskSource->SetOuterRadius(sphereRadius);	
			ss << "- picked point size" << sphereRadius<< endl;		
			print(ss.str().c_str());
			break;
		case ']' : 
			sphereRadius = sphereRadius*1.1;
			pickSource->SetRadius(sphereRadius*pickRadiusRadio);
			sphereSource->SetRadius(sphereRadius);
			diskSource->SetOuterRadius(sphereRadius);	
			ss << "+ picked point size " << sphereRadius << endl;
			print(ss.str().c_str());
			break;
		case 'u': // clear picked points
			for(int i=0; i<npts; ++i){
				bool tmp = pickPtsIdx[i];
				if(tmp) ren1->RemoveActor(pickActors[i]);
				pickActors[i] =0;
				pickPtsIdx[i] = 0.0;
				//sphereAssembly->AddPart(sa);
			}
			break;
	}	
	style->OnKeyPress();//OnKeyDown

	rwi->Render();
}
void pick(int x, int y, vtkRenderer* ren1, vtkRenderWindowInteractor* rwi){
	vtkPointPicker *picker = vtkPointPicker::SafeDownCast(rwi->GetPicker());
	picker->Pick((double)x, (double)y, 0.0, ren1);

	double* globalCoordinate = picker->GetPickPosition(); 
	stringstream ss;
	ss << "CTRL pressed: " << "Picked id: " << picker->GetPointId() << endl;
	ss << "Global position: " << globalCoordinate[0] << ", " << globalCoordinate[1] << ", " << globalCoordinate[2] << endl;
	print(ss.str().c_str());

	int id = picker->GetPointId();
	if (id>-1){
		bool bPick = !pickPtsIdx[id];
		pickPtsIdx[id] = bPick;
		if (bPick)	{
			VTK_CREATE(vtkActor, sa);	
			sa->SetMapper(pickMapper);
			vtkProperty *pp = sa->GetProperty();
			pp->SetColor(pointPickColor);
			sa->SetPosition(pts[id], pts[npts + id], pts[2*npts + id]);
			sa->SetPickable(0);
			pickActors[id] = sa;
			ren1->AddActor(sa);
		}else{
			ren1->RemoveActor(pickActors[id]);
			pickActors[id] = 0;
			//pickActors[id]->SetVisibility(! pickActors[id]->GetVisibility() );
		}
	}
}
static void mouseCallback( vtkObject* object,
						unsigned long event,
						void* clientdata,
						void* vtkNotUsed(calldata) )
{
	vtkInteractorStyleTrackballCamera * style = (vtkInteractorStyleTrackballCamera*)clientdata;
	vtkRenderWindowInteractor* rwi=style->GetInteractor();
	vtkRenderer* ren1 = rwi->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

	vtkProperty* prop = actor->GetProperty();
	stringstream ss;
	switch( event ){	
		case vtkCommand::MouseMoveEvent:
			if(rwi->GetAltKey() && light->GetLightType()==VTK_LIGHT_TYPE_SCENE_LIGHT){
				vtkCamera* camera = ren1->GetActiveCamera();
				light->SetPosition(camera->GetPosition());
				light->SetFocalPoint(camera->GetFocalPoint());
			}
			style->OnMouseMove();
			break;
		//end of vtkCommand::MouseMoveEvent
		case vtkCommand::LeftButtonPressEvent:
			int x,y;
			rwi->GetEventPosition(x,y);
			//ss << "left button clicked: " << x << ", " << y << endl;
			if (rwi->GetControlKey()) {		
				pick(x,y, ren1, rwi);
			}
			style->OnLeftButtonDown();
			break;
		   // end vtkCommand::LeftButtonPressEvent
	}	
	rwi->Render();
}
void initScene(vtkSmartPointer<vtkRenderer> render, vtkSmartPointer<vtkRenderWindow> renWin, 
			   vtkSmartPointer<vtkRenderWindowInteractor> rwi, vtkSmartPointer<vtkInteractorStyleTrackballCamera> style){
   splatAssembly = 0;
   tranparencyEnabled = false;

	render->SetBackground(clear_color);		
	renWin->GetRenderers()->RemoveAllItems();
	renWin->AddRenderer(render);
	renWin->SetSize(window_width, window_height );	

	textActor->SetInput("Init");
	textActor->SetDisplayPosition(window_width*0.5, window_height*0.05);
	vtkTextProperty *tprop = textActor->GetTextProperty();
	tprop->SetColor(0.0,0.0,1.0);
	tprop->SetJustificationToCentered();
	render->AddViewProp(textActor);

	VTK_CREATE(vtkCallbackCommand, keyCommand);
	VTK_CREATE(vtkCallbackCommand, mouseCommand);
	keyCommand->SetClientData(style);
	keyCommand->SetCallback(keyCallback);
	mouseCommand->SetClientData(style);
	mouseCommand->SetCallback(mouseCallback);
	style->AddObserver(vtkCommand::KeyPressEvent, keyCommand);
	style->AddObserver(vtkCommand::MouseMoveEvent, mouseCommand);	
	style->AddObserver(vtkCommand::LeftButtonPressEvent, mouseCommand);	
	VTK_CREATE(vtkPointPicker, picker);
	rwi->SetPicker(picker);
	//picker->SetTolerance	(0.025);//0.025 is the defaul fraction of rendering window size. (Rendering window size is measured across diagonal.)

	// for drawing picked points
	pickSource->SetRadius(sphereRadius*pickRadiusRadio);	
	pickMapper->SetInputConnection(pickSource->GetOutputPort());	
	pickActors.clear();
	pickActors.reserve(npts);
	for(int i=0; i<npts; ++i){
		pickActors.push_back(0);
		pickPtsIdx[i] = 0.0;
	}
}
void mexFunction(int nlhs, mxArray*plhs[], 
				 int nrhs, const mxArray*prhs[])
{
	/* Check for proper number of arguments. */
	if (nrhs < 1) {
		mexErrMsgTxt("At least one parameter required.");
	}

	/* Handle parameters and outputs. */
	int scalars_dim(0);
	double *scalars(0), *normals(0);

	npts = mxGetM(prhs[0]);
	pts = mxGetPr(prhs[0]);
	if (nrhs > 1) {
		scalars = mxGetPr(prhs[1]);
		scalars_dim = mxGetN(prhs[1]);
	}
	if (nrhs > 2) {
		normals = mxGetPr(prhs[2]);
	}
	if (nrhs > 3) {
		radius = mxGetPr(prhs[3]);
		radius_dim = mxGetN(prhs[3]);
	}

	plhs[0] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	pickPtsIdx = mxGetPr(plhs[0]);

	/* initScene. */	
	VTK_CREATE(vtkRenderer, render);//have to be a local variable, so the mex can exist normally
	VTK_CREATE(vtkRenderWindow, renWin);//have to be a local variable, so the mex can exist normally
	VTK_CREATE(vtkRenderWindowInteractor, rwi);//have to be a local variable, so the mex can exist normally
	VTK_CREATE(vtkInteractorStyleTrackballCamera, style);	//if it is not a local variable, after GUI appeared, there are some delay.
	initScene(render, renWin, rwi, style);

	 /******************** prepare data. ****************************/
	matlabDataToVtkPolyData(pts, npts, scalars, scalars_dim, normals, radius, radius_dim);	
	mapper->SetInput(polyData);
	actor->SetMapper(mapper);
	vtkProperty* prop = actor->GetProperty();
	prop->SetColor(surfelcolor);
	prop->SetPointSize(5);
	render->AddActor(actor);

  /******************** show. ****************************/
   render->ResetCamera();
	rwi->SetRenderWindow(renWin);	
	rwi->SetInteractorStyle(style);	
	vtkCamera* camera = render->GetActiveCamera();
	light->SetPosition(camera->GetPosition());
	light->SetFocalPoint(camera->GetFocalPoint());
	light->SetSwitch(0);
	//light->SetLightType(VTK_LIGHT_TYPE_HEADLIGHT);//default is VTK_LIGHT_TYPE_SCENE_LIGHT
	//light->SetLightType(VTK_LIGHT_TYPE_CAMERA_LIGHT);
	//light->SetAmbientColor(1, 0.1, 0.1); // no use?
	//light->SetSpecularColor(1, 0, 0);// no use?
	render->AddLight(light);
	rwi->Initialize();
	rwi->Start();

	mexPrintf("show_pts_vtk terminated correctly.\n");
}
