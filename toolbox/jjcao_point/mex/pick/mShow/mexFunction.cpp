/*=================================================================
Show\Pick\Draw points of surface in MATLAB with VTK

Hui Wang & Shengfa Wang 

intput->(V,F);
output->I---the pick points index 

example: I = mShow(V,F);
*=================================================================*/

#include "stdafx.h"
#include "cgalPointToVtk.h"

#include <math.h>
#include <mex.h>
//
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkProperty.h>
#include <vtkLODActor.h>
#include <vtkSphereSource.h>

#include <vtkPointPicker.h>
#include <vtkCellPicker.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleTrackballActor.h>
#include <vtkCallbackCommand.h>


class vtkPickEvent :
  public vtkCommand
{
public:
 static vtkPickEvent* New(){return new vtkPickEvent;}
 void Delete(){delete this;};
 void SetParameters(vtkRenderer* r){m_ren = r;}
 virtual void Execute(vtkObject * caller, unsigned long l, void *callData)
 {  
 // if(!m_picker) return;
  vtkRenderWindowInteractor* iren = reinterpret_cast<vtkRenderWindowInteractor*>(caller);
  vtkAbstractPicker * picker = iren->GetPicker();
  int *pickPos = iren->GetEventPosition();
  picker->Pick((double)pickPos[0], (double)pickPos[1], 0.0, m_ren);
 }
protected:
 vtkRenderer* m_ren;
};

class vtkKeyPress :
	public vtkCommand
{
public:
	static vtkKeyPress* New(){return new vtkKeyPress;}
	void Delete(){delete this;};
	virtual void Execute(vtkObject * caller, unsigned long l, void *callData)
	{  
		// if(!m_picker) return;
		vtkRenderWindowInteractor* iren = reinterpret_cast<vtkRenderWindowInteractor*>(caller);

		vtkInteractorStyleTrackballActor *style = vtkInteractorStyleTrackballActor::New();
		iren->SetInteractorStyle(style);
		style->Delete();
	}
	
};

class vtkKeyRelease :
	public vtkCommand
{
public:
	static vtkKeyRelease* New(){return new vtkKeyRelease;}
	void Delete(){delete this;};
	virtual void Execute(vtkObject * caller, unsigned long l, void *callData)
	{  
		// if(!m_picker) return;
		vtkRenderWindowInteractor* iren = reinterpret_cast<vtkRenderWindowInteractor*>(caller);

		vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
		iren->SetInteractorStyle(style);
		style->Delete();
	}

};


class vtkPrintPointID :
 public vtkCommand
{
public:
 static vtkPrintPointID* New(){return new vtkPrintPointID;}
 void Delete(){delete this;};
 void SetParameters(vtkRenderer* r){m_ren = r;}
 void SetDate(double *v){V = v;}
 void SetOutputDate(double *pointID){PointID = pointID;}
 void SetLength(int l){L = l;}
 virtual void Execute(vtkObject *caller, unsigned long l, void *callData)
 {
  vtkPointPicker*  ptPicker = reinterpret_cast<vtkPointPicker*>(caller);
//  ptPicker->Pick((double)pickPos[0], (double)pickPos[1], 0.0, g_ren1);//触发了EndPickEvent：在CHoleRepairDlg::OnLandmark()中的g_picker->AddObserver( vtkCommand::EndPickEvent, landmarkObserver );

  //std::string fileName("w.txt");
  //std::ofstream stream(fileName.c_str());
  //stream<<ptPicker->GetPointId() <<"\n";
  //stream.close();

   int numOfVertices=L;
   int Id=ptPicker->GetPointId();
   PointID[Id]=Id + 1;
   Point_3 p(V[Id],V[numOfVertices + Id],V[2 * numOfVertices + Id]);
  double RGB[3] = {0.1 ,0.1,1.0};
  cgalPointToVtk(m_ren,p,0.005,RGB);
 }
protected:
	 vtkRenderer* m_ren;
	 double *V,*PointID;
	 int L;
};
void matlabDataToVtkPolyData(double* V,double* F,int numOfVertices,int numOfFacets,vtkPolyData* polyData)
{
	vtkPoints *points = vtkPoints::New();
	int i;
    
    for(i = 0;i < numOfVertices;i++)
	  points->InsertNextPoint(V[i],V[numOfVertices + i],V[2 * numOfVertices + i]);

    vtkCellArray *polys = vtkCellArray::New();
	for(i = 0; i < numOfFacets;i++)
	{
       polys->InsertNextCell(3);

       polys->InsertCellPoint(F[i] - 1);
	   polys->InsertCellPoint(F[numOfFacets + i] - 1);
	   polys->InsertCellPoint(F[2 * numOfFacets + i] - 1);
	}

	polyData->SetPoints(points);
    polyData->SetPolys(polys);

	points->Delete();
	polys->Delete();
}

void mexFunction( int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] )
{ 
	int numOfVertices,numOfFacets;
	double *V,*F,*PointID;
	//std::vector<int> *PointID;


	//Change data
	numOfVertices = mxGetM(prhs[0]);
	numOfFacets = mxGetM(prhs[1]);

	if(mxGetN(prhs[1]) != 3)
		mexErrMsgTxt("The mesh must be triangle mesh!");
	V = mxGetPr(prhs[0]);
	F = mxGetPr(prhs[1]);

	plhs[0] = mxCreateDoubleMatrix(numOfVertices,1,mxREAL);
	PointID = mxGetPr(plhs[0]);

    //The display of VTK
	vtkPolyData *polyData = vtkPolyData::New();
    matlabDataToVtkPolyData(V,F,numOfVertices,numOfFacets,polyData);

	vtkPolyDataMapper *polyDataMapper = vtkPolyDataMapper::New();
	polyDataMapper->SetInput(polyData);
    
	vtkProperty *proper = vtkProperty::New();
    proper->SetInterpolationToPhong();

	vtkLODActor *actor = vtkLODActor::New();
	actor->SetProperty(proper);
	actor->SetMapper(polyDataMapper);

	vtkRenderer *render = vtkRenderer::New();
	render->SetBackground( 1.0,1.0,1.0);
	render->AddActor(actor);

	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(render);
	
	vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);
	

	vtkPointPicker * pointPicker = vtkPointPicker::New();
    
	//pick points
	//vtkPickEvent * lbp = vtkPickEvent::New();
	//lbp->SetParameters(render);
	//iren->AddObserver(vtkCommand::PickEvent, lbp);
	
	//move actor
	vtkKeyPress * rbp = vtkKeyPress::New();
	iren->AddObserver(vtkCommand::KeyPressEvent, rbp);
	
	//move camera
	//vtkKeyRelease * rbr = vtkKeyRelease::New();
	//iren->AddObserver(vtkCommand::KeyReleaseEvent, rbr);

	vtkPrintPointID *ppID = vtkPrintPointID::New();
	ppID->SetParameters(render);
	ppID->SetDate(V);
	ppID->SetLength(numOfVertices);
	ppID->SetOutputDate(PointID);
	pointPicker->AddObserver(vtkCommand::EndPickEvent, ppID);

	iren->SetPicker(pointPicker); //设置picker
	iren->Initialize();
	iren->Start();


	
	
	polyData->Delete();
	polyDataMapper->Delete();
	proper->Delete();
	actor->Delete();
	render->Delete();
	renWin->Delete();
	iren->Delete();
	 pointPicker->Delete();
	 //lbp->Delete();
	 ppID->Delete();
}
