/*=================================================================
Pick points and move them, return the index and the new positions
This is matlab dll program 

Shengfa Wang, Hui Wang 

intput->(V,F);
output->I---the pick points index 
        newPosition---the new positions

example: [I,newPosition] = movePoint(V,F);

Ari. 7, 2010
*=================================================================*/

#include "stdafx.h"
#include "cgalPointToVtk.h"

#include <math.h>
#include <mex.h>

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
#include <vtkInteractorObserver.h>

#include <fstream>
#include <vector>

typedef std::vector<int> ConIndex;
typedef std::vector<double> Point_3D;
typedef std::vector<Point_3D> Points;



//For pick point on the mesh
class PickCommand : public vtkCommand
{
public:
	static PickCommand* New() {return new PickCommand;}
	void Delete() {m_index.clear();delete this;}
	void SetRender(vtkRenderer* render) {m_render = render;}
	int GetNumber() {return m_index.size();}
	void GetIndex(double *index)
	{
      for(int i = 0;i < GetNumber();i++)
		 index[i] = m_index.at(i);
	}
	virtual void Execute(vtkObject *caller, unsigned long l, void *callData)
	{
		vtkPointPicker*  pick = reinterpret_cast<vtkPointPicker*>(caller);
		
		int id = pick->GetPointId();
		if(id != -1 && alreadyPicked(id + 1) == false)
		{
           m_index.push_back(id + 1);
		   double pos[3] = {0,0,0};
		   pick->GetPickPosition(pos);
		   double RGB[3] = {0.1 ,0.1,1.0};
		   cgalPointToVtk(m_render,pos,0.01,RGB);
		}
	}

private:
	vtkRenderer* m_render;
	ConIndex m_index;
	bool alreadyPicked(int id)
	{  
		bool picked = false;

		for(int i = 0;i < GetNumber();i++)
		  if(m_index.at(i) == id)
		  {
			  picked = true;
			  return picked;
		  }

	   return picked;
	}
};


/*
//For move, change the camera style to interact style
class MoveCommand :public vtkCommand
{
public:
	static MoveCommand* New(){return new MoveCommand;}
	void Delete(){delete this;};
	virtual void Execute(vtkObject * caller, unsigned long l, void *callData)
	{  
		vtkRenderWindowInteractor* iren = reinterpret_cast<vtkRenderWindowInteractor*>(caller);
		
	    vtkInteractorStyleTrackballActor *style = vtkInteractorStyleTrackballActor::New();
		iren->SetInteractorStyle(style);

		style->Delete();
	}
};
*/

//For pick new position
class middleRelease :public vtkCommand
{
public:
	static middleRelease* New(){return new middleRelease;}
	void Delete()
	{
		for(int i = 0;i < GetNumber();i++)
			m_newPostions.at(i).clear();
		m_newPostions.clear();

		delete this;
	}
	void SetRender(vtkRenderer *render){m_render = render;}
	int  GetNumber() {return m_newPostions.size();}
	void GetNewPoistions(double* V)
	{
      int i(0),num(GetNumber());

	  for(i = 0; i < num;i++)
	  {
        V[i] = m_newPostions.at(i).at(0);
        V[num + i] = m_newPostions.at(i).at(1);
		V[2 * num + i] = m_newPostions.at(i).at(2);
	  }
	}

	virtual void Execute(vtkObject * caller, unsigned long l, void *callData)
	{  
		vtkRenderWindowInteractor* iren = reinterpret_cast<vtkRenderWindowInteractor*>(caller);

		double pOrigin[3],dOrigin[3]; 
		iren->GetPicker()->GetPickPosition(pOrigin);
		vtkInteractorObserver::ComputeWorldToDisplay(m_render,pOrigin[0],pOrigin[1],pOrigin[2],dOrigin);

		double pNew[4];
		int dNew[2];
		iren->GetEventPosition(dNew);

		vtkInteractorObserver::ComputeDisplayToWorld(m_render,(double)dNew[0],(double)dNew[1],dOrigin[2],pNew); 
		double RGB[3] = {1,0,0};
		cgalPointToVtk(m_render,pNew,0.01,RGB);

		Point_3D P;
		P.resize(3,0);
        P.at(0) = pNew[0];
		P.at(1) = pNew[1];
		P.at(2) = pNew[2];
        m_newPostions.push_back(P);
	    P.clear();
	}

private:
	vtkRenderer *m_render;
	Points m_newPostions;
};


//For key release,change the interact style to camera style
/*
class vtkKeyRelease :public vtkCommand
{
public:
	static vtkKeyRelease* New(){return new vtkKeyRelease;}
	void Delete(){delete this;};
	virtual void Execute(vtkObject * caller, unsigned long l, void *callData)
	{  
		vtkRenderWindowInteractor* iren = reinterpret_cast<vtkRenderWindowInteractor*>(caller);

		vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
		iren->SetInteractorStyle(style);
		style->Delete();
	}

};
*/

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
	double *V,*F,*PointID,*PointNewPosition;


	//Change data
	numOfVertices = mxGetM(prhs[0]);
	numOfFacets = mxGetM(prhs[1]);

	if(mxGetN(prhs[1]) != 3)
		mexErrMsgTxt("The mesh must be triangle mesh!");
	V = mxGetPr(prhs[0]);
	F = mxGetPr(prhs[1]);



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
    


	//Pick point on the mesh
	PickCommand *pickCommand = PickCommand::New(); 
	pickCommand->SetRender(render);
	vtkPointPicker *picker = vtkPointPicker::New();
	picker->AddObserver(vtkCommand::EndPickEvent,pickCommand);
    
    
	//Define the new position by the right button
	middleRelease *midRelease = middleRelease::New();
	midRelease->SetRender(render);
	iren->AddObserver(vtkCommand::RightButtonPressEvent,midRelease);

    

	//Start render and interact 
	iren->SetPicker(picker);
	iren->Initialize();
	iren->Start();
    


	//Return the pick point index and the new positions
	plhs[0] = mxCreateDoubleMatrix(pickCommand->GetNumber(),1,mxREAL);
	PointID = mxGetPr(plhs[0]);
    pickCommand->GetIndex(PointID);
    
	plhs[1] = mxCreateDoubleMatrix(midRelease->GetNumber(),3,mxREAL);
	PointNewPosition = mxGetPr(plhs[1]);
	midRelease->GetNewPoistions(PointNewPosition);


	//Delete data
	polyData->Delete();
	polyDataMapper->Delete();
	proper->Delete();
	actor->Delete();
	render->Delete();
	renWin->Delete();
	iren->Delete();
	picker->Delete();
	pickCommand->Delete();
	midRelease->Delete();
}
