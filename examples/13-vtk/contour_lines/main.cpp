// jjcao, 2013

#include <iostream>
#include <vtkQuadric.h>
#include <vtkSampleFunction.h>
#include <vtkExtractVOI.h>
#include <vtkContourFilter.h>
//#include <vtkBandedPolyDataContourFilter.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>

#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyleTrackballCamera.h>

using namespace std;

#pragma comment( lib, "vtkRendering.lib" )
#pragma comment( lib, "vtkGraphics.lib" )
#pragma comment( lib, "vtkFiltering.lib" )
#pragma comment( lib, "vtkImaging.lib" )
#pragma comment( lib, "vtkCommon.lib" )
//#pragma comment( lib, "vtksys.lib" )

int main()
{
	vtkQuadric *quadric = vtkQuadric::New();
	quadric->SetCoefficients(.5, 1, .2, 0, .1, 0, 0, .2, 0, 0);

	vtkSampleFunction *sample = vtkSampleFunction::New();
	sample->SetSampleDimensions(30,30,30);
	sample->SetImplicitFunction(quadric);
	sample->ComputeNormalsOff();
  
	vtkExtractVOI *extract = vtkExtractVOI::New();
	extract->SetInputConnection(sample->GetOutputPort());
	extract->SetVOI(0, 29, 0, 29, 15, 15);
	extract->SetSampleRate(1, 2, 3);

	vtkContourFilter *contours = vtkContourFilter::New();
	//vtkBandedPolyDataContourFilter *contours = vtkBandedPolyDataContourFilter::New();
	//contours->GenerateContourEdgesOn();
	contours->SetInputConnection(extract->GetOutputPort());
	contours->GenerateValues(3, 0.0, 1.2);
  
  
	vtkPolyDataMapper *contMapper = vtkPolyDataMapper::New();
	contMapper->SetInputConnection(contours->GetOutputPort());
	//contMapper->SetInput(contours->GetContourEdgesOutput());
  
	contMapper->SetScalarRange(0.0, 1.2);

	vtkActor *contActor = vtkActor::New();
	contActor->SetMapper( contMapper );

	vtkRenderer *ren1= vtkRenderer::New();
	ren1->AddActor( contActor );
	ren1->SetBackground( 0.1, 0.2, 0.4 );


	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer( ren1 );
	renWin->SetSize( 300, 300 );

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	vtkInteractorStyleTrackballCamera *style = 
	vtkInteractorStyleTrackballCamera::New();
	iren->SetInteractorStyle(style);

	//
	iren->Initialize();
	iren->Start();
  
	vtkPolyData* ct = contours->GetOutput();
	//vtkPolyData* ct = contours->GetContourEdgesOutput();
	vtkFloatArray *scalars = (vtkFloatArray *)(ct->GetPointData()->GetScalars());
	vtkPoints* points = ct->GetPoints();
	int numP = points->GetNumberOfPoints();
	int numS = scalars->GetSize();
	double scalarRange[2];	
	scalars->GetRange(scalarRange);

	for ( int i = 0; i < numS; ++i)
	{
		cout << scalars->GetValue(i) << ", ";
	}	

	vtkCellArray* lines = ct->GetLines();
	int nLines = ct->GetNumberOfLines();

	quadric->Delete();
	sample->Delete();
	extract->Delete();
	contours->Delete();
	contMapper->Delete();
	contActor->Delete();
	ren1->Delete();
	renWin->Delete();
	iren->Delete();
	style->Delete();

	return 0;
}
