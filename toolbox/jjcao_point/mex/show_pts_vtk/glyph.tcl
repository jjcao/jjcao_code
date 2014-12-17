package require vtk
package require vtkinteraction

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
  renWin AddRenderer ren1
vtkRenderWindowInteractor iren
  iren SetRenderWindow renWin

# create the piplinee, ball and spikes
#
vtkSphereSource sphere
  sphere SetThetaResolution 10 ; sphere SetPhiResolution 10
vtkPlaneSource plane
  plane SetResolution 10 10;
vtkConeSource cone
  cone SetResolution 5; cone SetHeight 0.1 ; cone SetRadius 0.1
vtkDiskSource disk
  disk SetInnerRadius 0; disk SetOuterRadius 0.1
vtkTransform trans
  trans RotateY 90
vtkTransformFilter diskTF
  diskTF SetInputConnection [disk GetOutputPort]
  diskTF SetTransform trans

# glyph for sphere and cone
vtkGlyph3D glyphSC
  glyphSC SetInputConnection [sphere GetOutputPort]
  glyphSC SetSource [cone GetOutput]
  glyphSC SetVectorModeToUseNormal 
vtkPolyDataMapper mapperSC
  mapperSC SetInputConnection [glyphSC GetOutputPort]
vtkActor actorSC
  actorSC SetMapper mapperSC

# glyph for plane and cone
vtkGlyph3D glyphPC
  glyphPC SetInputConnection [plane GetOutputPort]
  glyphPC SetSource [cone GetOutput]
  glyphPC SetVectorModeToUseNormal 
vtkPolyDataMapper mapperPC
  mapperPC SetInputConnection [glyphPC GetOutputPort]
vtkActor actorPC
  actorPC SetMapper mapperPC

# glyph for sphere and disk
vtkGlyph3D glyphSD
  glyphSD SetInputConnection [sphere GetOutputPort]
  glyphSD SetSource [diskTF GetOutput]
  glyphSD SetVectorModeToUseNormal 
vtkPolyDataMapper mapperSD
  mapperSD SetInputConnection [glyphSD GetOutputPort]
vtkActor actorSD
  actorSD SetMapper mapperSD

# glyph for plane and disk
vtkGlyph3D glyphPD
  glyphPD SetInputConnection [plane GetOutputPort]
  glyphPD SetSource [diskTF GetOutput]
  glyphPD SetVectorModeToUseNormal 
vtkPolyDataMapper mapperPD
  mapperPD SetInputConnection [glyphPD GetOutputPort]
vtkActor actorPD
  actorPD SetMapper mapperPD

#
actorSC SetPosition -0.8 -0.8 0; ren1 AddActor actorSC;
actorPC SetPosition -0.8  0.8 0; ren1 AddActor actorPC;
actorSD SetPosition  0.8 -0.8 0; ren1 AddActor actorSD;
actorPD SetPosition  0.8  0.8 0; ren1 AddActor actorPD;

renWin SetSize 800 800
renWin Render;
iren AddObserver UserEvent {wm deiconify .vtkInteract}
iren Initialize
wm withdraw .