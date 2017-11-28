#include <vtkSmartPointer.h>
 
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkCubeSource.h>
#include <vtkConeSource.h>
#include <vtkSphereSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkCubeSource.h>

int main(int, char *[])
{
  // Creer un cube.
  vtkCubeSource *cube = vtkCubeSource::New();
  vtkSphereSource *sphere = vtkSphereSource::New();
  //cube->SetXLength(2.);
  //vtkConeSource *cube = vtkConeSource::New();
  cube->SetCenter(0., 0., 0.);
  sphere->SetCenter(1., 1., 1.);
  sphere->SetThetaResolution(100);
  sphere->SetPhiResolution(50);

  // Creer un mapper et un actor.
  vtkPolyDataMapper *cubeMapper = vtkPolyDataMapper::New();
  cubeMapper->SetInputConnection( cube->GetOutputPort() );
  
  vtkActor *cubeActor = vtkActor::New();
  cubeActor->SetMapper( cubeMapper );
 
   // Creer un mapper et un actor.
  vtkPolyDataMapper *sphereMapper = vtkPolyDataMapper::New();
  sphereMapper->SetInputConnection( sphere->GetOutputPort() );
  
  vtkActor *sphereActor = vtkActor::New();
  sphereActor->SetMapper( sphereMapper );
/*
  sphereActor->GetProperty()->SetColor(1,0,0);
  sphereActor->GetProperty()->SetAmbient(0.125);
  sphereActor->GetProperty()->SetDiffuse(0.0);
  sphereActor->GetProperty()->SetSpecular(0.0);
*/
  sphereActor->GetProperty()->SetColor(0.23,1,0);
  sphereActor->GetProperty()->SetAmbient(0.50);
  sphereActor->GetProperty()->SetDiffuse(0.1);
  sphereActor->GetProperty()->SetSpecular(0.40);
  sphereActor->GetProperty()->SetOpacity(0.5);
  sphereActor->AddPosition(1.25, 1.25, 0);

  // Creer un renderer, un render window et un render window interactor
  vtkRenderer *ren1= vtkRenderer::New();
  // ajouter les actors à la scène
  ren1->AddActor( cubeActor );
  ren1->AddActor( sphereActor );
  ren1->SetBackground( 0.1, 0.2, 0.4 );
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer( ren1 );
  renWin->SetSize( 300, 300 );
  
 
 
  // Render et interact
  int i;
  for (i = 0; i < 360; ++i)
  {
    // render the image
    renWin->Render();
    // rotate the active camera by one degree
    ren1->GetActiveCamera()->Azimuth( 1 );
  }
  
  cube->Delete();
  sphere->Delete();
  cubeMapper->Delete();
  sphereMapper->Delete();
  cubeActor->Delete();
  sphereActor->Delete();
  ren1->Delete();
  renWin->Delete();
    
 
  return EXIT_SUCCESS;
}
