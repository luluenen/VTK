#include <vtkActor.h>
#include <vtkAssembly.h>
#include <vtkCamera.h>
#include <vtkConeSource.h>
#include <vtkContourFilter.h>
#include <vtkCubeSource.h>
#include <vtkDataSetMapper.h>
#include <vtkDICOMImageReader.h>
#include <vtkMarchingCubes.h>
#include <vtkOutlineFilter.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkProgrammableFilter.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkTIFFWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkVectorText.h>
#include <vtkVolume16Reader.h>
#include <vtkWindowToImageFilter.h>


//#include <Math.h>
//#include <String.h>
#include <time.h>
#include <unistd.h>
#include <ctime>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

struct params
{
    vtkSmartPointer<vtkProgrammableFilter> change_coords_filter;
    double Lmax;
    float **X0;
    float **DX;
    double z1,z2;
    vtkSmartPointer<vtkUnstructuredGrid> output ;
    vtkProgrammableFilter* filter;
};

vtkSmartPointer<vtkActor> make_cone(double *x,double H,double R) {
    //## Creates the cone
    double xx[3];
    vtkSmartPointer<vtkConeSource> cone = vtkSmartPointer<vtkConeSource>::New();
    cone->SetHeight(H);
    cone->SetRadius(R);
    cone->SetResolution(20);
    cone->SetDirection(0,1,0);
    
    //## Creates the mapper
    vtkSmartPointer<vtkPolyDataMapper> coneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    coneMapper->SetInputConnection(cone->GetOutputPort() );
    
    //## Creates the actor
    vtkSmartPointer<vtkActor> cone_actor = vtkSmartPointer<vtkActor>::New();
    cone_actor->SetMapper( coneMapper );
    cone_actor->GetProperty()->SetOpacity(0.5);
    xx[0]=x[0];
    xx[1]=0.5*H;
    xx[2]=x[1];
    cone_actor->SetPosition(xx);
    return cone_actor;
}

//## This function is the kernel of the `vtkProgrammableFilter'
//## It changes the coordinates of the worm points
void  change_coords_fun(void* arguments)
{
    params* inputs = static_cast<params*>(arguments);
    
    int numPts;
    double L;
    
    vtkSmartPointer<vtkUnstructuredGrid> input;
    input = inputs->change_coords_filter->GetUnstructuredGridInput();
    
    numPts = input->GetNumberOfPoints();
    vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
    
    L = inputs->z2 - inputs->z1;
    
    //## were fixed.
    double aux;
    aux = (pow(inputs->Lmax,2)-pow(L,2))/2;
    if (aux<0) {
        aux = 0.0;
    }
    double H;
    H = 0.25*sqrt(aux);
    
    float xx[2000];
    float yy[2000];
    float zz[2000];
    
    //## Change the coordinates point by point
    for (int i=0 ; i<=numPts; i++) {
        double *x;
        x=input->GetPoint(i);
        xx[i] = inputs->X0[i][0] + H*inputs->DX[i][0];
        yy[i] = inputs->X0[i][1] + H*inputs->DX[i][1];
        zz[i] = inputs->X0[i][2] + H*inputs->DX[i][2];
        zz[i] = inputs->z1+(L / inputs->Lmax)*zz[i];
        newPts->InsertPoint(i, xx[i], yy[i], zz[i]);
    }
    inputs->output = inputs->change_coords_filter->GetUnstructuredGridOutput();
    inputs->output->CopyStructure(input);
    inputs->output->SetPoints(newPts);
}

int main (int argc, char *argv[])
{
    float **DX;
    float **X0;
    FILE* fichier1;
    double L;
    int nLignes;
    int nColonnes;
    double x,z,R,H;
    double xx[2];
    double z1,z2;
    double Dmin;
    double Rworm;
    double Lmin, Lmax;
    double speed;

    
    nLignes = 3000;
    nColonnes = 3;
    DX = new float* [ nLignes ];
    for (int i=0; i < nLignes; i++)
        DX[i] = new float[ nColonnes ];
    
    X0 = new float* [ nLignes ];
    for (int i=0; i < nLignes; i++)
        X0[i] = new float[ nColonnes ];
    
    vtkSmartPointer<vtkActor> cones[300] ;

    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName("worm.vtk");
    
    vtkSmartPointer<vtkProgrammableFilter> change_coords_filter =  vtkSmartPointer<vtkProgrammableFilter>::New();
    change_coords_filter->SetInputConnection(reader->GetOutputPort());
    
    // A number of cones are added to create a landscape where the worm
    // moves. The following parameters are defined.
    //
    // The size of the landscape (the part of the domain where cones are
    // defined).
    int Lscape = 150;
    // Number of cones
    int ncones = 300;
    // Radius of cones bases are generated randomly in
    // this range
    double Rmaxcone = 5.0;
    double Rmincone = 1;
    // Radius of the worm
    Rworm = 1;
    // Cones are not added in |x|<2*Dmin
    Dmin = Rworm + Rmaxcone;
    

    
    //## All the cones are gathered in an `assembly'
    vtkSmartPointer<vtkAssembly> assembly = vtkSmartPointer<vtkAssembly>::New();
    //## Also the are gathered in a Python array
    //## I think that this may be avoided we know how
    //## to iterate on the assembly.
    srand(time(0));
    
    for (int k=0;k<ncones;k++){
        
        //## Generate a ramdom point in |x|,|z|<Lscape
        x = Lscape*(2.0*(rand()/(double)(RAND_MAX))-1.0);
        //## Exclude a strip along the path of the worm (the z axis)
        if (abs(x)<Dmin) {
            if (x<0) {
                x = -x;
            }
        }
        z = Lscape*(2.0*(rand()/(double)(RAND_MAX))-1.0);
        //## Generate random value of cone radius
        R = Rmincone+(Rmaxcone-Rmincone)*(rand()/(double)(RAND_MAX));
        //## Aspect ratio of cones is fixed
        H = 3*R;
        //## Make the cone
        xx[0] = x;
        xx[1] = z;
        cones[k] = make_cone(xx,H,R);
        //## Add cone to assembly and to `cones' Python array
        assembly->AddPart(cones[k]);
    }
    
    //## The algorithm of movement of the worm is as follows.  The state of
    //## the worm is defined by the positions of the extremes z1,z2.  At
    //## each time one of its extremmities are moving with velocity `speed'.
    //## Which of the extremmities is moving is stored in variable `moving'.
    //## When the length of the worm passes a threshold `Lmax' then we
    //## switch to moving=1 (the tail), and conversely if `L<Lmin' then we
    //## switch to moving the head (moving=2)
    
    //## Initial position of worm extremmities
    z1 = 0.0;
    z2 = 9.0;
    //## Min/Max length of worm
    Lmin = 4.0;
    Lmax = 8.0;
    //## Advance speed of one extremmity at a time
    speed = 3.0;
    //## Which extremmity is advancing now
    int moving;
    moving = 1;
    
    //## Time step
    double Dt;
    Dt = 0.2;
    //## Initial time
    double t;
    t = 0.0;
    //## Average speed of worm
    double mspeed;
    mspeed = speed/2;
    
    //## Position of worm node coordinates at undisturbed position
    //## (i.e. a plain cylinder)
    
    FILE* fichier;
    fichier = fopen("worm.nod.tmp","r");
    int k=0;
    while(fichier != NULL && k < 1701) {
        fscanf(fichier,"%f %f %f",&X0[k][0],&X0[k][1],&X0[k][2]);
        k++;
    }
    
    //    X0 = read_array(file("worm.nod.tmp"));
    //## Displacements for the maximum distorted
    //## position of the worm, i.e. the position
    //## of the coordinates of the worm are X = X0+alpha*DX
    //## with 0<=alpha<=1
    
    fichier1 = fopen("worm.dx-nod.tmp","r");
    k=0;
    while(fichier1 != NULL && k < 1701) {
        fscanf(fichier1,"%f %f %f\n",&DX[k][0],&DX[k][1],&DX[k][2]);
        k++;
    }
    
    params myParams;
    
    myParams.change_coords_filter = change_coords_filter;
    
    
    myParams.Lmax = Lmax;
    
    myParams.X0 = X0;
    myParams.DX = DX;
    myParams.z1 = z1;
    myParams.z2 = z2;
    
    change_coords_filter->SetExecuteMethod(change_coords_fun, &myParams);
    
    //# Create mapper
    
    vtkSmartPointer<vtkDataSetMapper> ugmapper = vtkSmartPointer<vtkDataSetMapper>::New();
    ugmapper->SetInputConnection(change_coords_filter->GetOutputPort());
    ugmapper->ScalarVisibilityOff();
    
    //# Create actor
    vtkSmartPointer<vtkActor> ugactor = vtkSmartPointer<vtkActor>::New();
    ugactor->SetMapper(ugmapper);
    ugactor->GetProperty()->SetDiffuseColor(1, 1, 1);
    
    //# Create the usual rendering stuff.
    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren);
    renWin->SetSize(640,480);
    
    ren->SetBackground(.1, .2, .4);
    ren->AddActor(ugactor);
    ren->AddActor(assembly);
    
    //## This is used to store the frames
    //## for creating a movie
    vtkSmartPointer<vtkWindowToImageFilter> w2i =  vtkSmartPointer<vtkWindowToImageFilter>::New();
    w2i->SetInput(renWin);
    w2i->Update();
    
    //## The TIFF writer
    vtkSmartPointer<vtkTIFFWriter> writer = vtkSmartPointer<vtkTIFFWriter>::New();
    writer->SetInputConnection(w2i->GetOutputPort());
    writer->SetCompressionToJPEG();
    
    //## Initialize camera
    vtkSmartPointer<vtkCamera> cam;
    cam = ren->GetActiveCamera();
    cam->SetFocalPoint(0.,0.,Lmax/2);
    cam->SetViewUp(0.,1.,0.);
    cam->Azimuth(30.0);
    
    //## The flying camera is as follows.
    //## It rotates around the focal point but also the distance to the
    //## object varies sinusoidally between Rcammin and Rcammax
    double Rcammax = 6*Lmax;
    double Rcammin = 2*Lmax;
    //## The speed of rotation of the camera
    double wcam = 0.05;
    //## The speed with which the radius varies
    double wRcam = 1.33*wcam;
    
    
    //# Initial frame
    int frame = 0;
    //## Total nbr of frames
    int nframes = 3000;
    //## Store frames as TIFF images in ./YUV directory
    int save_frames = 0;
    
    //# Loop while: rotating the camera and modify
    //# node coordinates
    int flag = 1;
    while (flag == 1) {
        t += Dt;
        //## time.sleep(0.3)
        if (frame%10==0) {
            printf("frame: %d ",frame);
        }
        
        //## This recomputes the clipping plane Otherwise, since the camera
        //## is rotating, some objects may disappear
        ren->ResetCameraClippingRange();
        change_coords_filter->Modified();
        
        //## Actual length of worm
        L = z2-z1;
        //## Recompute which extreme is moving
        if (L>Lmax) {
            moving = 1;
        }
        if (L<Lmin) {
            moving = 2;
        }
        //## Update positions of extremes
        if (moving == 1) {
            z1 += speed*Dt;
        }
        else {
            z2 += speed*Dt;
        }
        //## Substract mean velocity
        //## We describe the world from a system moving with the
        //## mean velocity of the worm, i.e. mspeed = speed/2
        z1 -= mspeed*Dt;
        z2 -= mspeed*Dt;
        myParams.z1 = z1;
        myParams.z2 = z2;
        
        //## Update camera
        double phi = wcam*t;
        double Rcamav =(Rcammin+Rcammax)/2.0;
        double DRcam = (Rcammin-Rcammax)/2.0;
        double Rcam = Rcamav + DRcam*sin(wRcam*t);
        double Hcam = 0.25*Rcam;
        //cam->SetPosition(Rcam*cos(phi),Hcam,Rcam*sin(phi));
        //cam->SetPosition(100,100,100);
        renWin->Render();
        
        char buffer[100];
        sprintf (buffer, "./YUV/frame.%d.tiff", frame);
        //## Save current frame
        if (save_frames ) {
            w2i->Modified();
            writer->SetFileName(buffer);
            writer->Write();
        }
        
        if (nframes>0 && frame == nframes) {
            flag =0;
        }
        else {
            //## Update the position of the cones.
            //## They move with velocity [0,0, -mspeed].
            //## In order to keep the cones inside the landscape
            //## when a cone exits the landscape (i.e. z_cone<-Lscape)
            //## then it is moved to the same position but at the inlet (i.e.
            //## z_cone = +Lscape
            for (int k =0;k<ncones;k++){
                vtkSmartPointer<vtkActor> cone;
                cone = cones[k];
                cone->AddPosition(0,0,-mspeed*Dt);
                double *xx, x, y, z;
                xx = cone->GetPosition();
                x = xx[0];
                y = xx[1];
                z = xx[2];
                if (z<-Lscape) {
                    z = Lscape;
                    cone->SetPosition(x,y,z);
                }
            }
            //## Update frame counter
            frame += 1;
        }
    }
}

