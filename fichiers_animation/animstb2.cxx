#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVolume16Reader.h>
#include <vtkDICOMImageReader.h>
#include <vtkActor.h>
#include <vtkOutlineFilter.h>
#include <vtkCamera.h>
#include <vtkCubeSource.h>
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>

#include <vtkPolyDataNormals.h>
#include <vtkContourFilter.h>
#include <vtkMarchingCubes.h>
#include <vtkSmartPointer.h>
#include <vtkWindowToImageFilter.h>
#include <vtkTIFFWriter.h>
#include <vtkVectorText.h>
//#include <Math.h>
//#include <String.h>
#include <time.h>
#include <unistd.h>
#include <ctime>
void  replot(int i,double z, vtkActor **cub_actors, vtkActor **txt_actors, int L, int Lz, double *xcubes) {
    cub_actors[i]->SetPosition(xcubes[i],0.0,z);
    txt_actors[i]->SetPosition(xcubes[i]-L/2,-0.25*L,z+Lz/2+0.1);
}


void swap1(vtkActor **v,int i, int j){
    vtkActor *tmp = v[i];
    v[i] = v[j];
    v[j] = tmp;
}

void swap2(double *v,int i, int j){
    double tmp = v[i];
    v[i] = v[j];
    v[j] = tmp;
}

void swap3(int *v,int i, int j){
    int tmp = v[i];
    v[i] = v[j];
    v[j] = tmp;
}

void exchange(int i,int j, int exc, double *xcubes,vtkActor **cub_actors, vtkActor **txt_actors, int L, int Lz, vtkSmartPointer<vtkRenderWindow> renWin, double *vals) {
    double fac;
    int nfx = 1;
    double xcubei0 = xcubes[i];
    double xcubej0 = xcubes[j];
    for (int k=0;k<nfx;k++) {
        fac = exc*float(k+1)/float(nfx);
        xcubes[i] = xcubei0 + fac*(xcubej0-xcubei0);
        xcubes[j] = xcubej0 - fac*(xcubej0-xcubei0);
        replot(i, 0.1, cub_actors, txt_actors,  L,  Lz, xcubes);
        replot(j, 0.2, cub_actors, txt_actors, L, Lz, xcubes);
        renWin->Render();
        usleep(200000);
    }
    replot(i,0, cub_actors, txt_actors,  L,  Lz, xcubes);
    replot(j,0, cub_actors, txt_actors,  L,  Lz, xcubes);
    renWin->Render();
    usleep(200000);
    if (exc) {
        swap1(cub_actors,i,j);
    }
    swap1(txt_actors,i,j);
    swap2(xcubes,i,j);
    swap2(vals,i,j);
}

void bubble(int nnum,double *vals,double *xcubes,vtkActor **cub_actors, vtkActor **txt_actors, int L, int Lz, vtkSmartPointer<vtkRenderWindow> renWin) {
    // # Bubble sort
    for (int k =0;k<nnum;k++){
        for (int l=nnum-1; l>k;l--){
            if (vals[l]<vals[l-1]) {
                exchange(l,l-1,1, xcubes, cub_actors, txt_actors,  L,  Lz, renWin, vals);
            }
        }
    }
}

int generer_bornes(int, int);
void initialiser_aleat(unsigned int);

int appele_srand = 0;

int generer_bornes(int min, int max)
{
    if(appele_srand != 1)
        initialiser_aleat((unsigned)time(NULL));
    return rand()%(max-min+1) + min;
}

void initialiser_aleat(unsigned int n)
{
    srand(n);
    appele_srand = 1;
}

int partition(int j1,int j2, double *vals, double *xcubes,vtkActor **cub_actors, vtkActor **txt_actors, int L, int Lz, vtkSmartPointer<vtkRenderWindow> renWin)
{
    int j, l, r;
    int n ;
    double pivot;
    
    n = j2-j1;
    
    //print "j1 %d, j2 %d, n %d" % (j1,j2,n)
    initialiser_aleat((unsigned)time(NULL));
    j = generer_bornes(j1,j2-1);
    
    pivot = vals[j];
    
    //   ## verify there is at least one element lower
    int ok = 0;
    int k;
    k=j1;
    int flag =0;
    while (flag ==0 && k<=j2) {
        if (vals[k]<pivot){
            ok = 1;
            flag = 1;
        }
        else if (vals[k]>pivot) {
            ok = 1;
            pivot = vals[k];
            flag = 1;
        }
        k++;
    }
    
    if (!ok) {
        //print "(%d,%d) done" % (j1,j2)
        return j2;
    }
    //print vals
    l=j1;
    r=j2-1;
    flag = 0;
    while (flag == 0) {
        while (vals[l]<pivot) {
            l += 1;
        }
        while (vals[r]>=pivot) {
            r -= 1;
        }
        //        # print "vals[l=%d] %d, vals[r=%d] = %d" % (l,vals[l],r,vals[r])
        if (l>r) {
            flag = 1;
        }
        else {
            exchange(r,l,1,xcubes,cub_actors, txt_actors,  L,  Lz,  renWin, vals);
        }
    }
    
    return l;
}


void  quick(int j1, int j2, double *vals, double *xcubes,vtkActor **cub_actors, vtkActor **txt_actors, int L, int Lz, vtkSmartPointer<vtkRenderWindow> renWin) {
    // # print "entering quick(%d,%d)" % (j1,j2)
    int l;
    if (j2-j1<=1) {
        
    }
    else {
        l = partition(j1,j2,vals, xcubes,cub_actors, txt_actors,  L,  Lz,  renWin);
        if (l==j2) {
            
        }
        else {
            quick(j1, l, vals, xcubes, cub_actors, txt_actors,  L,  Lz,  renWin);
            quick(l, j2, vals, xcubes, cub_actors, txt_actors,  L,  Lz,  renWin);
        }
    }
}


int main (int argc, char *argv[])
{
    
    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren);
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    
    renWin->SetSize(640, 480);
    
    ren->SetBackground(0.7,0.7,0.7);
    double L, Lz, Lx;
    int nnum;
    L = 1;
    Lz = 0.5;
    nnum = 50;
    Lx = 10.0/nnum;
    
    //# Create the Renderer, RenderWindow, and RenderWindowInteractor
    
    
    vtkActor *cub_actors[nnum] ;
    vtkActor *txt_actors[nnum] ;
    double xcubes[nnum] ;
    double xcube = -Lx*nnum/2;
    double vals[nnum] ;
    int maxval = nnum;
    srand((unsigned)time(0));
    int val;
    for (int k=0 ;k<nnum;k++) {
        
        val = rand()%maxval;
        vals[k]= val;
        
        vtkSmartPointer<vtkCubeSource> cube = vtkSmartPointer<vtkCubeSource>::New();
        cube->SetXLength(Lx);
        cube->SetYLength(L);
        cube->SetZLength(Lz);
        
        vtkSmartPointer<vtkPolyDataMapper> cubMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        cubMapper->SetInputConnection(cube->GetOutputPort() );
        cubMapper->ScalarVisibilityOff();
        
        cub_actors[k] = vtkActor::New();
        cub_actors[k]->SetMapper(cubMapper);
        cub_actors[k]->SetPosition(xcube,0.0,0.0);
        
        double alpha = double(val)/double(maxval);
        vtkSmartPointer<vtkProperty> prop = vtkSmartPointer<vtkProperty>::New();
        prop = cub_actors[k]->GetProperty();
        if (alpha<0.5) {
            alpha *= 2.0;
            prop->SetColor(0,alpha,1-alpha);
        }
        else
        {
            alpha = 2*(alpha-0.5);
            prop->SetColor(alpha,1-alpha,0.0);
        }
        prop->SetOpacity(1.0);
        prop->SetSpecular(1.0);
        prop->SetDiffuse(1.);
        prop->EdgeVisibilityOn();
        
        vtkSmartPointer<vtkVectorText> text = vtkSmartPointer<vtkVectorText>::New();
        char s[100];
        sprintf(s,"%d",val);
        text->SetText(s);
        
        vtkSmartPointer<vtkPolyDataMapper> txtMapper  = vtkSmartPointer<vtkPolyDataMapper>::New();
        txtMapper->SetInputConnection(text->GetOutputPort() );
        
        txt_actors[k] = vtkActor::New();
        txt_actors[k]->SetMapper(txtMapper);
        txt_actors[k]->SetPosition(xcube-L/2,-0.25*L,Lz/2+0.1);
        
        prop = txt_actors[k]->GetProperty();
        prop->SetColor(0,0,0);
        
        //# Add the actors to the renderer
        ren->AddActor(cub_actors[k]);
        //ren->AddActor(txt_actors[k]);
        
        xcubes[k]=xcube;
        xcube += Lx;
    }
    
    
    vtkSmartPointer<vtkCamera> cam = vtkSmartPointer<vtkCamera>::New();
    cam = ren->GetActiveCamera();
    cam->SetViewUp(0.,1.,0.);
    cam->SetFocalPoint(-1,0,0);
    
    //## This is used to store the frames
    //## for creating a movie
    vtkSmartPointer<vtkWindowToImageFilter>w2i = vtkSmartPointer<vtkWindowToImageFilter>::New();
    w2i->SetInput(renWin);
    //w2i->Update();
    
    //## The TIFF writer
    vtkSmartPointer<vtkTIFFWriter> writer = vtkSmartPointer<vtkTIFFWriter>::New();
    writer->SetInputConnection(w2i->GetOutputPort());
    writer->SetCompressionToJPEG();
    writer->SetFileName("toto.jpg");

    
    double xcam[3] = {-5*1.5,3*1.5,8*1.5};
    
    cam->SetPosition(xcam);
    // Set a background color for the renderer and set the size of the
    // render window (expressed in pixels).
    ren->SetBackground(.2, .3, .4);
    renWin->SetSize(1000, 800);
    
    // Note that when camera movement occurs (as it does in the Dolly()
    // method), the clipping planes often need adjusting. Clipping planes
    // consist of two planes: near and far along the view direction. The
    // near plane clips out objects in front of the plane; the far plane
    // clips out objects behind the plane. This way only what is drawn
    // between the planes is actually rendered.
    ren->ResetCameraClippingRange ();
     renWin->Render();
  //  writer->Write();
    // Initialize the event loop and then start it.
    bubble(nnum, vals, xcubes, cub_actors, txt_actors,  L,  Lz,  renWin);
    //quick(0,nnum, vals, xcubes, cub_actors, txt_actors,  L,  Lz,  renWin);
    iren->Initialize();
    iren->Start();
    
    return EXIT_SUCCESS;
}

