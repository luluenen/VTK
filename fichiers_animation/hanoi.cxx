#include <vtkActor.h>
#include <vtkAssembly.h>
#include <vtkCamera.h>
#include <vtkConeSource.h>
#include <vtkContourFilter.h>
#include <vtkCubeSource.h>
#include <vtkCylinderSource.h>
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



//#include <String.h>
#include <time.h>
#include <unistd.h>
#include <ctime>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>



void set_rand_color(vtkSmartPointer<vtkActor> actor) {
    
    //## This is one possibility. Take a random triple
    //## and divide by the maximum
    double rgb[3];
    double maxx;
    srand(time(0));
    rgb[0] = rand();
    rgb[1] = rand();
    rgb[2] = rand();
    
    maxx = std::max(rgb[0],rgb[1]);
    maxx = std::max(maxx,rgb[2]);
    
    rgb[0] /= maxx;
    rgb[1] /= maxx;
    rgb[2] /= maxx;
    
    
    //## Set the color for the actor
    vtkSmartPointer<vtkProperty> prop ;
    prop = actor->GetProperty();
    prop->SetColor(rgb);
    
}


class disc
{
    public :
    int i,j;
    double pos, peg;
    int ndisc;
 //   double Hmove;
 //   double pegsep;
 //   double Hdisc;
    vtkSmartPointer<vtkActor> cyl_actor;
    
    disc(int j,double peg, double pos,vtkSmartPointer<vtkRenderer> ren, double Hdisc, double Rmin, double Rmax, int ndisc, double Hmove, double pegsep);
    void repos(double pegsep,double Hdisc, double Hmove);
    void repos(int pegto, int posto,int alpha,double pegsep,double Hdisc, double Hmove);
};


disc::disc(int j,double peg, double pos,vtkSmartPointer<vtkRenderer> ren, double Hdisc, double Rmin, double Rmax, int ndisc, double Hmove, double pegsep )
{
    double r;
    this->j = j;
    this->pos = pos;
    this->peg = peg;
  //  this->Hmove = Hmove;
  //  this->pegsep = pegsep;
  //  this->Hdisc = Hdisc;
    //## create te actor
    vtkSmartPointer<vtkCylinderSource> cylinder = vtkSmartPointer<vtkCylinderSource>::New();
    cylinder->SetHeight(Hdisc);
    r = Rmin+(Rmax-Rmin)/(ndisc-1)*j;
    cylinder->SetRadius(r);
    cylinder->SetResolution(20);
    
    vtkSmartPointer<vtkPolyDataMapper> cylMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cylMapper->SetInputConnection(cylinder->GetOutputPort() );
    cylMapper->ScalarVisibilityOff();
    
    cyl_actor = vtkSmartPointer<vtkActor>::New();
    cyl_actor->SetMapper(cylMapper);
    
    //## cyl_actor.SetPosition(0,self.pos*Hdisc,0.0)
    set_rand_color(cyl_actor);
    //  set_color(cyl_actor,j);
    repos(pegsep, Hdisc,  Hmove);
    ren->AddActor(cyl_actor);
}

//## Reposition the disc according to its rod and peg.
//## If (pegto,posto,alpha) are passed then it is
//## set to the position (1-alpha)*(rod,pos)+alpha*(rodto,posto)
//## and a correction for the height (int the `y' direction)
void disc::repos(int pegto, int posto,int alpha,double pegsep, double Hdisc, double Hmove) {
    double x[3];
    float xto[3];
    double a2;
    x[0] = (1-this->peg)*pegsep;
    x[1] = this->pos*Hdisc;
    x[2] =0;
    if (alpha != 0.0) {
        xto[0]= (1-pegto)*pegsep;
        xto[1]=posto*Hdisc;
        xto[2]=0;
        a2 = alpha+alpha*(1-alpha)*(2*alpha-1);
        for (int k=0;k<3;k++) {
            x[k] = x[k]*(1-a2)+xto[k]*a2;
        }
        x[1] += Hmove*alpha*(1-alpha);
    }
    cyl_actor->SetPosition(x);
}


void disc::repos(double pegsep,double Hdisc, double Hmove) {
    int pegto=-1;
    int posto=-1;
    int alpha=0.0;
    double x[3];
    float xto[3];
    double a2;
    x[0] = (1-peg)*pegsep;
    x[1] = pos*Hdisc;
    x[2] =0;
    if (alpha != 0.0) {
        xto[0]= (1-pegto)*pegsep;
        xto[1]=posto*Hdisc;
        xto[2]=0;
        a2 = alpha+alpha*(1-alpha)*(2*alpha-1);
        for (int k=0;k<3;k++) {
            x[k] = x[k]*(1-a2)+xto[k]*a2;
        }
        x[1] += Hmove*alpha*(1-alpha);
    }
    cyl_actor->SetPosition(x);
    
}

//## Set a random color for an actor

void set_color(vtkSmartPointer<vtkActor> actor, int j, double **colors){
    vtkSmartPointer<vtkProperty> prop ;
    prop = actor->GetProperty();
    prop->SetColor(colors[j]);
}

//## Add a peg actor
void addpeg(double x, double Hpeg, double Rpeg,vtkSmartPointer<vtkRenderer> ren){
    vtkSmartPointer<vtkCylinderSource> cylinder = vtkSmartPointer<vtkCylinderSource>::New();
    cylinder->SetHeight(Hpeg);
    cylinder->SetRadius(Rpeg);
    cylinder->SetResolution(20);
    
    vtkSmartPointer<vtkPolyDataMapper> cylMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cylMapper->SetInputConnection(cylinder->GetOutputPort() );
    cylMapper->ScalarVisibilityOff();
    
    vtkSmartPointer<vtkActor> cyl_actor = vtkSmartPointer<vtkActor>::New();
    cyl_actor->SetMapper(cylMapper);
    cyl_actor->SetPosition(x,Hpeg/2.0,0);
    vtkSmartPointer<vtkProperty> prop ;
    prop = cyl_actor->GetProperty();
    prop->SetColor(0.8,0.2,0.2);
    ren->AddActor(cyl_actor);
}

//## Prints one peg
void hprint1(char *peg,disc *discs[15]) {
    printf( "%s: ", peg);
    int taille;
    for (int d=0;d<15;d++){
        if (discs[d] != NULL) {
            printf("%d ",discs[d]->j);
        }
    }
}

//## prints all pegs
void hprint(disc *discs[3][15]) {
    hprint1("A",discs[0]);
    hprint1("B",discs[1]);
    hprint1("C",discs[2]);
    
}







void move1(disc *discs[3][15],int taillelist[3], int pegfrom, int pegto, int nframe, vtkSmartPointer<vtkRenderWindow>renWin, int render,double pegsep, double Hdisc, double Hmove) {
    int frame;
    //# take the top of `pegfrom'
    disc *d;
    double alpha;
    int posto;
    d = discs[pegfrom][taillelist[pegfrom]-1];
    discs[pegfrom][taillelist[pegfrom]] = NULL;
    taillelist[pegfrom] += -1;
    int tsleep = 1;
    //# Comput the height position in the new peg
    posto = taillelist[pegto];
    for (int k=0;k<nframe;k++) {
        //## Compute the weight for this frame
        alpha = float(k)/(nframe-1);
        //## Reposition the disc
        d->repos((double)(pegto),double(posto),alpha, pegsep, Hdisc, Hmove);
        if (render) {
            renWin->Render();
            //sleep(tsleep);
        }
    }
    
    //## Reset the internal position
    d->pos = posto;
    d->peg = pegto;
    //## Put the disc in the new peg
    discs[pegto][taillelist[pegto]]=d;
    taillelist[pegto] += 1;
    //# print "after move ",
    hprint(discs);
}

//## This is the `divide-and-conquer' algorithm for
//## solving the Tower of Hanoi problem.
//## Moves `nd' discs from peg `pegfrom' to `pegto'
void  move(disc *discs[3][15],int taillelist[3],int nd, int pegfrom,int pegto, int nframe, vtkSmartPointer<vtkRenderWindow>renWin, int render,double pegsep, double Hdisc, double Hmove) {
    int other;
    if (nd==1) {
        //## This ends the recursion
        move1(discs, taillelist, pegfrom, pegto, nframe, renWin, render, pegsep,  Hdisc,  Hmove);
    }
    else {
        //## Compute the remaining peg, uses
        //## the fact that the sum of the pegs gives 3
        other = 3-pegfrom-pegto;
        //## Move `nd-1' smaller discs from `pegfrom' to `other'
        move(discs,taillelist,nd-1,pegfrom,other, nframe, renWin, render, pegsep,  Hdisc,  Hmove);
        //## Move the larger disc to `pegto'
        move(discs,taillelist,1,pegfrom,pegto, nframe, renWin, render, pegsep,  Hdisc,  Hmove);
        //## Move back the `nd-1' discs from `other' to `pegto'
        move(discs,taillelist,nd-1,other,pegto, nframe, renWin, render, pegsep,  Hdisc,  Hmove);
    }
}

int main (int argc, char *argv[])
{
    double Rmax, Rmin, Rpeg, Hdiscs, Hpeg, pegsep, tsleep, Hdisc;
    int Hmove;
    int ndisc;
    int nframe, mkvideo, render;
    int taillelist[3];
    
    
    //## Show graphically the Towers of Hanoi problem
    
    Rmax = 1.0  ;    //                        # max radius of disc
    Rmin = 0.2   ;     //                      # min radius of disc
    Rpeg = 0.1    ;      //                    # radius of peg
    Hdiscs = 1.0   ;       ///                  # height of the total stack of discs
    Hpeg = 1.3*Hdiscs;        //               # Height of the peg
    Hmove = 5   ;               //             # height of the parabole that makes
    //# a disc during the move
    ndisc = 8  ;              //               # number of discs
    pegsep = 2.2*Rmax   ;       //             # separation of pegs
    tsleep = 0.003    ;           //           # sleep time while performing the
    //# visualization (but not video)
    nframe = 10  ;       //                    # Nbr of frames during a move
    mkvideo = 1 ;          //                  # make video or not
    render = 1;              //                # render or not (for debugging)
    
    //## Height a single disc
    Hdisc = std::min(Hdiscs/ndisc,0.3*Rmax);
    //# Create the Renderer, RenderWindow, and RenderWindowInteractor
    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren);
    renWin->SetSize(640, 480);
    
    //## This is an alternative to set_rand_color()
    //## First we generate a list of colors by taking
    //## first RGB and then interpolations between them , and so on
    /*                colors = np.zeros((3,3))
     colors[0] = np.array([1,0,0])
     colors[1] = np.array([0,1,0])
     colors[2] = np.array([0,0,1])
     while True:
     n = colors.shape[0]
     if n>=ndisc:
     break
     colors = np.resize(colors,[2*n,3])
     for k in range(n):
     k1 = (k+1)%n
     colors[n+k] = 0.5*(colors[k]+colors[k1])*/
    
    
    //## Add the 3 pegs ad x=0, and +-pegsep
    addpeg(-pegsep,  Hpeg,  Rpeg, ren);
    addpeg(0,  Hpeg,  Rpeg, ren);
    addpeg(+pegsep,  Hpeg,  Rpeg, ren);
    
    //## The discs are stored in a list of 3 lists
    disc *discs[3][15];
    for (int k=0;k<15; k++){
        discs[0][k]= NULL;
        discs[1][k]= NULL;
        discs[2][k]= NULL;
    }
    
    //## Initially all discs are in peg 0
    for (int k=0;k<ndisc; k++){
        //disc discobj
        discs[0][k]= new disc(ndisc-1-k,0,k, ren, Hdisc, Rmin, Rmax, ndisc, Hmove, pegsep);
    }
    
    //## color the background
    ren->SetBackground(0.7,0.7,0.7);
    
    //## prepare the camera
    vtkSmartPointer<vtkCamera> cam;
    cam = ren->GetActiveCamera();
    cam->SetViewUp(0.,1.,0.);
    cam->SetPosition(-3,3.5,7);
    cam->SetFocalPoint(0,0,0);
    
    //## This is used to store the frames
    //## for creating a movie
    
    vtkSmartPointer<vtkWindowToImageFilter> w2i = vtkSmartPointer<vtkWindowToImageFilter>::New();
    w2i->SetInput(renWin);
    w2i->Update();
    
    //## The TIFF writer
    vtkSmartPointer<vtkTIFFWriter> writer = vtkSmartPointer<vtkTIFFWriter>::New();
    writer->SetInputConnection(w2i->GetOutputPort());
    writer->SetCompressionToJPEG();
    
    if (render) {
        renWin->Render();
        hprint(discs);
    }
    
    //## `frame' tracks the number of frame for generating
    //## the frame image;  filename
    int frame = 0;
    
    //## Moves the top disc from peg `pegfrom' to `pegto'
    
    //## This is to give the user time
    //## to clean the VTK window
    //raw_input("press <enter>");
    int flagrender = 1;
    taillelist[0] = 8;
    taillelist[1] = 0;
    taillelist[2] = 0;
    move(discs,taillelist, ndisc,0,1, nframe, renWin, flagrender, pegsep,  Hdisc,  Hmove);
    //raw_input("press <enter>");
}

/* #!/usr/bin/env python

import time
import random
import sys
import os
# from np import *
import numpy as np
import vtk

## Show graphically the Towers of Hanoi problem

Rmax = 1.0                              # max radius of disc
Rmin = 0.2                              # min radius of disc
Rpeg = 0.1                              # radius of peg
Hdiscs = 1.0                            # height of the total stack of discs
Hpeg = 1.3*Hdiscs                       # Height of the peg
Hmove = 5                               # height of the parabole that makes
                                        # a disc during the move
ndisc = 8                               # number of discs
pegsep = 2.2*Rmax                       # separation of pegs
tsleep = 0.003                          # sleep time while performing the
                                        # visualization (but not video)
nframe = 10                             # Nbr of frames during a move
mkvideo = 1                             # make video or not
render = 1                              # render or not (for debugging)

## Height a single disc
Hdisc = min(Hdiscs/ndisc,0.3*Rmax)
# Create the Renderer, RenderWindow, and RenderWindowInteractor
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(640, 480)

## Set a random color for an actor
def set_rand_color(actor):
    if 0:
        ## This is one possibility. Take a random triple
        ## and divide by the maximum
        r = random.random()
        g = random.random()
        b = random.random()
        maxx = max(r,g,b)
        r /= maxx
        g /= maxx
        b /= maxx
    else:
        ## Colors are of the form
        ## x,1-x,0 and permutations
        color = np.random.rand(3)
        k = random.randint(0,2)
        color[k] = 0
        color /= color.sum()
        print color
    ## Set the color for the actor
    prop = actor.GetProperty()
    prop.SetColor(color)

## This is an alternative to set_rand_color()
## First we generate a list of colors by taking
## first RGB and then interpolations between them , and so on
colors = np.zeros((3,3))
colors[0] = np.array([1,0,0])
colors[1] = np.array([0,1,0])
colors[2] = np.array([0,0,1])
while True:
    n = colors.shape[0]
    if n>=ndisc:
        break
    colors = np.resize(colors,[2*n,3])
    for k in range(n):
        k1 = (k+1)%n
        colors[n+k] = 0.5*(colors[k]+colors[k1])

## After, the `set_color()' for `j' takes the `j'-th
## color from the generated list
def set_color(actor,j):
    prop = actor.GetProperty()
    prop.SetColor(colors[j])

## The disc contains its index (in [0,ndisc)) and
## its position (peg and position)
class disc:
    def __init__(self,j,peg,pos):
        self.j = j
        self.pos = pos
        self.peg = peg
        ## create te actor
        cylinder = vtk.vtkCylinderSource()
        cylinder.SetHeight(Hdisc)
        r = Rmin+(Rmax-Rmin)/(ndisc-1)*j
        cylinder.SetRadius(r)
        cylinder.SetResolution(20)
        
        cylMapper = vtk.vtkPolyDataMapper()
        cylMapper.SetInputConnection(cylinder.GetOutputPort() )
        cylMapper.ScalarVisibilityOff()
        
        cyl_actor = vtk.vtkActor()
        cyl_actor.SetMapper(cylMapper)
        self.actor = cyl_actor
        ## cyl_actor.SetPosition(0,self.pos*Hdisc,0.0)
        ## set_rand_color(cyl_actor)
        set_color(cyl_actor,j)
        self.repos()
        ren.AddActor(cyl_actor)

    ## Reposition the disc according to its rod and peg.
    ## If (pegto,posto,alpha) are passed then it is
    ## set to the position (1-alpha)*(rod,pos)+alpha*(rodto,posto)
    ## and a correction for the height (int the `y' direction)
    def repos(self,pegto=-1,posto=-1,alpha=0.0):
        x = np.array([(1-self.peg)*pegsep,self.pos*Hdisc,0])
        if alpha!=0.0:
            xto = np.array([(1-pegto)*pegsep,posto*Hdisc,0])
            a2 = alpha+alpha*(1-alpha)*(2*alpha-1)
            x = x*(1-a2)+xto*a2
            x[1] += Hmove*alpha*(1-alpha)
        self.actor.SetPosition(x)

## Add a peg actor
def addpeg(x):
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetHeight(Hpeg)
    cylinder.SetRadius(Rpeg)
    cylinder.SetResolution(20)

    cylMapper = vtk.vtkPolyDataMapper()
    cylMapper.SetInputConnection(cylinder.GetOutputPort() )
    cylMapper.ScalarVisibilityOff()

    cyl_actor = vtk.vtkActor()
    cyl_actor.SetMapper(cylMapper)
    cyl_actor.SetPosition(x,Hpeg/2.0,0)
    prop = cyl_actor.GetProperty()
    prop.SetColor(0.8,0.2,0.2)
    ren.AddActor(cyl_actor)

## Add the 3 pegs ad x=0, and +-pegsep
addpeg(-pegsep)
addpeg(0)
addpeg(+pegsep)

## The discs are stored in a list of 3 lists
discs = [[],[],[]]
## Initially all discs are in peg 0
for k in range(ndisc):
    discs[0].append(disc(ndisc-1-k,0,k))

## Prints one peg
def hprint1(peg,discs):
    print "%s: " % peg,
    for d in discs:
        print "%d " % d.j,

## prints all pegs
def hprint(discs):
    hprint1("A",discs[0])
    hprint1("B",discs[1])
    hprint1("C",discs[2])
    print

## color the background
ren.SetBackground(0.7,0.7,0.7)

## prepare the camera
cam = ren.GetActiveCamera()
cam.SetViewUp(0.,1.,0.);
cam.SetPosition(-3,3.5,7)
cam.SetFocalPoint(0,0,0)

## This is used to store the frames
## for creating a movie
w2i = vtk.vtkWindowToImageFilter()
w2i.SetInput(renWin)
w2i.Update()

## The TIFF writer
writer = vtk.vtkTIFFWriter()
writer.SetInputConnection(w2i.GetOutputPort())
writer.SetCompressionToJPEG()

if render:
    renWin.Render()
hprint(discs)

## `frame' tracks the number of frame for generating
## the frame image filename
frame = 0

## Moves the top disc from peg `pegfrom' to `pegto'
def move1(discs,pegfrom,pegto):
    global frame
    # take the top of `pegfrom'
    d = discs[pegfrom].pop()
    # Comput the height position in the new peg
    posto = len(discs[pegto])
    for k in range(nframe):
        ## Compute the weight for this frame
        alpha = float(k)/(nframe-1)
        ## Reposition the disc
        d.repos(pegto,posto,alpha)
        if render:
            renWin.Render()
            if mkvideo:
                w2i.Modified()
                ## generate the fram filename
                assert os.path.isdir("./YUV")
                tiff = "./YUV/frame.%d.tiff" % frame
                yuv = "./YUV/frame.%d.yuv" % frame
                writer.SetFileName(tiff)
                writer.Write()
                ## Convert to YUV and gzip
                os.system("convert %s %s ; gzip -f %s" % (tiff,yuv,yuv))
                os.unlink(tiff)
                frame += 1
            else:
                time.sleep(tsleep)
    ## Reset the internal position
    d.pos = posto
    d.peg = pegto
    ## Put the disc in the new peg
    discs[pegto].append(d);
    # print "after move ",
    hprint(discs)

## This is the `divide-and-conquer' algorithm for
## solving the Tower of Hanoi problem.
## Moves `nd' discs from peg `pegfrom' to `pegto'
def move(discs,nd,pegfrom,pegto):
    if nd==1:
        ## This ends the recursion
        move1(discs,pegfrom,pegto)
    else:
        ## Compute the remaining peg, uses
        ## the fact that the sum of the pegs gives 3
        other = 3-pegfrom-pegto
        ## Move `nd-1' smaller discs from `pegfrom' to `other'
        move(discs,nd-1,pegfrom,other)
        ## Move the larger disc to `pegto'
        move(discs,1,pegfrom,pegto)
        ## Move back the `nd-1' discs from `other' to `pegto'
        move(discs,nd-1,other,pegto)

## This is to give the user time
## to clean the VTK window 
raw_input("press <enter>")
move(discs,ndisc,0,1)
raw_input("press <enter>")*/
