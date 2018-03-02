/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: Cylinder.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This simple example shows how to do basic rendering and pipeline
// creation using C++.
//
#include "vtkCylinderSource.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"

#include "showpath.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

double thickness = 1;

void makeScene(void)
{
	double fibreColor[] = {1.0, 0.5, 0.3};
	double cellColor[] = {1.0, 0.3, 0.0};
//	double cellSize = 0.05;
	int i, j, k;
	CELL_POS *cp1, *cp2;
	double fpos[3], v[3];
	double minpos[3], scalepos[3], tickpos[3];
	double Pi = 3.15159;
	double r, g, b;
	double dr=0.15;
	double dg=0.10;
	double db=0.15;
	double dtick;
	double ticklen=10;
	vtkActor *cylinderActor;
	vtkActor *cellActor;
	
	ren1->SetBackground(1.0, 1.0, 1.0);
	minpos[0] = minpos[1] = minpos[2] = 1.0e10;
	for (i=1; i<=npaths; i++) {
		PATH_POS *p = &path[i];
		for (k=0; k<p->nlinks; k++) {
			cp1 = &p->cell[k];
			for (j=0;j<3;j++) {
				minpos[j] = MIN(minpos[j],cp1->pos[j]);
			}
		}
	}
	printf("minpos: %f %f %f\n",minpos[0],minpos[1],minpos[2]);
	scalepos[0] = 0;
	scalepos[1] = 1.5*minpos[1];
	scalepos[2] = 0;	//1.5*minpos[2];
	// Make a scale bar
	r = 0.0;
	g = 0.0;
	b = 0.0;
	fibreColor[0] = r;
	fibreColor[1] = g;
	fibreColor[2] = b;
	/*
	cylinderActor = vtkActor::New();
	cylinderActor->SetMapper(cylinderMapper);
	cylinderActor->GetProperty()->SetColor(fibreColor);
	double s[] = {1.5, 50.0, 1.5};
	cylinderActor->SetScale(s);
	cylinderActor->SetPosition(scalepos);
	cylinderActor->RotateWXYZ(90.0,0,0,1);
	ren1->AddActor(cylinderActor);
	printf("made scale bar\n");
*/
	// Make axes
	double axisshiftfactor = 1.0;
	minpos[0] = axisshiftfactor*minpos[0];
	minpos[1] = axisshiftfactor*minpos[1];
	minpos[2] = axisshiftfactor*minpos[2];
	double axislen;	// = 300;
    cout << "Enter the axis length (multiple of 100) : "; 
    cin >> axislen; 
	int nticks = axislen/100;
	for (k=0;k<3;k++) {
		cylinderActor = vtkActor::New();
		cylinderActor->SetMapper(cylinderMapper);
		cylinderActor->GetProperty()->SetColor(fibreColor);
		double s[] = {1.5, axislen, 1.5};
		cylinderActor->SetScale(s);
		if (k == 0) {
			minpos[0] += axislen/2;
			cylinderActor->SetPosition(minpos);
			cylinderActor->RotateWXYZ(90.0,0,0,1);	// x axis
			ren1->AddActor(cylinderActor);
			minpos[0] -= axislen/2;
		} else if (k == 1) {
			minpos[1] += axislen/2;
			cylinderActor->SetPosition(minpos);
			cylinderActor->RotateWXYZ(90.0,0,1,0);	// y axis
			ren1->AddActor(cylinderActor);
			minpos[1] -= axislen/2;
		} else if (k == 2) {
			minpos[2] += axislen/2;
			cylinderActor->SetPosition(minpos);
			cylinderActor->RotateWXYZ(90.0,1,0,0);	// z axis
			ren1->AddActor(cylinderActor);
			minpos[2] -= axislen/2;
		}
		// Make ticks
		dtick = axislen/nticks - 0.5;
		for (int ktick=0; ktick<nticks; ktick++) {
			cylinderActor = vtkActor::New();
			cylinderActor->SetMapper(cylinderMapper);
			cylinderActor->GetProperty()->SetColor(fibreColor);
			double s[] = {1.5, ticklen, 1.5};
			cylinderActor->SetScale(s);
			if (k == 0) {	// x axis
				tickpos[0] = minpos[0] + (ktick+1)*dtick;
				tickpos[1] = minpos[1] + ticklen/2;
				tickpos[2] = minpos[2] +  0;
				cylinderActor->SetPosition(tickpos);
				cylinderActor->RotateWXYZ(90.0,0,1,0);	// y axis
				ren1->AddActor(cylinderActor);
			} else if (k == 1) {	// y axis
				tickpos[0] = minpos[0] + ticklen/2;
				tickpos[1] = minpos[1] + (ktick+1)*dtick;
				tickpos[2] = minpos[2] +  0;
				cylinderActor->SetPosition(tickpos);
				cylinderActor->RotateWXYZ(90.0,0,0,1);	// x axis
				ren1->AddActor(cylinderActor);
			} else if (k == 2) {	// z axis
				tickpos[0] = minpos[0] +  0;
				tickpos[1] = minpos[1] + ticklen/2;
				tickpos[2] = minpos[2] + (ktick+1)*dtick;
				cylinderActor->SetPosition(tickpos);
				cylinderActor->RotateWXYZ(90.0,0,1,0);	// y axis
				ren1->AddActor(cylinderActor);
			}
		}
	}

	r = 1.0;
	g = 0.3;
	b = 0.2;
	for (i=1; i<=npaths; i++) {
		if (i%3 == 0) r -= dr;
		if (r < 0) r = 1;
		if (i%3 == 1) g += dg;
		if (g > 1) g = 0;
		if (i%3 == 2) b += db;
		if (b > 1) b = 0;
		fibreColor[0] = r;
		fibreColor[1] = g;
		fibreColor[2] = b;
		PATH_POS *p = &path[i];
		printf("path: %d nlinks: %d\n",i,p->nlinks);
		for (k=0; k<p->nlinks; k++) {
			cylinderActor = vtkActor::New();
			cylinderActor->SetMapper(cylinderMapper);
			cylinderActor->GetProperty()->SetColor(fibreColor);
			cp1 = &p->cell[k];
			cp2 = &p->cell[k+1];
//			printf("link: %d cp1: %f %f %f cp2: %f %f %f\n",k,cp1->pos[0],cp1->pos[1],cp1->pos[2],cp2->pos[0],cp2->pos[1],cp2->pos[2]);
			for (j=0; j<3; j++) {
				fpos[j] = (cp1->pos[j] + cp2->pos[j])/2;
				v[j] = cp1->pos[j] - cp2->pos[j];
			}
			double v_mod = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
			double s[] = {0.5, v_mod, 0.5};
			cylinderActor->SetScale(s);
			for (j=0; j<3; j++)
				v[j] = v[j]/v_mod;

			double sina = sqrt(v[0]*v[0] + v[2]*v[2]);
			double cosa = v[1];
			double theta = asin(sina)*(180.0/Pi);
			if (cosa < 0)
				theta = 180 - theta;

			cylinderActor->SetPosition(fpos);
			cylinderActor->RotateWXYZ(theta,v[2],0,-v[0]);

			ren1->AddActor(cylinderActor);
			cellActor = vtkActor::New();
			cellActor->SetMapper(cellMapper);
			cellActor->GetProperty()->SetColor(fibreColor);	//(cellColor);
			cellActor->SetPosition(cp1->pos);
			ren1->AddActor(cellActor);			
		}
	}
// Mark start point
	cellActor = vtkActor::New();
	cellActor->SetMapper(cellMapper);
	cellActor->GetProperty()->SetColor(1.0,1.0,0.0);	//(cellColor);
	cellActor->SetPosition(0,0,0);
	cellActor->SetScale(5);
	ren1->AddActor(cellActor);			


//  cylinderActor->GetProperty()->SetColor(1.0,0.0,0.0);
	/*
	for (i=1; i<=ncells; i++) {
		CELL_POS *p = &cell[i];

		cellActor = vtkActor::New();
		cellActor->SetMapper(cellMapper);
		cellActor->GetProperty()->SetColor(cellColor);
		cellActor->SetPosition(p->pos);
		ren1->AddActor(cellActor);
//		printf("Cell pos: %d %f %f %f\n",i,p->pos[0],p->pos[1],p->pos[2]);
	}
	*/
//  cellActor->SetScale(cellSize);

}

// I have forgotten how to read/parse a line with integer and float.
// For now using integer (x,y,z), no problem.
void readGeometry(char *pathfile)
{
	char line[128];

	fp = fopen(pathfile,"r");
	fgets(line, 128, fp);
	sscanf (line, "%d", &npaths);
	printf("npaths: %d\n",npaths);
	path = (PATH_POS*)malloc((npaths+1)*sizeof(PATH_POS));		// 1-based 
	for (int i=1; i<=npaths; i++) {
		PATH_POS *p = &path[i];
		fgets(line, 128, fp);
		sscanf (line, "%d", &p->nlinks);
//		printf("path: %d nlinks: %d\n",i,p->nlinks);
		p->cell = (CELL_POS *)malloc((p->nlinks+1)*sizeof(CELL_POS));
		for (int k=1; k<=p->nlinks; k++) {
			fgets(line, 128, fp);
			sscanf (line, "%lf %lf %lf", &p->cell[k].pos[0],&p->cell[k].pos[1],&p->cell[k].pos[2]);
//			printf("link: %d cell: %f %f %f\n",k,p->cell[k].pos[0],p->cell[k].pos[1],p->cell[k].pos[2]);
		}
	}
	fclose(fp);
}

int main(int argc, char**argv)
{
	char *pathfile;
	if (argc != 2) {
		printf("Usage: showpath pathfile\n");
		return 1;
	}
	pathfile = argv[1];
  	vtkCylinderSource *cylinder = vtkCylinderSource::New();
  	cylinder->SetResolution(12);
	cylinder->SetRadius(thickness);
    cylinder->SetHeight(1);

  // The mapper is responsible for pushing the geometry into the graphics
  // library. It may also do color mapping, if scalars or other attributes
  // are defined.
  //

//  vtkPolyDataMapper *cylinderMapper = vtkPolyDataMapper::New();
  cylinderMapper = vtkPolyDataMapper::New();
  cylinderMapper->SetInputConnection(cylinder->GetOutputPort());
	
	vtkSphereSource *cell = vtkSphereSource::New();
	cell->SetThetaResolution(12);
	cell->SetPhiResolution(12);
	cell->SetRadius(thickness/2);
	cellMapper = vtkPolyDataMapper::New();
	cellMapper->SetInputConnection(cell->GetOutputPort());
	
  // Create the graphics structure. The renderer renders into the
  // render window. The render window interactor captures mouse events
  // and will perform appropriate camera or actor manipulation
  // depending on the nature of the events.
  //

//  vtkRenderer *ren1 = vtkRenderer::New();
  ren1 = vtkRenderer::New();
//  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren1);
//  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);
  iren->LightFollowCameraOn();


	readGeometry(pathfile);
	printf("Did readgeometry\n");
  // Add the actors to the renderer, set the background and size
  makeScene();
  renWin->SetSize(900, 900);
  // We'll zoom in a little by accessing the camera and invoking a "Zoom"
  // method on it.
  ren1->ResetCamera();
  //ren1->GetActiveCamera()->Zoom(1.5);
//  ren1->GetActiveCamera()->SetParallelProjection(1);
  renWin->Render();

  // This starts the event loop and as a side effect causes an initial render.
  iren->Start();

  // Exiting from here, we have to delete all the instances that
  // have been created.
  cylinder->Delete();
  cylinderMapper->Delete();
 // cylinderActor->Delete();
  cell->Delete();
  cellMapper->Delete();
  ren1->Delete();
  renWin->Delete();
  iren->Delete();

  return 0;
}




