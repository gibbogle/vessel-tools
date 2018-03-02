/* showFRC.h */

#include <stdio.h>

vtkRenderer *ren1;
vtkRenderWindow *renWin;
vtkRenderWindowInteractor *iren;
vtkPolyDataMapper *cylinderMapper;
vtkPolyDataMapper *cellMapper;

FILE* fp;

struct cell_pos {
	double pos[3];
};
typedef cell_pos CELL_POS;

//struct fibre_pos {
//	int tag;			// 1-based
//	int cell1, cell2;	// 1-based
//	double length;
//};
//typedef fibre_pos FIBRE_POS;

struct path_pos {
	int nlinks;
	CELL_POS *cell;
};
typedef path_pos PATH_POS;

int npaths;
PATH_POS *path;

//FIBRE_POS *fibre;
//CELL_POS *cell;
//int ncells, nfibres;
