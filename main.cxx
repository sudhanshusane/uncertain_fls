#include <iostream>
#include <math.h>
#include <vtkSmartPointer.h>
#include <vtkRectilinearGridReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkAbstractArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkDataSetWriter.h>


using namespace std;

float compute_distance(float *p, float *q)
{
	return sqrt((pow((p[0]-q[0]),2.0) + pow((p[1]-q[1]),2.0) + pow((p[2]-q[2]),2.0)));
}

int GetNumberOfPoints(const int *dims)
{
	return dims[0]*dims[1]*dims[2];
}

int GetNumberOfCells(const int *dims)
{
	return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
}

int GetPointIndex(const int *idx, const int *dims)
{
	return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
}

int GetCellIndex(const int *idx, const int *dims)
{
	return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
}

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
	idx[0] = pointId%dims[0];
	idx[1] = (pointId/dims[0])%dims[1];
	idx[2] = pointId/(dims[0]*dims[1]);
}

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
	idx[0] = cellId%(dims[0]-1);
	idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
	idx[2] = cellId/((dims[0]-1)*(dims[1]-1));
}

int main(int argc, char* argv[])
{
  
	vtkSmartPointer<vtkRectilinearGridReader> reader = 
    vtkSmartPointer<vtkRectilinearGridReader>::New();

	reader->SetFileName(argv[1]);
	reader->Update();

	vtkSmartPointer<vtkRectilinearGrid> mesh = 
    vtkSmartPointer<vtkRectilinearGrid>::New();
	mesh = reader->GetOutput();
	int num_pts = mesh->GetNumberOfPoints();

	vtkAbstractArray* a1 = mesh->GetPointData()->GetArray(argv[7]);
	vtkAbstractArray* a2 = mesh->GetPointData()->GetArray(argv[8]);
//	vtkAbstractArray* a3 = mesh->GetPointData()->GetArray("streamvort");
//	vtkAbstractArray* a4 = mesh->GetPointData()->GetArray("vortmag");

	vtkFloatArray* att1 = vtkFloatArray::SafeDownCast(a1);
	vtkFloatArray* att2 = vtkFloatArray::SafeDownCast(a2);
//	vtkFloatArray* att3 = vtkFloatArray::SafeDownCast(a3);
//	vtkFloatArray* att4 = vtkFloatArray::SafeDownCast(a4);

	int dims[3];

  mesh->GetDimensions(dims);

	cout << num_pts << endl;

// Output the range of each attribute 

	double range[2]; 
	att1->GetRange(range);
	cout << "Att1 range " << range[0] << " to " << range[1] << endl;
	
	att2->GetRange(range);

	cout << "Att2 range " << range[0] << " to " << range[1] << endl;

	float att1_min = atof(argv[3]);
	float att1_max = atof(argv[4]);
	
	float att2_min = atof(argv[5]);
	float att2_max = atof(argv[6]);

	cout << "Interval start: " << att1_min << endl;
	cout << "Interval end: " << att1_max << endl;
	cout << "Interval start: " << att2_min << endl;
	cout << "Interval end: " << att2_max << endl;

	/*
 * Get the coordinates into three arrays - done
 * Write a function to get the x,y,z, indices based on index - done
 * Write a function to get the index based on x,y,z - done
 * Set up a basic way to specify a trait for each dimension. - done
 * Loop over grid points and set 1 or 0 depending on analysis of each dimension. - done
 * Write to a vtk file - done
 * Implement EDT algorithm
 * Write distance field to a vtk file
 */

	float *xc = (float*) mesh->GetXCoordinates()->GetVoidPointer(0);
	float *yc = (float*) mesh->GetYCoordinates()->GetVoidPointer(0);
	float *zc = (float*) mesh->GetZCoordinates()->GetVoidPointer(0);

	int *binary_image = (int*)malloc(sizeof(float)*num_pts);
	float* stddev = (float*)malloc(sizeof(float)*num_pts);

	vtkSmartPointer<vtkFloatArray> binaryArray = vtkSmartPointer<vtkFloatArray>::New();
	binaryArray->SetName("binary");
	binaryArray->SetNumberOfComponents(1); 
	
	vtkSmartPointer<vtkFloatArray> distArray = vtkSmartPointer<vtkFloatArray>::New();
	distArray->SetName("distance");
	distArray->SetNumberOfComponents(1); 
	
	vtkSmartPointer<vtkFloatArray> distArray_99 = vtkSmartPointer<vtkFloatArray>::New();
	distArray->SetName("distance99");
	distArray->SetNumberOfComponents(1); 
	vtkSmartPointer<vtkFloatArray> distArray_95 = vtkSmartPointer<vtkFloatArray>::New();
	distArray->SetName("distance95");
	distArray->SetNumberOfComponents(1); 
	vtkSmartPointer<vtkFloatArray> distArray_68 = vtkSmartPointer<vtkFloatArray>::New();
	distArray->SetName("distance68");
	distArray->SetNumberOfComponents(1); 

	for(int n = 0; n < num_pts; n++)
	{
		stddev[n] = 1.0;
	}

	int cells = 0;
	int cells_99 = 0;
	int cells_95 = 0;
	int cells_68 = 0;
	
	for(int n = 0; n < num_pts; n++)
	{
		double c1 = att1->GetTuple1(n);
		double c2 = att2->GetTuple1(n);
		if(c1 > att1_min && c1 < att1_max && c2 > att2_min && c2 < att2_max)
		{
			binaryArray->InsertNextTuple1(0.0);
			binary_image[n] = 0;
			cells++;
		}
		else
		{
			binaryArray->InsertNextTuple1(1.0);
			binary_image[n] = 1;
		}
	}

	float* g = (float*)malloc(sizeof(float)*num_pts);
	float* h = (float*)malloc(sizeof(float)*num_pts);
	float* s = (float*)malloc(sizeof(float)*num_pts);	
	
	#pragma omp parallel for
	for(int n = 0; n < num_pts; n++)
	{
		// For each point. Scan the x-axis. 		
		int idx[3];
		GetLogicalPointIndex(idx, n, dims);	
		// we are going to scan x-axis for fixed values of idx[1], idx[2]	
		float pos1;
		pos1 = xc[idx[0]];	


		float min_dist = pow((xc[dims[0]-1] - xc[0]), 2.0);		

		for(int i = 0; i < dims[0]; i++)
		{
			int index = idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+i;
			if(binary_image[index] == 0)
			{
				float pos2;
				pos2 = xc[i];

				float dist = pow((pos2 - pos1), 2.0);
//				float dist = compute_distance(pos1, pos2);
				if(dist < min_dist)
				{
					min_dist = dist;
				}
			}
		}
		
		g[n] = min_dist;
	}
		
	#pragma omp parallel for
	for(int n = 0; n < num_pts; n++)
	{
		// For each point. Scan the x-axis. 		
		int idx[3];
		GetLogicalPointIndex(idx, n, dims);	
		// we are going to scan x-axis for fixed values of idx[1], idx[2]	
		float pos1;
		pos1 = yc[idx[1]];	

		float min_dist; // = pow((yc[dims[1]-1] - yc[0]), 2.0) + pow((xc[dims[0]-1] - xc[0]), 2.0);		

		for(int i = 0; i < dims[1]; i++)
		{
			int index = idx[2]*dims[0]*dims[1]+i*dims[0]+idx[0];
			float pos2;
			pos2 = yc[i];

			float dist = pow((pos2 - pos1), 2.0) + g[index];
			if(i == 0)
			{
				min_dist = dist;
			}
//				float dist = compute_distance(pos1, pos2);
			if(dist < min_dist)
			{
				min_dist = dist;
			}
		}
		
		h[n] = min_dist;
	}
		
	#pragma omp parallel for
	for(int n = 0; n < num_pts; n++)
	{
		// For each point. Scan the x-axis. 		
		int idx[3];
		GetLogicalPointIndex(idx, n, dims);	
		// we are going to scan x-axis for fixed values of idx[1], idx[2]	
		float pos1;
		pos1 = zc[idx[2]];	

		float min_dist; // = pow((yc[dims[1]-1] - yc[0]), 2.0) + pow((yc[dims[1]-1] - yc[0]), 2.0) + pow((xc[dims[0]-1] - xc[0]), 2.0);		

		for(int i = 0; i < dims[2]; i++)
		{
			int index = i*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
			float pos2;
			pos2 = zc[i];

			float dist = pow((pos2 - pos1), 2.0) + h[index];
//			float dist = compute_distance(pos1, pos2);
      if(i == 0)
      {
        min_dist = dist;
      }
			if(dist < min_dist)
			{
				min_dist = dist;
			}
		}
		
		s[n] = min_dist;
	}


	for(int n = 0; n < num_pts; n++)
	{
		distArray->InsertNextTuple1(sqrt(s[n]));
	}

	cout << "Number of points within interval: " << cells << endl;
	cout << "Number of points within 99% confidence interval: " << cells_99 << endl;
	cout << "Number of points within 95% confidence interval: " << cells_95 << endl;
	cout << "Number of points within 68% confidence interval: " << cells_68 << endl;
	/* Write binary_image data as a scalar field to a vtk data set and output it. */
  vtkSmartPointer<vtkDoubleArray> xCoords =
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> yCoords =
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> zCoords =
    vtkSmartPointer<vtkDoubleArray>::New();

  for(int i = 0; i < dims[0]; i++)
  xCoords->InsertNextValue(xc[i]);

  for(int i = 0; i < dims[1]; i++)
  yCoords->InsertNextValue(yc[i]);

  for(int i = 0; i < dims[2]; i++)
  zCoords->InsertNextValue(zc[i]);

  vtkSmartPointer<vtkDataSetWriter> writer = 
    vtkSmartPointer<vtkDataSetWriter>::New();

  vtkSmartPointer<vtkRectilinearGrid> outputGrid =
      vtkSmartPointer<vtkRectilinearGrid>::New();

  outputGrid->SetDimensions(dims[0], dims[1], dims[2]);
  outputGrid->SetXCoordinates(xCoords);
  outputGrid->SetYCoordinates(yCoords);
  outputGrid->SetZCoordinates(zCoords);

//  outputGrid->GetPointData()->AddArray(binaryArray);
  outputGrid->GetPointData()->AddArray(distArray);
  outputGrid->GetPointData()->AddArray(distArray_99);
  outputGrid->GetPointData()->AddArray(distArray_95);
  outputGrid->GetPointData()->AddArray(distArray_68);

  writer->SetFileName(argv[2]);
  writer->SetInputData(outputGrid);
  writer->SetFileTypeToASCII();
  writer->Write();

}
