#include <iostream>
#include <math.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGrid.h>
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
	vtkSmartPointer<vtkStructuredGridReader> mean_reader = 
    vtkSmartPointer<vtkStructuredGridReader>::New();

	mean_reader->SetFileName(argv[1]);
	mean_reader->Update();

	vtkSmartPointer<vtkStructuredGrid> mesh = 
    vtkSmartPointer<vtkStructuredGrid>::New();
	mesh = mean_reader->GetOutput();
	
	int num_pts = mesh->GetNumberOfPoints();
		
	vtkSmartPointer<vtkStructuredGridReader> dev_reader = 
    vtkSmartPointer<vtkStructuredGridReader>::New();
	dev_reader->SetFileName(argv[2]);
	dev_reader->Update();

	vtkSmartPointer<vtkStructuredGrid> std_mesh = 
    vtkSmartPointer<vtkStructuredGrid>::New();
	std_mesh = dev_reader->GetOutput();

	cout << "Has array : " << mesh->GetPointData()->GetArrayName(0) << endl;
	cout << "Has array : " << std_mesh->GetPointData()->GetArrayName(0) << endl;

	vtkAbstractArray* a1 = mesh->GetPointData()->GetAbstractArray(argv[6]);
	vtkAbstractArray* a2 = std_mesh->GetPointData()->GetAbstractArray(argv[7]);

	vtkFloatArray* mean = vtkFloatArray::SafeDownCast(a1);
	vtkFloatArray* dev = vtkFloatArray::SafeDownCast(a2);

	int dims[3];

  mesh->GetDimensions(dims);

	cout << num_pts << endl;

// Output the range of each attribute 

	double range[2]; 
	mean->GetRange(range);
	cout << "Att1 range " << range[0] << " to " << range[1] << endl;
	
	float att1_min = atof(argv[4]);
	float att1_max = atof(argv[5]);
	
	cout << "Interval start: " << att1_min << endl;
	cout << "Interval end: " << att1_max << endl;

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
	
	float *xc = (float*)malloc(sizeof(float)*dims[0]);
	float *yc = (float*)malloc(sizeof(float)*dims[1]);
	float *zc = (float*)malloc(sizeof(float)*dims[2]);

	for(int i = 0; i < dims[0]; i++)
	{
		double pt[3];
		mesh->GetPoint(i, pt);
		xc[i] = pt[0];
	}
	
	for(int i = 0; i < dims[1]; i++)
	{
		int index = dims[0]*i;
		double pt[3];
		mesh->GetPoint(index, pt);
		yc[i] = pt[1];
	}

	for(int i = 0; i < dims[2]; i++)
	{
		int index = dims[0]*dims[1]*i;
		double pt[3];
		mesh->GetPoint(index, pt);
		zc[i] = pt[2];
	}

	int *binary_image = (int*)malloc(sizeof(int)*num_pts);
	float *minvals = (float*)malloc(sizeof(float)*num_pts);
	float *maxvals = (float*)malloc(sizeof(float)*num_pts);

	vtkSmartPointer<vtkFloatArray> binaryArray = vtkSmartPointer<vtkFloatArray>::New();
	binaryArray->SetName("binary");
	binaryArray->SetNumberOfComponents(1); 
	
	vtkSmartPointer<vtkFloatArray> distArray = vtkSmartPointer<vtkFloatArray>::New();
	distArray->SetName("distance");
	distArray->SetNumberOfComponents(1); 
	
	vtkSmartPointer<vtkFloatArray> distArray_50 = vtkSmartPointer<vtkFloatArray>::New();
	distArray_50->SetName("distance50");
	distArray_50->SetNumberOfComponents(1); 
	
	vtkSmartPointer<vtkFloatArray> distArray_95 = vtkSmartPointer<vtkFloatArray>::New();
	distArray_95->SetName("distance95");
	distArray_95->SetNumberOfComponents(1); 
	
	vtkSmartPointer<vtkFloatArray> distArray_68 = vtkSmartPointer<vtkFloatArray>::New();
	distArray_68->SetName("distance68");
	distArray_68->SetNumberOfComponents(1); 

	int cells = 0;
	int cells_50 = 0;
	int cells_95 = 0;
	int cells_68 = 0;

	for(int n = 0; n < num_pts; n++)
	{
		double c1 = mean->GetTuple1(n);
		if(c1 > att1_min && c1 < att1_max) 
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

// Computing the confidence interval isosurfaces // 

// Confidence interval = 50. Z = 0.67449
  double Z = 0.67449;
	for(int n = 0; n < num_pts; n++)
  {
    double c1 = mean->GetTuple1(n);

		minvals[n] = c1 - (Z*dev->GetTuple1(n));
		maxvals[n] = c1 + (Z*dev->GetTuple1(n));

	}
		
	
	for(int n = 0; n < num_pts; n++)
  {
		if((minvals[n] <= att1_max && maxvals[n] >= att1_min))
    {
      binary_image[n] = 0;
      cells_50++;
    }
    else
    {
      binary_image[n] = 1;
    }
  }

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
		distArray_50->InsertNextTuple1(sqrt(s[n]));
	}

// Confidence interval = 95. Z = 1.95996
  Z = 1.95996;
	for(int n = 0; n < num_pts; n++)
  {
    double c1 = mean->GetTuple1(n);

		minvals[n] = c1 - Z*dev->GetTuple1(n);
		maxvals[n] = c1 + Z*dev->GetTuple1(n);
	}
		
	
	for(int n = 0; n < num_pts; n++)
  {
		if((minvals[n] <= att1_max && maxvals[n] >= att1_min)) 
    {
      binary_image[n] = 0;
      cells_95++;
    }
    else
    {
      binary_image[n] = 1;
    }
  }

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
		distArray_95->InsertNextTuple1(sqrt(s[n]));
	}

// Confidence interval = 68. Z = 1.0
  Z = 1.0;
	for(int n = 0; n < num_pts; n++)
  {
    double c1 = mean->GetTuple1(n);

		minvals[n] = c1 - Z*dev->GetTuple1(n);
		maxvals[n] = c1 + Z*dev->GetTuple1(n);
	}
		
	
	for(int n = 0; n < num_pts; n++)
  {
		if((minvals[n] <= att1_max && maxvals[n] >= att1_min)) 
    {
      binary_image[n] = 0;
      cells_68++;
    }
    else
    {
      binary_image[n] = 1;
    }
  }

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
		distArray_68->InsertNextTuple1(sqrt(s[n]));
	}


	cout << "Number of points within interval: " << cells << endl;
	cout << "Number of points within 50% confidence interval: " << cells_50 << endl;
	cout << "Number of points within 95% confidence interval: " << cells_95 << endl;
	cout << "Number of points within 68% confidence interval: " << cells_68 << endl;
	

	/* Write binary_image data as a scalar field to a vtk data set and output it. */

  vtkSmartPointer<vtkFloatArray> xCoords =
    vtkSmartPointer<vtkFloatArray>::New();
  vtkSmartPointer<vtkFloatArray> yCoords =
    vtkSmartPointer<vtkFloatArray>::New();
  vtkSmartPointer<vtkFloatArray> zCoords =
    vtkSmartPointer<vtkFloatArray>::New();

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
  outputGrid->GetPointData()->AddArray(distArray_50);
  outputGrid->GetPointData()->AddArray(distArray_95);
  outputGrid->GetPointData()->AddArray(distArray_68);

  writer->SetFileName(argv[3]);
  writer->SetInputData(outputGrid);
  writer->SetFileTypeToBinary();
  writer->Write();
}
