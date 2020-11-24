#include <iostream>
#include <math.h>
#include <vtkSmartPointer.h>
#include <vtkRectilinearGridReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
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
	int dims[3];
  mesh->GetDimensions(dims);

	int num_fields = atoi(argv[2]);
	int neighborhood_size = atoi(argv[3]);

	// Compute a standard deviation field for each input scalar.
	for(int n = 0; n < num_fields; n++)
	{
		vtkAbstractArray* a1 = mesh->GetPointData()->GetArray(n);
		vtkDoubleArray* att1 = vtkDoubleArray::SafeDownCast(a1);
		std::cout << att1->GetTuple1(0) << std::endl;
		std::string att1_name(mesh->GetPointData()->GetArrayName(n));
		std::stringstream s;
		s << "stddev_" << att1_name;	
	
		vtkSmartPointer<vtkDoubleArray> stddevArray = vtkSmartPointer<vtkDoubleArray>::New();
 		stddevArray->SetName(s.str().c_str());
  	stddevArray->SetNumberOfComponents(1); 
		
		for(int p = 0; p < num_pts; p++)
		{
			double mean = att1->GetTuple1(p);
			double sum = 0;
			int count = 0;	
			int idx[3];
			GetLogicalPointIndex(idx, p, dims);	
			for(int i = idx[0]-neighborhood_size; i <= idx[0]+neighborhood_size; i++)
			{
				for(int j = idx[1]-neighborhood_size; j <= idx[1]+neighborhood_size; j++)
				{
					for(int k = idx[2]-neighborhood_size; k <= idx[2]+neighborhood_size; k++)
					{
						int p_idx[3];
						p_idx[0] = i;
						p_idx[1] = j;
						p_idx[2] = k;
						int index = GetPointIndex(p_idx, dims);
						if(index >= 0 && index < num_pts)
						{
							sum += pow((mean - att1->GetTuple1(index)), 2.0);
							count += 1;
						}
					}
				}
			}
			stddevArray->InsertNextTuple1(sqrt(sum/count));	
		}
		mesh->GetPointData()->AddArray(stddevArray);
	}

  vtkSmartPointer<vtkDataSetWriter> writer = 
    vtkSmartPointer<vtkDataSetWriter>::New();

	writer->SetFileName(argv[4]);
	writer->SetInputData(mesh);
	writer->SetFileTypeToBinary();
  writer->Write();
}
