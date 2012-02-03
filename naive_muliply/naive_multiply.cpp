/* Naive_Multiply
   CSCI 5576 High Performance Scientific Computing
   Dave Coleman | david.t.coleman@colorado.edu
   2/2/2012

   You will generate three executables that all take three arguments. The first two
   arguments are matrices that you will multiply together. The third will be the name of the
   file that you will write your results to. The block_multiply takes an optional third
   argument, the block size.
*/

#include <iostream>
#include <fstream>
#include <istream>
#include <iterator>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <vector>
#include "hdf5.h"

using namespace std;

struct matrix{
	vector< vector<int> > data;
	int rows;
	int cols;
};
	
/* Prototypes -------------------------------------------------------- */
namespace hdf5
{
	void read_hdf5(const string &filename, matrix &data);
	void read_txt(const string &filename, matrix &data);
	void write_hdf5(const string &filename, matrix &data);	
	void write_txt(const string &filename, matrix &data);
};
void ioMatrix(bool in, const string &filename, matrix &data );

int main(int argc, char ** argv)
{
	// Check that 3 file names have been inputted
	if(argc < 4)
	{
		cout << "Insufficient file name inputs provided" << endl;
		return EXIT_FAILURE;
	}

	// Read in values
	string file_in1 = argv[1];
	string file_in2 = argv[2];	
	string file_out = argv[3];
	
	// Store Matrix 1
	matrix matrix1;
	ioMatrix(true, file_in1, matrix1);

	// Store Matrix 2
	matrix matrix2;
	ioMatrix(true, file_in2, matrix2);

	// Now naive multiply the two matricies
	// First check that # of colums of matrix1 is = to # rows of matrix2
	if( matrix1.cols != matrix2.rows )
	{
		cout << "Error: number of columns of the left matrix is not equal to the number of rows of the right matrix."
			 << endl;
		return EXIT_FAILURE;
	}

	// Create the output matrix
	matrix matrix_out;
	// Rows = rows of left matrix
	matrix_out.rows = matrix1.rows;
	matrix_out.cols = matrix2.cols;
	matrix_out.data.resize(matrix_out.rows);
	
	// Loop through every cell of the output matrix
	for(int i = 0; i < matrix_out.rows; ++i)
	{
		// Resize the # of columns for every row
		matrix_out.data[i].resize(matrix_out.cols);
		for(int j = 0; j < matrix_out.cols; ++j)
		{
			// Do the multiplication:
			
			int total = 0;
			// Loop through every col of left matrix
			for(int k = 0; k < matrix1.cols; ++k)
			{
				// Calculate and add to current cell
				total += matrix1.data[i][k] * matrix2.data[k][j];

				cout << matrix1.data[i][k] << " times " << matrix2.data[k][j] <<
					" - " << i << " " << k << " x " << k << " " << j << endl;
			}

			cout << i << " " << j << " equals " << total << endl;
			
			// Output answer to matrix
			matrix_out.data[i][j] = total;
		}
	}

	// Now write output matrix to chosen filetype
	ioMatrix(false, file_out, matrix_out);
	
    return EXIT_SUCCESS;
}

void ioMatrix(bool in, const string &filename, matrix &data )
{
	// Determine what file type is the input file
	if(filename.substr(filename.find_last_of(".") + 1) == "txt") {

		if(in)
			hdf5::read_txt(filename, data);
		else
			hdf5::write_txt(filename, data);
		
	}else if(filename.substr(filename.find_last_of(".") + 1) == "hdf5") {

		if(in)
			hdf5::read_hdf5(filename, data);
		else
			hdf5::write_hdf5(filename, data);
	}
}

namespace hdf5
{

	void read_hdf5(const string &filename, matrix &data )
	{
		hid_t file_id, dataset_id, space_id; 
		herr_t status;
	
		file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		dataset_id = H5Dopen(file_id, "DATASET", H5P_DEFAULT);
		space_id = H5Dget_space(dataset_id);

		int length = H5Sget_simple_extent_npoints(space_id);
		
		// Get actual dimentions:
		hsize_t dims, maxdims;
		int dimentions = H5Sget_simple_extent_dims(space_id, &dims, &maxdims);
		cout << endl << "Dimensions are " << dimentions << " with dims " << dims
			 << " and maxdims " << maxdims << " and total length " << length << endl;
		
		int image[dims][dims];
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &image);

		data.data.resize(dims);
		data.rows = dims;
		data.cols = dims;
		
		// Save to file
		for (int i=0; i<dims; i++) {
			data.data[i].resize(dims);
			for (int j=0; j<dims; j++) {
      			cout << image[i][j] << " ";
				data.data[i][j] = image[i][j];
			}
			cout << endl;
		}
		cout << endl;

		// Close hdf5 stuff
		//delete image;
		status = H5Sclose(space_id);
		status = H5Dclose(dataset_id);
		status = H5Fclose(file_id);
		
	}
	
	void read_txt(const string &filename, matrix &data)
	{
		// Read the txt file into a vector
		std::ifstream indata(filename.c_str());
		std::string line;
		std::string cell;

		// The first line should have the number of rows and columns
		std::getline(indata, line);
		std::stringstream lineStream(line);
		
		// Number of rows:		
		getline(lineStream,cell,' ');		
		int rows = atoi(cell.c_str());

		// Number of columns:
		getline(lineStream,cell,' ');
		int cols = atoi(cell.c_str());

		// Create vector
		data.data.resize(rows);
		data.rows = rows;
		data.cols = cols;
		
		// Loop through data
		for(int m = 0; m < rows; ++m)
		{
			getline(indata,line);
			std::stringstream  lineStream(line);

			data.data[m].resize(cols);
			
			for(int n = 0; n < cols; ++n)
			{
				// Get next number
				getline(lineStream,cell,' ');
				//cout << cell << " ";

				// Save next number to matrix
				data.data[m][n] = atoi(cell.c_str());
			}
			//cout << endl;
		}
		cout << rows << " --- " << cols << endl;
		
		cout << "data = " << data.data[1][1] << endl;

	}

	void write_txt(const string &filename, matrix &data)	
	{
		// Open a file to save to
		ofstream outfile;
		outfile.open (filename.c_str());
		outfile << data.rows << " " << data.cols << endl;

		// Save to file
		for (int i=0; i < data.rows; i++) {
			for (int j=0; j < data.cols; j++)
      			outfile << data.data[i][j] << " ";
			outfile << endl;
		}
		
		// Close the txt file
		outfile.close();			   
	}

	void write_hdf5(const string &filename, matrix &matrix)
	{
		// HDF5 handles
		hid_t file_id, dataset_id, space_id, property_id; 
		herr_t status;

		int data[matrix.rows][matrix.cols];		
		
		// Loop through data
		for(int m = 0; m < matrix.rows; ++m)
		{
			for(int n = 0; n < matrix.cols; ++n)
			{
				// Save next number to matrix
				data[m][n] = matrix.data[m][n];
			}
			//cout << endl;
		}
		hsize_t  dims[2] = {matrix.rows, matrix.cols};

		cout << "DIMS ARE " << dims[0] << " " << dims[1] << endl;
    
		//Create a new file using the default properties.
		file_id = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

		//Create dataspace.  Setting maximum size to NULL sets the maximum
		//size to be the current size.
		space_id = H5Screate_simple (2, dims, NULL);

		//Create the dataset creation property list, set the layout to compact.
		property_id = H5Pcreate (H5P_DATASET_CREATE);
		//status = H5Pset_layout (property_id, H5D_COMPACT);
		status = H5Pset_layout (property_id, H5D_CONTIGUOUS);

		// Create the dataset. 
		dataset_id = H5Dcreate (file_id, "DATASET", H5T_STD_I32LE, space_id, H5P_DEFAULT, property_id, H5P_DEFAULT);
   
		//Write the data to the dataset.
		status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);

		status = H5Sclose(space_id);
		status = H5Dclose(dataset_id);
		status = H5Fclose(file_id);
		status = H5Pclose(property_id);
	}
	
};


