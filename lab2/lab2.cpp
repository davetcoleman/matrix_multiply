/* Lab2
   CSCI 5576 High Performance Scientific Computing
   Dave Coleman | david.t.coleman@colorado.edu

   2/2/2012

   Note - this one file includes all the source code for all the executables. The Make file generates
   seperate versions for each function.

   You will generate three executables that all take three arguments. The first two
   arguments are matrices that you will multiply together. The third will be the name of the
   file that you will write your results to. The block_multiply takes an optional third
   argument, the block size.
*/

#include <iostream>
#include <fstream>
#include <istream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <sys/time.h>
#include "hdf5.h"
#include <math.h>

#if defined USE_INTEL
#include "mkl_cblas.h"
#include <stdio.h>
#endif

using namespace std;

//-------------------------------------------------------------------------------------------
// Matrix Structure
//-------------------------------------------------------------------------------------------
struct matrix{
	double *data;
	int rows;
	int cols;
};
	
//-------------------------------------------------------------------------------------------
// Prototypes
//-------------------------------------------------------------------------------------------
void ioMatrix(bool in, const string &filename, matrix &data );
double get_time();
void naive_multiply( matrix &matrix1, matrix &matrix2, matrix &matrix_out);
void block_multiply(matrix &matrix1, matrix &matrix2, matrix &matrix_out, int block_size);
void intel_multiply(matrix &matrix1, matrix &matrix2, matrix &matrix_out);
void printMatrix(matrix &data);

namespace hdf5
{
	void read_hdf5(const string &filename, matrix &data);
	void read_txt(const string &filename, matrix &data);
	void write_hdf5(const string &filename, matrix &data);	
	void write_txt(const string &filename, matrix &data);
};

//-------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------
int main(int argc, char ** argv)
{
	// Check that 3 file names have been inputted
	if(argc < 4)
	{
		if(argc < 3)
		{
			cout << "Insufficient file name inputs provided" << endl;
			return EXIT_FAILURE;
		}
		else // do a converstion!
		{
			cout << "Converting files" << endl;
			// Read in values
			string file_in1 = argv[1];
			string file_out = argv[2];
	
			// Read and write matrix
			matrix matrix1;
			ioMatrix(true, file_in1, matrix1);
			ioMatrix(false, file_out, matrix1);

			// Done with program
			return EXIT_SUCCESS;
		}
	}

	// Read in values
	string file_in1 = argv[1];
	string file_in2 = argv[2];	
	string file_out = argv[3];

	int block_size = 32;
	// Check for 4th parameter
	if(argc > 4)
		block_size = atoi(argv[4]);
	
	// Store Matrix 1
	matrix matrix1;
	ioMatrix(true, file_in1, matrix1);
	//printMatrix(matrix1);
	
	// Store Matrix 2
	matrix matrix2;
	ioMatrix(true, file_in2, matrix2);

	
	// Now naive multiply the two matricies -----------------------------------
	
	// First check that # of colums of matrix1 is = to # rows of matrix2
	if( matrix1.cols != matrix2.rows )
	{
		cout << "Error: number of columns of the left matrix is not equal to the number of "
			 << "rows of the right matrix."
			 << endl;
		return EXIT_FAILURE;
	}

	// Create the output matrix
	matrix matrix_out;

	// Do the multiplication:
	int number_runs = 1; // increase this to do benchmarking
	double start, end; // for benchmarking
	start = get_time();
	
	for(int test = 0; test < number_runs; ++test)
	{
		if(number_runs > 1)
			cout << "Starting run " << number_runs << endl;

		// The Makefile tells the compiler which matrix multiply to use.
		// Naive method is chosen by default
#if defined USE_BLOCK
		cout << "Block Matrix Multiply" << endl;
   		block_multiply(matrix1, matrix2, matrix_out, block_size);
#elif defined USE_INTEL
		cout << "Intel Matrix Multiply" << endl;
   		intel_multiply(matrix1, matrix2, matrix_out);
#else
		cout << "Naive Matrix Multiply" << endl;		
   		naive_multiply(matrix1, matrix2, matrix_out);
#endif
		
	}
	end = get_time();

	//	if(number_runs > 1)
	//{
		cout << "Matrix Multiply ran in " << (end-start) << " seconds " << endl;
		cout << endl;
		//}	

	//printMatrix(matrix_out);
	
	// Now write output matrix to chosen filetype
	//ioMatrix(false, file_out, matrix_out);
	
    return EXIT_SUCCESS;
}
//-------------------------------------------------------------------------------------------
// Multiply Method 1
//-------------------------------------------------------------------------------------------
void naive_multiply(matrix &matrix1, matrix &matrix2, matrix &matrix_out)
{
	// Rows = rows of left matrix
	matrix_out.rows = matrix1.rows;
	matrix_out.cols = matrix2.cols;
	matrix_out.data = new double[matrix_out.rows * matrix_out.cols];

	// Loop through every cell of the output matrix
	for(int i = 0; i < matrix_out.rows; ++i)
	{
		for(int j = 0; j < matrix_out.cols; ++j)
		{
			// initalize result value to zero
			//cout << "writing to location " << i*matrix_out.rows + j << endl;
			matrix_out.data[i*matrix_out.rows + j] = 0;
					
			// Loop through every col of left matrix
			for(int k = 0; k < matrix1.cols; ++k)
			{
				//cout << "data " << matrix1.data[i*matrix1.rows + k] << " - " << matrix2.data[k*matrix2.rows + j];
				//cout << " ========== " << i << " " << k << " x " << k << " " << j;
				
				// Calculate and add to current cell
				matrix_out.data[i*matrix_out.rows + j] += matrix1.data[i*matrix1.rows + k] *
					matrix2.data[k*matrix2.rows + j];
				
				//cout << " ========== " << matrix_out.data[i*matrix_out.rows + j] << endl;
			}
		}
	}
}
//-------------------------------------------------------------------------------------------
// Multiply Method 2
//-------------------------------------------------------------------------------------------
void block_multiply(matrix &matrix1, matrix &matrix2, matrix &matrix_out, int block_size)
{
	//cout << "Starting Block Multiply" << endl;
	//cout << "Block size: " << block_size << endl;

	// Rows = rows of left matrix
	matrix_out.rows = matrix1.rows;
	matrix_out.cols = matrix2.cols;
	matrix_out.data = new double[matrix_out.rows * matrix_out.cols];
	
	// Initialize all values to zero
	for(int i = 0; i < matrix_out.rows * matrix_out.cols; ++i)
   		matrix_out.data[i] = 0;
		   
	int row_blocks = int(ceil( double(matrix_out.rows) / block_size ));
	int col_blocks = int(ceil( double(matrix_out.cols) / block_size ));	

	//cout << "Block rows: " << row_blocks << endl;
	//cout << "Block cols: " << col_blocks << endl;

	// Loop through every cell of the output matrix
	for(int ii = 0; ii < row_blocks; ++ii)
	{
		//cout << "for ii= " << ii<< endl;
		for(int jj = 0; jj < col_blocks; ++jj)
		{
			//cout << "for jj= " << jj<< endl;
			for(int kk = 0; kk < row_blocks; ++kk)
			{
				//cout << "for kk= " << kk<< endl;
				for(int i = ii*block_size; i < min(ii*block_size + block_size, matrix_out.rows); ++i)
				{
					//cout << "for i= " << i<< endl;
					for(int j = jj*block_size; j < min(jj*block_size + block_size, matrix_out.cols); ++j)
					{
						//cout << "for j= " << j<< endl;
						
						// Loop through every col of left matrix
						for(int k = kk*block_size; k < min(kk*block_size + block_size, matrix_out.cols); ++k)
						{
							/*cout << "for k= " << k << endl;
							
							  cout << "data " << matrix1.data[i*matrix1.rows + k] << " - ";
							  cout            << matrix2.data[k*matrix2.rows + j];
							  cout << " ========== " << i << " " << k << " x " << k << " " << j;
							  cout << " -> " << i << " " << j;
							*/
							
							// Calculate and add to current cell
							matrix_out.data[i*matrix_out.rows + j] += matrix1.data[i*matrix1.rows + k] *
								matrix2.data[k*matrix2.rows + j];
				
							//cout << " ========== " << matrix_out.data[i*matrix_out.rows + j] << endl;
						
						}
					}
				}
			}
		}
	}
	
	//printMatrix(matrix_out);
	
}
//-------------------------------------------------------------------------------------------
// Multiply Method 3
//-------------------------------------------------------------------------------------------
void intel_multiply(matrix &matrix1, matrix &matrix2, matrix &matrix_out)
{
#if defined USE_INTEL	

	// Rows = rows of left matrix
	matrix_out.rows = matrix1.rows;
	matrix_out.cols = matrix2.cols;
	matrix_out.data = new double[matrix_out.rows * matrix_out.cols];
	
	// Initialize all values to zero
	for(int i = 0; i < matrix_out.rows * matrix_out.cols; ++i)
   		matrix_out.data[i] = 0;
	
	int lda = matrix1.rows;
	int ldb = matrix2.rows;
	int ldc = matrix_out.rows; 
	
	//                                                      M             N             K           alpha
	cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, matrix1.rows, matrix2.cols, matrix1.cols, 1.0,
				 matrix1.data, lda, matrix2.data, ldb, 0.0, matrix_out.data, ldc);
	//           A             lda  B             ldb  beta C                ldc
	
#endif
}
//-------------------------------------------------------------------------------------------
// Read and write matrix depending on file type
//-------------------------------------------------------------------------------------------
void ioMatrix(bool in, const string &filename, matrix &data )
{
	int number_runs = 1; // increase this to do benchmarking
	double start, end; // for benchmarking
	start = get_time();
	
	for(int test = 0; test < number_runs; ++test)
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
	end = get_time();

	if(number_runs > 1)
	{
		cout << "IO Test ran in " << (end-start) << " seconds " << endl;
		cout << "File type: " << filename.substr(filename.find_last_of(".") + 1) << endl;
		if(in)
			cout << "Input" << endl;
		else
			cout << "Output" << endl;
		cout << endl;
	}
	
}

//-------------------------------------------------------------------------------------------
// For benchmarking
//-------------------------------------------------------------------------------------------
double get_time()
{
	struct timeval t;
	struct timezone tzp;
	gettimeofday(&t, &tzp);
	return t.tv_sec + t.tv_usec*1e-6;
}
//-------------------------------------------------------------------------------------------
// For Testing
//-------------------------------------------------------------------------------------------
void printMatrix(matrix &data)
{
	cout << endl << "Printing Matrix with " << data.rows << " rows and " << data.cols << " cols." << endl;
	for(int i = 0; i < data.rows; ++i)
	{
		for(int j = 0; j < data.cols; ++j)
		{
			cout << data.data[i*data.rows + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}




namespace hdf5
{
	//-------------------------------------------------------------------------------------------
	// Read hdf5
	//-------------------------------------------------------------------------------------------
	void read_hdf5(const string &filename, matrix &data )
	{
		hid_t file_id, dataset_id, space_id; 
		herr_t status;
	
		file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		dataset_id = H5Dopen(file_id, "DATASET", H5P_DEFAULT);
		space_id = H5Dget_space(dataset_id);

		//int length = H5Sget_simple_extent_npoints(space_id);
		
		// Get actual dimentions:
		hsize_t dims;
   		hsize_t maxdims;
		//int dimentions = H5Sget_simple_extent_dims(space_id, &dims, &maxdims);
		cout << "here" << endl;		
		H5Sget_simple_extent_dims(space_id, &dims, &maxdims);
		cout << "here2" << endl;
		cout << dims << " " << maxdims << " output" << endl;
		
		//cout << endl << "Dims " << dims[0] << " " << dims[1] << " and maxdims " << &maxdims << endl;
		
		int image[ dims ][ dims ];  // TODO: make this work for non-square matricies
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &image);

		// TODO: row col differenciation
		data.data = new double[dims*dims];	
		data.rows = dims;
		data.cols = dims;
		
		// Convert to basic array
		for (int i=0; i < data.rows; i++) {
			for (int j=0; j < data.cols; j++) {
				data.data[i*data.rows + j] = image[i][j];
			}
		}
		
		
		// Close hdf5 stuff
		//for(unsigned int i = 0; i < dims; ++i)
		//	delete[] image[i];
		//		delete[] image;
		status = H5Sclose(space_id);
		status = H5Dclose(dataset_id);
		status = H5Fclose(file_id);
		
	}
	
	//-------------------------------------------------------------------------------------------
	// Read txt
	//-------------------------------------------------------------------------------------------	
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
		data.data = new double[rows*cols];
		
		//data.data.resize(rows);
		data.rows = rows;
		data.cols = cols;
		
		// Loop through data
		for(int i = 0; i < rows; ++i)
		{
			getline(indata,line);
			std::stringstream  lineStream(line);

			//data.data[m].resize(cols);
			
			for(int j = 0; j < cols; ++j)
			{
				// Get next number
				getline(lineStream,cell,' ');
				//cout << cell << " ";

				// Save next number to matrix
				data.data[i*rows + j] = atoi(cell.c_str());
			}
			//cout << endl;
		}

	}
	//-------------------------------------------------------------------------------------------
	// Write txt
	//-------------------------------------------------------------------------------------------
	void write_txt(const string &filename, matrix &data)	
	{
		// Open a file to save to
		ofstream outfile;
   		outfile.open (filename.c_str());
		outfile << data.rows << " " << data.cols << endl;

		// Save to file
		for (int i=0; i < data.rows; i++)
		{
			for (int j=0; j < data.cols; j++)
			{
				//cout << i << " " << j << endl;
				outfile << data.data[i*data.rows + j] << " ";
			}
			outfile << endl;
		}
		
		// Close the txt file
		outfile.close();			   
	}
	//-------------------------------------------------------------------------------------------
	// Write hdf5
	//-------------------------------------------------------------------------------------------
	void write_hdf5(const string &filename, matrix &matrix)
	{
		// HDF5 handles
		hid_t file_id, dataset_id, space_id, property_id; 
		herr_t status;

		int image[matrix.rows][matrix.cols];		

		// Loop through data
		for(int i = 0; i < matrix.rows; ++i)
		{
			for(int j = 0; j < matrix.cols; ++j)
			{
				// Save next number to matrix
				//cout <<  matrix.data[i*matrix.rows + j] << " ";
				image[i][j] = int(matrix.data[i*matrix.rows + j]);
			}
			//cout << endl;
		}
		
		hsize_t  dims[2] = {matrix.rows, matrix.cols};

		//cout << "DIMS ARE " << dims[0] << " " << dims[1] << endl;
    
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
   
		// Write the data to the dataset.
		status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &image);

		// Close H5 stuff
		status = H5Sclose(space_id);
		status = H5Dclose(dataset_id);
		status = H5Fclose(file_id);
		status = H5Pclose(property_id);

		// Delete variable
		//for(int i = 0; i < matrix.cols; ++i)
		//	delete[] image[i];
		//		delete[] image;		
	}
};


