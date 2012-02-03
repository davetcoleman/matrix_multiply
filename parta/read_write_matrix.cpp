/* Read_Write_Matrix Driver
   CSCI 5576 High Performance Scientific Computing
   Dave Coleman | david.t.coleman@colorado.edu
   2/2/2012

   Each matrix will be described in a file that is one of two formats: either a flat text file or
   an HDF5 file. You will need to be able to read and write both file types. Write a driver
   program that accepts two arguments and call this program read_write_matrix. We will
   pass two filenames to your program. The first file contains a matrix in one of the two
   formats (*.txt or *.hdf5). Read this matrix and write it out to the second file. We will mix
   and match formats.
*/

#include <iostream>
#include <fstream>
#include <istream>
#include <iterator>
#include <algorithm>
#include <vector>
#include "hdf5.h"

using namespace std;

/* Prototypes -------------------------------------------------------- */
namespace hdf5
{
	void hdf5_to_txt(const std::string &file_in, const std::string &file_out);	
	void read(const std::string &filename, std::vector<int> &data);
	void write(const std::string &filename, std::vector<int> &data);
};

int main(int argc, char ** argv)
{
	// Check that 2 file names have been inputted
	if(argc < 3)
	{
		std::cout << "Insufficient file name inputs provided" << std::endl;
		return EXIT_FAILURE;
	}

	// Read in values
	std::string file_in = argv[1];
	std::string file_out = argv[2];

	// Determine what file type is the input file
	// TODO: this assumes the file has a period in it. safe for this assignment
	if(file_in.substr(file_in.find_last_of(".") + 1) == "txt") {
		std::cout << file_in << " is a TXT file " << std::endl;
	}else if(file_in.substr(file_in.find_last_of(".") + 1) == "hdf5") {

		// INPUT = HDF5		
		if(file_out.substr(file_out.find_last_of(".") + 1) == "txt") {
			// OUTPUT = TXT
			std::cout << "HDF5 to TXT" << std::endl;
			hdf5::hdf5_to_txt(file_in, file_out);
		}else if(file_out.substr(file_out.find_last_of(".") + 1) == "hdf5") {
			// OUTPUT = HDF5
			std::cout << "HDF5 to HDF5" << std::endl;
			
			fstream f(file_in.c_str(), fstream::in|fstream::binary);
			f << noskipws;
			istream_iterator<unsigned char> begin(f);
			istream_iterator<unsigned char> end;

			fstream f2(file_out.c_str(), fstream::out|fstream::trunc|fstream::binary);
			ostream_iterator<char> begin2(f2);

			copy(begin, end, begin2);
		}
	}
	
	/*
	  std::vector<int> values(30,10);
	  std::fill (values.begin(),values.begin()+10,5); 
	
	  std::string filename = "data-1d.hdf5";
	
	  hdf5::write(filename,values);


	*/
    return 0;
}


namespace hdf5
{

	void hdf5_to_txt(const std::string &file_in, const std::string &file_out)
	{
		hid_t file_id, dataset_id, space_id; 
		herr_t status;
	
		file_id = H5Fopen(file_in.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		dataset_id = H5Dopen(file_id, "DATASET", H5P_DEFAULT);
		space_id = H5Dget_space(dataset_id);

		int length = H5Sget_simple_extent_npoints(space_id);
		
		// Get actual dimentions:
		hsize_t dims, maxdims;
		int dimentions = H5Sget_simple_extent_dims(space_id, &dims, &maxdims);
		std::cout << "Dimensions are " << dimentions << " with dims " << dims
				  << " and maxdims " << maxdims << " and total length " << length << std::endl;
		
		int image[dims][dims];
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &image);

		std::cout << "here" << std::endl;
		
		// Open a file to save to
		std::ofstream outfile;
		outfile.open (file_out.c_str());
		outfile << dims << " " << dims << std::endl;

		// Save to file
		for (int i=0; i<dims; i++) {
			for (int j=0; j<dims; j++)
      			outfile << image[i][j] << " ";
			outfile << std::endl;
		}
		

		// Close hdf5 stuff
		//delete [] image;
		status = H5Sclose(space_id);
		status = H5Dclose(dataset_id);
		status = H5Fclose(file_id);
		
		// Close the txt file
		outfile.close();			   
	}

	void read(const std::string &filename, std::vector<int> &data)
	{
		hid_t file_id, dataset_id, space_id; 
		herr_t status;
	
		file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		dataset_id = H5Dopen(file_id, "DATASET", H5P_DEFAULT);
		space_id = H5Dget_space(dataset_id);

		int length = H5Sget_simple_extent_npoints(space_id);
		int  * image = new int[length];
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, image);

		// Copy back to vector
		data.resize(length);
		for(int i=0; i<length; ++i){
			data[i] = image[i];
		}
		delete [] image;
		
		status = H5Sclose(space_id);
		status = H5Dclose(dataset_id);
		status = H5Fclose(file_id);
	}

	void write(const std::string &filename, std::vector<int> &data)
	{
		// HDF5 handles
		hid_t file_id, dataset_id, space_id, property_id; 
		herr_t status;

		hsize_t  dims[1] = {data.size()};
    
		//Create a new file using the default properties.
		file_id = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

		//Create dataspace.  Setting maximum size to NULL sets the maximum
		//size to be the current size.
		space_id = H5Screate_simple (1, dims, NULL);

		//Create the dataset creation property list, set the layout to compact.
		property_id = H5Pcreate (H5P_DATASET_CREATE);
		status = H5Pset_layout (property_id, H5D_COMPACT);

		// Create the dataset. 
		dataset_id = H5Dcreate (file_id, "DATASET", H5T_STD_I32LE, space_id, H5P_DEFAULT, property_id, H5P_DEFAULT);
   
		//Write the data to the dataset.
		status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

		status = H5Sclose(space_id);
		status = H5Dclose(dataset_id);
		status = H5Fclose(file_id);
		status = H5Pclose(property_id);
	}

};


