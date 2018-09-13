/**
 * created on Mar 1 2018
 * author: Christian Kehl
 */
#include "FAMS_ASCII.hpp"
#include <iomanip>
#include <fstream>
#include <limits>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/classification.hpp>

void FamsFile_ASCII::setFilename(std::string filename) {
	_filename = filename;
}

void FamsFile_ASCII::setDimensions(uint* dimensions) {
	if(dimensions!=NULL)
		_dim=dimensions;
}

void FamsFile_ASCII::SetDimensions(unsigned int* numpy_input_dimensions, int field_size) {
	_ndims = field_size;
	if(_dim==NULL) {
		_dim = new uint[_ndims];
		bzero(_dim, _ndims*sizeof(uint));
	} else {
		delete [] _dim;
		_dim = new uint[_ndims];
		bzero(_dim, _ndims*sizeof(uint));
	}
	memcpy(_dim, numpy_input_dimensions, _ndims*sizeof(uint));
}

void FamsFile_ASCII::setData(void* data, dtype datatype, bool memcopy) {
	if(!memcopy)
		_data=data;
	else {
		clearData();
		_datatype = USHORT;
		long len = 1;
		for(uint i=0; i<_ndims; i++)
			len*=long(_dim[i]);
		_data = new short[len];
		memcpy(_data, data, len*sizeof(short));
	}
}

void FamsFile_ASCII::setDataByUShort(ushort* numpy_input_data, long field_size) {
	clearData();
	_datatype = USHORT;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	_data = new ushort[fsize];
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(_data, numpy_input_data, fieldSize*sizeof(ushort));
}

void FamsFile_ASCII::setDataByUInt(uint* numpy_input_data, long field_size) {
	clearData();
	_datatype = USHORT;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	_data = new ushort[fsize];
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	//memcpy(_data, numpy_input_data, fieldSize*sizeof(uint));
	for(uint i=0; i<fieldSize; i++) {
		((ushort*)_data)[i]=numpy_input_data[i];
	}
}

void FamsFile_ASCII::getDataAsUShort(ushort* numpy_output_data, long field_size) {
	if(_dim==NULL)
		return;
	if(_data==NULL)
		return;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(numpy_output_data, (ushort*)_data, fieldSize*sizeof(ushort));
}

void FamsFile_ASCII::getDataAsUInt(uint* numpy_output_data, long field_size) {
	if(_dim==NULL)
		return;
	if(_data==NULL)
		return;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	//memcpy(numpy_output_data, (uint*)_data, fieldSize*sizeof(uint));
	for(uint i=0; i<fieldSize; i++) {
		numpy_output_data[i]= uint(((ushort*)_data)[i]);
	}
}

void FamsFile_ASCII::getDimensions(unsigned int* numpy_output_dimensions, int field_size) {
	if(_dim==NULL)
		return;
	int fieldSize = (field_size<=_ndims) ? field_size : _ndims;
	memcpy(numpy_output_dimensions, _dim, fieldSize*sizeof(unsigned int));
}

void FamsFile_ASCII::readFile(void) {
	//TODO: to be coded
	if(_filename.empty()) {
#ifdef DEBUG
		std::cerr << "No file information given." << std::endl;
#endif
	}
	int slashpos = _filename.find_last_of("/");
	int dotpos = _filename.find_last_of(".");
	std::string inputExtension = _filename.substr(dotpos+1, _filename.npos);
	if(inputExtension.find("txt")==inputExtension.npos) {
#ifdef DEBUG
		std::cout << "no spectral data file as input - EXITING ..." << std::endl;
#endif
		return;
	}
	std::string basepath = _filename.substr(0,slashpos);
	std::string basenameWend = _filename.substr(slashpos+1, _filename.npos);
	dotpos = basenameWend.find_last_of(".");
	//std::cout << basenameWend << std::endl;
	std::string basename = basenameWend.substr(0,dotpos);
	std::string dataFilename = basepath+"/"+basename+"."+"txt";
#ifdef DEBUG
	std::cout << "Basepath: " << basepath << std::endl;
	std::cout << "Basename (with ending): " << basenameWend << std::endl;
	std::cout << "Basename (no ending): " << basename << std::endl;
	std::cout << "filename: " << dataFilename << std::endl;
#endif

	std::ifstream dataFile;
	dataFile.open(dataFilename.c_str(), std::ios::in);
	std::vector<std::string> strs;
	strs.reserve(256);
	std::string line;
	int linenum=-1;
	long dataEntries = 0;
	ushort dataDims = 0;
	long dataEntry = 0;
	if(dataFile.is_open())
	{
		linenum++;
		std::getline(dataFile,line);
		while(dataFile.good()==true)
		{
			strs.clear();
			if(line.empty()) {
				if(dataFile.eof()==true)
					break;
				continue;
			}

			std::string trimmed = boost::algorithm::trim_copy(line);
			boost::split(strs, trimmed, boost::is_any_of(" "));
			for(uint k = 0; k < strs.size(); k++)
				strs.at(k) = boost::algorithm::trim_copy(strs.at(k));

			if(linenum==0) {
				assert(strs.size()>=2);
				dataEntries = long(atoi(strs.at(0).c_str()));
				dataDims = ushort(atoi(strs.at(1).c_str()));
				if(_dim!=NULL) {
					long fsize=1;
					for(int i=0; i<_ndims; i++)
						fsize*=long(_dim[i]);
					if(fsize!=(dataEntries*dataDims)) {
						std::cerr << "Stored- and defined dimensions don't match; using stored dimensions." << std::endl;
						_ndims=2;
						delete [] _dim;
						_dim = new uint[_ndims];
						_dim[0] = dataEntries;
						_dim[1] = dataDims;
					}
				} else {
					_ndims=2;
					_dim = new uint[_ndims];
					_dim[0] = dataEntries;
					_dim[1] = dataDims;
				}
				assert(_dim[_ndims-1]==dataDims);
				long fsize=1;
				for(int i=0; i<_ndims; i++)
					fsize*=long(_dim[i]);

				clearData();
				_data = new ushort[fsize];
			} else {
				assert(strs.size()==dataDims);
				long sIndex = dataEntry*dataDims;
				for(ushort dIdx = 0; dIdx < strs.size(); dIdx++) {
					((ushort*)_data)[sIndex+long(dIdx)]=ushort(atoi(strs[dIdx].c_str()));
				}
				dataEntry++;
			}
			std::getline(dataFile,line);
			if(dataFile.eof()==true)
				break;
		}
#ifdef DEBUG
		std::cout << "File reading operation done." << std::endl;
#endif
		dataFile.close();

	}
	else
	{
#ifdef DEBUG
		std::cout << "Unable to open multivariate data file." << std::endl;
#endif
		return;
	}
}

void FamsFile_ASCII::writeFile(void) {
	if(_dim==NULL) {
#ifdef DEBUG
		std::cerr << "No dimensionality information provided." << std::endl;
#endif
		return;
	}
	if(_data==NULL) {
#ifdef DEBUG
		std::cerr << "No data provided." << std::endl;
#endif
		return;
	}

	int slashpos = _filename.find_last_of("/");
	int dotpos = _filename.find_last_of(".");
	std::string inputExtension = _filename.substr(dotpos+1, _filename.npos);
	if(inputExtension.find("txt")==inputExtension.npos) {
		std::cout << "no volume .mhd file - EXITING ..." << std::endl;
		return;
	}
	std::string basepath = _filename.substr(0,slashpos);
	std::string basenameWend = _filename.substr(slashpos+1, _filename.npos);
	dotpos = basenameWend.find_last_of(".");
	//std::cout << basenameWend << std::endl;
	std::string basename = basenameWend.substr(0,dotpos);
	std::string dataFilename = basepath+"/"+basename+"."+"txt";
	//std::cout << basepath+"/"+basename << std::endl;


	std::ofstream dataFile;
	dataFile.open(dataFilename.c_str(), std::ios::out|std::ios::trunc);
	if(dataFile.is_open())
	{
		long dEntries = 1;
		for(int i=0; i < (_ndims-1); i++)
			dEntries*=_dim[i];
		uint dDims = _dim[_ndims-1];
		dataFile << dEntries << " " << dDims << std::endl;
		for(long i=0; i<dEntries; i++) {
			for(uint j=0; j<dDims; j++) {
				//TODO: make this dataorder-adaptive!!! TEST!!
				if(_dataorder==ROW_MAJOR)
					dataFile << ((ushort*)_data)[i*dDims+j];
				else if(_dataorder==COLUMN_MAJOR)
					dataFile << ((ushort*)_data)[j*dEntries+i];
				if(j<(dDims-1))
					dataFile << " ";
			}
			dataFile << std::endl;
		}
		dataFile.close();
	}
	else {
		std::cout << "Cannot open multivariate data file for writing." << std::endl;
	}
}






