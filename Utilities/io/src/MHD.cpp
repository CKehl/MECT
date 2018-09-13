/**
 * created on Jan 9 2018
 * author: Christian Kehl
 */
#include "MHD.hpp"
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

void MhdFile::setFilename(std::string filename) {
	_filename = filename;
}

void MhdFile::setDimensions(uint* dimensions) {
	if(dimensions!=NULL)
		_dim=dimensions;
}

void MhdFile::SetDimensions(unsigned int* numpy_input_dimensions, int field_size) {
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

void MhdFile::setSpacing(float* spacing) {
	if(spacing!=NULL)
		_spacing = spacing;

}

void MhdFile::SetSpacing(float* numpy_input_spacing, int field_size) {
	int fieldSize = (field_size<=_ndims) ? field_size : _ndims;
	if(_spacing==NULL) {
		_spacing = new float[_ndims];
		bzero(_spacing, _ndims*sizeof(float));
	}
	memcpy(_spacing, numpy_input_spacing, fieldSize*sizeof(float));
}

void MhdFile::SetSpacing(double* numpy_input_spacing, int field_size) {
	int fieldSize = (field_size<=_ndims) ? field_size : _ndims;
	if(_spacing==NULL) {
		_spacing = new float[_ndims];
		bzero(_spacing, _ndims*sizeof(float));
	}
	for(int i=0; i<fieldSize; i++)
		_spacing[i] = float(numpy_input_spacing[i]);
}

void MhdFile::setData(void* data, dtype datatype, bool memcopy) {
	if(!memcopy)
		_data=data;
	else {
		clearData();
		_datatype = datatype;
		long len = 1;
		for(uint i=0; i<(uint)_ndims; i++)
			len*=long(_dim[i]);
		switch(_datatype) {
			case CHAR: {
				_data = new char[len];
				memcpy(_data, data, len*sizeof(char));
				break;
			}
			case UCHAR: {
				_data = new uchar[len];
				memcpy(_data, data, len*sizeof(uchar));
				break;
			}
			case SHORT: {
				_data = new short[len];
				memcpy(_data, data, len*sizeof(short));
				break;
			}
			case USHORT: {
				_data = new ushort[len];
				memcpy(_data, data, len*sizeof(ushort));
				break;
			}
			case INT: {
				_data = new int[len];
				memcpy(_data, data, len*sizeof(int));
				break;
			}
			case UINT: {
				_data = new uint[len];
				memcpy(_data, data, len*sizeof(uint));
				break;
			}
			case LONG: {
				_data = new long[len];
				memcpy(_data, data, len*sizeof(long));
				break;
			}
			case ULONG: {
				_data = new ulong[len];
				memcpy(_data, data, len*sizeof(ulong));
				break;
			}
			case FLOAT: {
				_data = new float[len];
				memcpy(_data, data, len*sizeof(float));
				break;
			}
			case DOUBLE: {
				_data = new double[len];
				memcpy(_data, data, len*sizeof(double));
				break;
			}
			default: {
				break;
			}
		}
	}
}

void MhdFile::setDataAsChar(char* numpy_input_data, long field_size) {
	clearData();
	_datatype = CHAR;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	_data = new char[fsize];
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(_data, numpy_input_data, fieldSize*sizeof(char));
}

void MhdFile::setDataAsUChar(uchar* numpy_input_data, long field_size) {
	clearData();
	_datatype = UCHAR;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	_data = new uchar[fsize];
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(_data, numpy_input_data, fieldSize*sizeof(uchar));
}

void MhdFile::setDataAsShort(short* numpy_input_data, long field_size) {
	clearData();
	_datatype = SHORT;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	_data = new short[fsize];
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(_data, numpy_input_data, fieldSize*sizeof(short));
}

void MhdFile::setDataAsUShort(ushort* numpy_input_data, long field_size) {
	clearData();
	_datatype = USHORT;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	_data = new ushort[fsize];
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(_data, numpy_input_data, fieldSize*sizeof(ushort));
}

void MhdFile::setDataAsInt(int* numpy_input_data, long field_size) {
	clearData();
	_datatype = INT;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	_data = new int[fsize];
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(_data, numpy_input_data, fieldSize*sizeof(int));
}

void MhdFile::setDataAsUInt(uint* numpy_input_data, long field_size) {
	clearData();
	_datatype = UINT;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	_data = new uint[fsize];
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(_data, numpy_input_data, fieldSize*sizeof(uint));
}

void MhdFile::setDataAsLong(long* numpy_input_data, long field_size) {
	clearData();
	_datatype = LONG;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	_data = new long[fsize];
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(_data, numpy_input_data, fieldSize*sizeof(long));
}

void MhdFile::setDataAsULong(ulong* numpy_input_data, long field_size) {
	clearData();
	_datatype = ULONG;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	_data = new ulong[fsize];
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(_data, numpy_input_data, fieldSize*sizeof(ulong));
}

void MhdFile::setDataAsFloat(float* numpy_input_data, long field_size) {
	clearData();
	_datatype = FLOAT;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	_data = new float[fsize];
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(_data, numpy_input_data, fieldSize*sizeof(float));
}

void MhdFile::setDataAsDouble(double* numpy_input_data, long field_size) {
	clearData();
	_datatype = DOUBLE;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	_data = new double[fsize];
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(_data, numpy_input_data, fieldSize*sizeof(double));
}

void MhdFile::getDataAsChar(char* numpy_output_data, long field_size) {
	if(_dim==NULL)
		return;
	if(_data==NULL)
		return;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(numpy_output_data, (char*)_data, fieldSize*sizeof(char));
}

void MhdFile::getDataAsUChar(uchar* numpy_output_data, long field_size) {
	if(_dim==NULL)
		return;
	if(_data==NULL)
		return;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(numpy_output_data, (uchar*)_data, fieldSize*sizeof(uchar));
}

void MhdFile::getDataAsShort(short* numpy_output_data, long field_size) {
	if(_dim==NULL)
		return;
	if(_data==NULL)
		return;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(numpy_output_data, (short*)_data, fieldSize*sizeof(short));
}

void MhdFile::getDataAsUShort(ushort* numpy_output_data, long field_size) {
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

void MhdFile::getDataAsInt(int* numpy_output_data, long field_size) {
	if(_dim==NULL)
		return;
	if(_data==NULL)
		return;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(numpy_output_data, (int*)_data, fieldSize*sizeof(int));
}

void MhdFile::getDataAsUInt(uint* numpy_output_data, long field_size) {
	if(_dim==NULL)
		return;
	if(_data==NULL)
		return;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(numpy_output_data, (uint*)_data, fieldSize*sizeof(uint));
}

void MhdFile::getDataAsLong(long* numpy_output_data, long field_size) {
	if(_dim==NULL)
		return;
	if(_data==NULL)
		return;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(numpy_output_data, (long*)_data, fieldSize*sizeof(long));
}

void MhdFile::getDataAsULong(ulong* numpy_output_data, long field_size) {
	if(_dim==NULL)
		return;
	if(_data==NULL)
		return;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(numpy_output_data, (ulong*)_data, fieldSize*sizeof(ulong));
}

void MhdFile::getDataAsFloat(float* numpy_output_data, long field_size) {
	if(_dim==NULL)
		return;
	if(_data==NULL)
		return;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(numpy_output_data, (float*)_data, fieldSize*sizeof(float));
}

void MhdFile::getDataAsDouble(double* numpy_output_data, long field_size) {
	if(_dim==NULL)
		return;
	if(_data==NULL)
		return;
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	long fieldSize = (field_size<=fsize) ? field_size : fsize;
	memcpy(numpy_output_data, (double*)_data, fieldSize*sizeof(double));
}

void MhdFile::GetDimensions(unsigned int* numpy_output_dimensions, int field_size) {
	if(_dim==NULL)
		return;
	int fieldSize = (field_size<=_ndims) ? field_size : _ndims;
	memcpy(numpy_output_dimensions, _dim, fieldSize*sizeof(unsigned int));
}

void MhdFile::GetSpacing(float* numpy_output_spacing, int field_size) {
	if(_spacing==NULL)
		return;
	int fieldSize = (field_size<=_ndims) ? field_size : _ndims;
	memcpy(numpy_output_spacing, _spacing, sizeof(float)*fieldSize);
}

void MhdFile::readFile(void) {
	if(_filename.empty()) {
#ifdef DEBUG
		std::cerr << "No file information given." << std::endl;
#endif
	}
	int slashpos = _filename.find_last_of("/");
	int dotpos = _filename.find_last_of(".");
	std::string inputExtension = _filename.substr(dotpos+1, _filename.npos);
	if(inputExtension.find("mhd")==inputExtension.npos && inputExtension.find("raw")==inputExtension.npos) {
#ifdef DEBUG
		std::cout << "no volume .mhd file as input - EXITING ..." << std::endl;
#endif
		return;
	}
	std::string basepath = _filename.substr(0,slashpos);
	std::string basenameWend = _filename.substr(slashpos+1, _filename.npos);
	dotpos = basenameWend.find_last_of(".");
	//std::cout << basenameWend << std::endl;
	std::string basename = basenameWend.substr(0,dotpos);
	std::string mhdFilename = basepath+"/"+basename+"."+"mhd";
	std::string rawFilename = basepath+"/"+basename+"."+"raw";
#ifdef DEBUG
	std::cout << "Basepath: " << basepath << std::endl;
	std::cout << "Basename (with ending): " << basenameWend << std::endl;
	std::cout << "Basename (no ending): " << basename << std::endl;
	std::cout << "mhd filename: " << mhdFilename << std::endl;
	std::cout << "raw filename: " << rawFilename << std::endl;
#endif

	std::ifstream mhdFile;
	mhdFile.open(mhdFilename.c_str(), std::ios::in);
	std::vector<std::string> strs;
	strs.reserve(10);
	std::string line;
	if(mhdFile.is_open())
	{
		std::getline(mhdFile,line);
		while(mhdFile.good()==true)
		{
			strs.clear();
			if(line.empty()) {
				if(mhdFile.eof()==true)
					break;
				continue;
			}

			std::string trimmed = boost::algorithm::trim_copy(line);
			boost::split(strs, trimmed, boost::is_any_of(" "));
			for(uint k = 0; k < strs.size(); k++)
				strs.at(k) = boost::algorithm::trim_copy(strs.at(k));

			if(strs.at(0).find("CenterOfRotation")!=strs.at(0).npos) {
				// always 3D - could also be 2D ?
				uint ncenter = strs.size()-2;
				if(_volumeCenter!=NULL)
					delete [] _volumeCenter;
				_volumeCenter = new float[std::min(ncenter, uint(3))];
				for(uint i=0; i<(std::min(ncenter, uint(3))); i++) {
					_volumeCenter[i] = float(atof(strs.at(i+2).c_str()));
				}
			} else if(strs.at(0).find("ElementSpacing")!=strs.at(0).npos) {
				// always 3D - could also be 2D ?
				uint nspacing = strs.size()-2;
				if(_spacing!=NULL)
					delete [] _spacing;
				_spacing = new float[nspacing];
				for(uint i=0; i<nspacing; i++) {
					_spacing[i] = float(atof(strs.at(i+2).c_str()));
				}
			} else if(strs.at(0).find("DimSize")!=strs.at(0).npos) {
				uint ndims = strs.size()-2;
				if(_dim!=NULL)
					delete [] _dim;
				_dim = new uint[ndims];
				for(uint i=0; i<ndims; i++) {
					_dim[i] = uint(atoi(strs.at(i+2).c_str()));
				}
				_ndims = ndims;
			} else if(strs.at(0).find("ElementType")!=strs.at(0).npos) {
				std::string indicator = strs.at(2);
				if(indicator.find("_MET_CHAR")!=indicator.npos) {
					_datatype = CHAR;
				} else if(indicator.find("MET_UCHAR")!=indicator.npos) {
					_datatype = UCHAR;
				} else if(indicator.find("MET_SHORT")!=indicator.npos) {
					_datatype = SHORT;
				} else if(indicator.find("MET_USHORT")!=indicator.npos) {
					_datatype = USHORT;
				} else if(indicator.find("MET_INT")!=indicator.npos) {
					_datatype = INT;
				} else if(indicator.find("MET_UINT")!=indicator.npos) {
					_datatype = UINT;
				} else if(indicator.find("MET_LONG")!=indicator.npos) {
					_datatype = LONG;
				} else if(indicator.find("MET_ULONG")!=indicator.npos) {
					_datatype = ULONG;
				} else if(indicator.find("MET_FLOAT")!=indicator.npos) {
					_datatype = FLOAT;
				} else if(indicator.find("MET_DOUBLE")!=indicator.npos) {
					_datatype = DOUBLE;
				} else {
					_datatype = NONE;
				}
			} else if(strs.at(0).find("ElementDataFile")!=strs.at(0).npos) {
				std::string tempName = strs.at(2);
				if(rawFilename.find(tempName.c_str())==rawFilename.npos) {
					rawFilename = tempName;
				}
			}
			std::getline(mhdFile,line);
			if(mhdFile.eof()==true)
				break;
		}
#ifdef DEBUG
		std::cout << "File reading operation done." << std::endl;
#endif
		mhdFile.close();

	}
	else
	{
#ifdef DEBUG
		std::cout << "Unable to open mhd file." << std::endl;
#endif
		return;
	}

	clearData();
	/*
	long fsize = 1;
	for(int i=0; i<_ndims; i++)
		fsize*=long(_dim[i]);
	switch(_datatype) {
		case CHAR: {
			_data = new char[fsize];
			break;
		}
		case UCHAR: {
			_data = new uchar[fsize];
			break;
		}
		case SHORT: {
			_data = new short[fsize];
			break;
		}
		case USHORT: {
			_data = new ushort[fsize];
			break;
		}
		case INT: {
			_data = new int[fsize];
			break;
		}
		case UINT: {
			_data = new uint[fsize];
			break;
		}
		case LONG: {
			_data = new long[fsize];
			break;
		}
		case ULONG: {
			_data = new ulong[fsize];
			break;
		}
		case FLOAT: {
			_data = new float[fsize];
			break;
		}
		case DOUBLE: {
			_data = new double[fsize];
			break;
		}
		default: {
			break;
		}
	}
	*/

	std::ifstream rawFile;
	rawFile.open(rawFilename.c_str(), std::ios::in|std::ios::binary);
	if(rawFile.is_open()) {
		long fsize = 1;
		for(int i=0; i<_ndims; i++)
			fsize*=long(_dim[i]);

		switch(_datatype) {
			case CHAR: {
				_data = new char[fsize];
				char dataRead=0;
				for(long i=0;i<fsize;i++) {
					rawFile.read(reinterpret_cast <char*> (&dataRead), sizeof(char));
					((char*)_data)[i] = dataRead;
				}
				break;
			}
			case UCHAR: {
				_data = new uchar[fsize];
				uchar dataRead=0;
				for(long i=0;i<fsize;i++) {
					rawFile.read(reinterpret_cast <char*> (&dataRead), sizeof(uchar));
					((uchar*)_data)[i] = dataRead;
				}
				break;
			}
			case SHORT: {
				_data = new short[fsize];
				short dataRead=0;
				for(long i=0;i<fsize;i++) {
					rawFile.read(reinterpret_cast <char*> (&dataRead), sizeof(short));
					((short*)_data)[i] = dataRead;
				}
				break;
			}
			case USHORT: {
				_data = new ushort[fsize];
				ushort dataRead=0;
				for(long i=0;i<fsize;i++) {
					rawFile.read(reinterpret_cast <char*> (&dataRead), sizeof(ushort));
					((ushort*)_data)[i] = dataRead;
				}
				break;
			}
			case INT: {
				_data = new int[fsize];
				int dataRead=0;
				for(long i=0;i<fsize;i++) {
					rawFile.read(reinterpret_cast <char*> (&dataRead), sizeof(int));
					((int*)_data)[i] = dataRead;
				}
				break;
			}
			case UINT: {
				_data = new uint[fsize];
				uint dataRead=0;
				for(long i=0;i<fsize;i++) {
					rawFile.read(reinterpret_cast <char*> (&dataRead), sizeof(uint));
					((uint*)_data)[i] = dataRead;
				}
				break;
			}
			case LONG: {
				_data = new long[fsize];
				long dataRead=0;
				for(long i=0;i<fsize;i++) {
					rawFile.read(reinterpret_cast <char*> (&dataRead), sizeof(long));
					((long*)_data)[i] = dataRead;
				}
				break;
			}
			case ULONG: {
				_data = new ulong[fsize];
				ulong dataRead=0;
				for(long i=0;i<fsize;i++) {
					rawFile.read(reinterpret_cast <char*> (&dataRead), sizeof(ulong));
					((ulong*)_data)[i] = dataRead;
				}
				break;
			}
			case FLOAT: {
				_data = new float[fsize];
				float dataRead=0;
				for(long i=0;i<fsize;i++) {
					rawFile.read(reinterpret_cast <char*> (&dataRead), sizeof(float));
					if(std::isinf(dataRead))
						dataRead = 100.0;
					if(std::isnan(dataRead))
						dataRead = 0.0;
					((float*)_data)[i] = dataRead;
				}
				break;
			}
			case DOUBLE: {
				_data = new double[fsize];
				double dataRead=0;
				for(long i=0;i<fsize;i++) {
					rawFile.read(reinterpret_cast <char*> (&dataRead), sizeof(double));
					if(std::isinf(dataRead))
						dataRead = 100.0;
					if(std::isnan(dataRead))
						dataRead = 0.0;
					((double*)_data)[i] = dataRead;
				}
				break;
			}
			default: {
				break;
			}
		}
		rawFile.close();
	} else {
#ifdef DEBUG
		std::cout << "Cannot open raw file for writing" << std::endl;
#endif
	}
}

void MhdFile::writeFile(void) {
	if(_dim==NULL) {
#ifdef DEBUG
		std::cerr << "No dimensionality information provided." << std::endl;
#endif
		return;
	}
	if(_spacing==NULL) {
#ifdef DEBUG
		std::cerr << "No spacing information provided." << std::endl;
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
	if(inputExtension.find("mhd")==inputExtension.npos && inputExtension.find("raw")==inputExtension.npos) {
		std::cout << "no volume .mhd file - EXITING ..." << std::endl;
		return;
	}
	std::string basepath = _filename.substr(0,slashpos);
	std::string basenameWend = _filename.substr(slashpos+1, _filename.npos);
	dotpos = basenameWend.find_last_of(".");
	//std::cout << basenameWend << std::endl;
	std::string basename = basenameWend.substr(0,dotpos);
	std::string mhdFilename = basepath+"/"+basename+"."+"mhd";
	std::string rawFilename = basepath+"/"+basename+"."+"raw";
	//std::cout << basepath+"/"+basename << std::endl;

	if(_volumeCenter==NULL)
		_volumeCenter = new float[3];
	bzero(_volumeCenter, 3*sizeof(float));

	std::ofstream mhdFile;
	mhdFile.open(mhdFilename.c_str(), std::ios::out|std::ios::trunc);
	if(mhdFile.is_open())
	{
		mhdFile << "ObjectType = Image" << std::endl;
		mhdFile << "NDims = " << _ndims << std::endl;
		mhdFile << "BinaryData = True"<< std::endl;
		mhdFile << "BinaryDataByteOrderMSB = False"<< std::endl;
		mhdFile << "CompressedData = False"<< std::endl;
		if(_ndims>2) {
			mhdFile << "TransformMatrix = 1 0 0 0 1 0 0 0 1"<< std::endl;
			mhdFile << "Offset = 0 0 0"<< std::endl;
			mhdFile << "CenterOfRotation = " << _volumeCenter[0] << " " << _volumeCenter[1] << " " << _volumeCenter[2] << std::endl;
		}
		mhdFile << "ElementSpacing = ";
		for(int i=0; i<_ndims; i++) {
			mhdFile <<_spacing[i];
			if(i<(_ndims-1))
				mhdFile << " ";
		}
		mhdFile << std::endl;
		mhdFile << "DimSize = ";
		for(int i=0; i<_ndims; i++) {
			mhdFile << _dim[i];
			if(i<(_ndims-1))
				mhdFile << " ";
		}
		mhdFile << std::endl;
		if(_ndims>2) {
			mhdFile << "AnatomicalOrientation = ???"<< std::endl;
		}
		switch(_datatype) {
			case CHAR: {
				mhdFile << "ElementType = " << "MET_CHAR" << std::endl;
				break;
			}
			case UCHAR: {
				mhdFile << "ElementType = " << "MET_UCHAR" << std::endl;
				break;
			}
			case SHORT: {
				mhdFile << "ElementType = " << "MET_SHORT" << std::endl;
				break;
			}
			case USHORT: {
				mhdFile << "ElementType = " << "MET_USHORT" << std::endl;
				break;
			}
			case INT: {
				mhdFile << "ElementType = " << "MET_INT" << std::endl;
				break;
			}
			case UINT: {
				mhdFile << "ElementType = " << "MET_UINT" << std::endl;
				break;
			}
			case LONG: {
				mhdFile << "ElementType = " << "MET_LONG" << std::endl;
				break;
			}
			case ULONG: {
				mhdFile << "ElementType = " << "MET_ULONG" << std::endl;
				break;
			}
			case FLOAT: {
				mhdFile << "ElementType = " << "MET_FLOAT" << std::endl;
				break;
			}
			case DOUBLE: {
				mhdFile << "ElementType = " << "MET_DOUBLE" << std::endl;
				break;
			}
			default: {
				break;
			}
		}
		mhdFile << "ElementDataFile = " << basename+"."+"raw" << std::endl;
		mhdFile.close();
	}
	else {
		std::cout << "Cannot open mhd file for writing." << std::endl;
	}

	std::ofstream rawFile;
	rawFile.open(rawFilename.c_str(), std::ios::out|std::ios::binary);
	if(rawFile.is_open()) {
		long fsize = 1;
		for(int i=0; i<_ndims; i++)
			fsize*=long(_dim[i]);
		//fixed 12-bit integer precision - with resize

		for(long i=0; i<fsize; i++) {
			switch(_datatype) {
				case CHAR: {
					char* mdata = (char*)_data;
					rawFile.write(reinterpret_cast <char*> (&mdata[i]), sizeof(char));
					break;
				}
				case UCHAR: {
					uchar* mdata = (uchar*)_data;
					rawFile.write(reinterpret_cast <char*> (&mdata[i]), sizeof(uchar));
					break;
				}
				case SHORT: {
					short* mdata = (short*)_data;
					rawFile.write(reinterpret_cast <char*> (&mdata[i]), sizeof(short));
					break;
				}
				case USHORT: {
					ushort* mdata = (ushort*)_data;
					rawFile.write(reinterpret_cast <char*> (&mdata[i]), sizeof(ushort));
					break;
				}
				case INT: {
					int* mdata = (int*)_data;
					rawFile.write(reinterpret_cast <char*> (&mdata[i]), sizeof(int));
					break;
				}
				case UINT: {
					uint* mdata = (uint*)_data;
					rawFile.write(reinterpret_cast <char*> (&mdata[i]), sizeof(uint));
					break;
				}
				case LONG: {
					long* mdata = (long*)_data;
					rawFile.write(reinterpret_cast <char*> (&mdata[i]), sizeof(long));
					break;
				}
				case ULONG: {
					ulong* mdata = (ulong*)_data;
					rawFile.write(reinterpret_cast <char*> (&mdata[i]), sizeof(ulong));
					break;
				}
				case FLOAT: {
					float* mdata = (float*)_data;
					rawFile.write(reinterpret_cast <char*> (&mdata[i]), sizeof(float));
					break;
				}
				case DOUBLE: {
					double* mdata = (double*)_data;
					rawFile.write(reinterpret_cast <char*> (&mdata[i]), sizeof(double));
					break;
				}
				default: {
					break;
				}
			}
		}

		rawFile.close();
	} else {
		std::cout << "Cannot open dat file for writing" << std::endl;
	}
}



