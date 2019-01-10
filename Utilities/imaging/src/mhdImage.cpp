/**
 * Created on 14 Jan 2018
 * author: christian
 */
#include "mhdImage.hpp"
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

template<typename T>
MhdImage2D<T>::MhdImage2D() : _data(NULL), _spacing(NULL), _filename(""), _externalData(false) {

}

template<typename T>
MhdImage2D<T>::~MhdImage2D() {
	if(_data!=NULL) {
		if(!_externalData) {
			_data->clean();
			delete _data;
		}
		_data=NULL;
	}
	if(_spacing!=NULL) {
		delete [] _spacing;
	}
}

template<typename T>
void MhdImage2D<T>::read() {
	if(_data==NULL) {
		std::cerr << "empty image in MhdImage2D reader" << std::endl;
		exit(1);
	}

	if(_filename.empty()) {
		std::cerr << "No filename provided to MhdImage2D reader." << std::endl;
		exit(1);
	}
	const int _ndims=2;

	int slashpos = _filename.find_last_of("/");
	int dotpos = _filename.find_last_of(".");
	std::string inputExtension = _filename.substr(dotpos+1, _filename.npos);
	if(inputExtension.find("mhd")==inputExtension.npos && inputExtension.find("raw")==inputExtension.npos) {
		std::cerr << "no volume .mhd file as input - EXITING ..." << std::endl;
		return;
	}
	std::string basepath = _filename.substr(0,slashpos);
	std::string basenameWend = _filename.substr(slashpos+1, _filename.npos);
	dotpos = basenameWend.find_last_of(".");
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
				// skip - that's a 2D reader giving out 2D images
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
				uint dim[_ndims];
				for(uint i=0; i<ndims; i++) {
					dim[i] = uint(atoi(strs.at(i+2).c_str()));
				}
				_data->setDimensions(dim, true);
			} else if(strs.at(0).find("ElementType")!=strs.at(0).npos) {
				std::string indicator = strs.at(2);
				dtype typeRead;
				if(indicator.find("_MET_CHAR")!=indicator.npos) {
					typeRead = CHAR;
				} else if(indicator.find("MET_UCHAR")!=indicator.npos) {
					typeRead = UCHAR;
				} else if(indicator.find("MET_SHORT")!=indicator.npos) {
					typeRead = SHORT;
				} else if(indicator.find("MET_USHORT")!=indicator.npos) {
					typeRead = USHORT;
				} else if(indicator.find("MET_INT")!=indicator.npos) {
					typeRead = INT;
				} else if(indicator.find("MET_UINT")!=indicator.npos) {
					typeRead = UINT;
				} else if(indicator.find("MET_LONG")!=indicator.npos) {
					typeRead = LONG;
				} else if(indicator.find("MET_ULONG")!=indicator.npos) {
					typeRead = ULONG;
				} else if(indicator.find("MET_FLOAT")!=indicator.npos) {
					typeRead = FLOAT;
				} else if(indicator.find("MET_DOUBLE")!=indicator.npos) {
					typeRead = DOUBLE;
				} else {
					typeRead = NONE;
				}

				if((_data->getDatatype()!=NONE) && (_data->getDatatype()!=typeRead)) {
					std::cerr << "Expected datatype and provided datatype do not match - please provide data in " << getDTypeString(_data->getDatatype()) << " format." << std::endl;
					return;
				}
				if(_data->getDatatype()==NONE)
					_data->setDatatype(typeRead);
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

	// now pre-defined in image: dataorder; read: datatype, dimensions; All ready to allocate image
	if(_data->getData()==NULL)
		_data->createImage();
	else {
		_data->clean();
	}

	std::ifstream rawFile;
	rawFile.open(rawFilename.c_str(), std::ios::in|std::ios::binary);
	if(rawFile.is_open()) {
		long fsize = 1;
		for(int i=0; i<_ndims; i++)
			fsize*=long(_data->getDimensions()[i]);

		T dataRead;
		for(long i=0;i<fsize;i++) {
			rawFile.read(reinterpret_cast<char*>(&dataRead), sizeof(T));
			_data->getData()[i]=dataRead;
		}

		rawFile.close();
	} else {
#ifdef DEBUG
		std::cout << "Cannot open raw file for writing" << std::endl;
#endif
	}
}

template<typename T>
void MhdImage2D<T>::write() {
	if(_data==NULL) {
		std::cerr << "No data provided." << std::endl;
		exit(1);
	}
	if(_spacing==NULL) {
		std::cerr << "No spacing information provided." << std::endl;
		exit(1);
	}
	if(_filename.empty()) {
		std::cerr << "empty filename." << std::endl;
		exit(1);
	}
	// check if image is valid
	if(isValid(_data)==false) {
		std::cerr << "Invalid image data to write." << std::endl;
		exit(1);
	}

	const int _ndims=2;

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
	std::string basename = basenameWend.substr(0,dotpos);
	std::string mhdFilename = basepath+"/"+basename+"."+"mhd";
	std::string rawFilename = basepath+"/"+basename+"."+"raw";

	//if(_volumeCenter==NULL)
	//	_volumeCenter = new float[3];
	//bzero(_volumeCenter, 3*sizeof(float));

	std::ofstream mhdFile;
	mhdFile.open(mhdFilename.c_str(), std::ios::out|std::ios::trunc);
	if(mhdFile.is_open())
	{
		mhdFile << "ObjectType = Image" << std::endl;
		mhdFile << "NDims = " << _ndims << std::endl;
		mhdFile << "BinaryData = True"<< std::endl;
		mhdFile << "BinaryDataByteOrderMSB = False"<< std::endl;
		mhdFile << "CompressedData = False"<< std::endl;
		//if(_ndims>2) {
		//	mhdFile << "TransformMatrix = 1 0 0 0 1 0 0 0 1"<< std::endl;
		//	mhdFile << "Offset = 0 0 0"<< std::endl;
		//	mhdFile << "CenterOfRotation = " << _volumeCenter[0] << " " << _volumeCenter[1] << " " << _volumeCenter[2] << std::endl;
		//}
		mhdFile << "ElementSpacing = ";
		for(int i=0; i<_ndims; i++) {
			mhdFile <<_spacing[i];
			if(i<(_ndims-1))
				mhdFile << " ";
		}
		mhdFile << std::endl;
		mhdFile << "DimSize = ";
		for(int i=0; i<_ndims; i++) {
			mhdFile << _data->getDimension(i);
			if(i<(_ndims-1))
				mhdFile << " ";
		}
		mhdFile << std::endl;
		if(_ndims>2) {
			mhdFile << "AnatomicalOrientation = ???"<< std::endl;
		}
		switch(_data->getDatatype()) {
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
			fsize*=long(_data->getDimensions()[i]);
		T* mdata = _data->getData();
		for(long i=0; i<fsize; i++) {
			rawFile.write(reinterpret_cast<char*>(&mdata[i]), sizeof(T));
		}
		rawFile.close();
	} else {
		std::cout << "Cannot open dat file for writing" << std::endl;
	}
}

template class MhdImage2D<char>;
template class MhdImage2D<unsigned char>;
template class MhdImage2D<short>;
template class MhdImage2D<unsigned short>;
template class MhdImage2D<int>;
template class MhdImage2D<unsigned int>;
template class MhdImage2D<long>;
template class MhdImage2D<unsigned long>;
template class MhdImage2D<float>;
template class MhdImage2D<double>;





template<typename T>
MhdImage3D<T>::MhdImage3D() : _data(NULL), _spacing(NULL), _volumeCenter(NULL), _filename(""), _externalData(false) {

}

template<typename T>
MhdImage3D<T>::~MhdImage3D() {
	if(_data!=NULL) {
		if(!_externalData) {
			_data->clean();
			delete _data;
		}
		_data=NULL;
	}
	if(_spacing!=NULL) {
		delete [] _spacing;
	}
	if(_volumeCenter!=NULL) {
		delete [] _volumeCenter;
	}
}

template<typename T>
void MhdImage3D<T>::read() {
	if(_data==NULL) {
		std::cerr << "empty image in MhdImage2D reader" << std::endl;
		exit(1);
	}

	if(_filename.empty()) {
		std::cerr << "No filename provided to MhdImage2D reader." << std::endl;
		exit(1);
	}
	const int _ndims=3;

	int slashpos = _filename.find_last_of("/");
	int dotpos = _filename.find_last_of(".");
	std::string inputExtension = _filename.substr(dotpos+1, _filename.npos);
	if(inputExtension.find("mhd")==inputExtension.npos && inputExtension.find("raw")==inputExtension.npos) {
		std::cerr << "no .mhd file as input - EXITING ..." << std::endl;
		return;
	}
	std::string basepath = _filename.substr(0,slashpos);
	std::string basenameWend = _filename.substr(slashpos+1, _filename.npos);
	dotpos = basenameWend.find_last_of(".");
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
				assert(ndims == _ndims);
				uint dim[_ndims];
				for(uint i=0; i<ndims; i++) {
					dim[i] = uint(atoi(strs.at(i+2).c_str()));
				}
				_data->setDimensions(dim, true);
			} else if(strs.at(0).find("ElementType")!=strs.at(0).npos) {
				std::string indicator = strs.at(2);
				dtype typeRead;
				if(indicator.find("_MET_CHAR")!=indicator.npos) {
					typeRead = CHAR;
				} else if(indicator.find("MET_UCHAR")!=indicator.npos) {
					typeRead = UCHAR;
				} else if(indicator.find("MET_SHORT")!=indicator.npos) {
					typeRead = SHORT;
				} else if(indicator.find("MET_USHORT")!=indicator.npos) {
					typeRead = USHORT;
				} else if(indicator.find("MET_INT")!=indicator.npos) {
					typeRead = INT;
				} else if(indicator.find("MET_UINT")!=indicator.npos) {
					typeRead = UINT;
				} else if(indicator.find("MET_LONG")!=indicator.npos) {
					typeRead = LONG;
				} else if(indicator.find("MET_ULONG")!=indicator.npos) {
					typeRead = ULONG;
				} else if(indicator.find("MET_FLOAT")!=indicator.npos) {
					typeRead = FLOAT;
				} else if(indicator.find("MET_DOUBLE")!=indicator.npos) {
					typeRead = DOUBLE;
				} else {
					typeRead = NONE;
				}

				if((_data->getDatatype()!=NONE) && (_data->getDatatype()!=typeRead)) {
					std::cerr << "Expected datatype and provided datatype do not match - please provide data in " << getDTypeString(_data->getDatatype()) << " format." << std::endl;
					return;
				}
				if(_data->getDatatype()==NONE)
					_data->setDatatype(typeRead);
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

	// now pre-defined in image: dataorder; read: datatype, dimensions; All ready to allocate image
	if(_data->getData()==NULL)
		_data->createImage();
	else {
		_data->clean();
	}

	std::ifstream rawFile;
	rawFile.open(rawFilename.c_str(), std::ios::in|std::ios::binary);
	if(rawFile.is_open()) {
		long fsize = 1;
		for(int i=0; i<_ndims; i++)
			fsize*=long(_data->getDimensions()[i]);

		T dataRead;
		for(long i=0;i<fsize;i++) {
			rawFile.read(reinterpret_cast<char*>(&dataRead), sizeof(T));
			_data->getData()[i]=dataRead;
		}

		rawFile.close();
	} else {
#ifdef DEBUG
		std::cout << "Cannot open raw file for writing" << std::endl;
#endif
	}
}

template<typename T>
void MhdImage3D<T>::write() {
	if(_data==NULL) {
		std::cerr << "No data provided." << std::endl;
		exit(1);
	}
	if(_spacing==NULL) {
		std::cerr << "No spacing information provided." << std::endl;
		exit(1);
	}
	if(_filename.empty()) {
		std::cerr << "empty filename." << std::endl;
		exit(1);
	}
	// check if image is valid
	if(isValid(_data)==false) {
		std::cerr << "Invalid image data to write." << std::endl;
		exit(1);
	}

	const int _ndims=3;

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
	std::string basename = basenameWend.substr(0,dotpos);
	std::string mhdFilename = basepath+"/"+basename+"."+"mhd";
	std::string rawFilename = basepath+"/"+basename+"."+"raw";

	if(_volumeCenter==NULL) {
		_volumeCenter = new float[3];
		//bzero(_volumeCenter, 3*sizeof(float));
		memset(_volumeCenter, 0, 3*sizeof(float));
	}

	std::ofstream mhdFile;
	mhdFile.open(mhdFilename.c_str(), std::ios::out|std::ios::trunc);
	if(mhdFile.is_open())
	{
		mhdFile << "ObjectType = Image" << std::endl;
		mhdFile << "NDims = " << _ndims << std::endl;
		mhdFile << "BinaryData = True"<< std::endl;
		mhdFile << "BinaryDataByteOrderMSB = False"<< std::endl;
		mhdFile << "CompressedData = False"<< std::endl;
		mhdFile << "TransformMatrix = 1 0 0 0 1 0 0 0 1"<< std::endl;
		mhdFile << "Offset = 0 0 0"<< std::endl;
		mhdFile << "CenterOfRotation = " << _volumeCenter[0] << " " << _volumeCenter[1] << " " << _volumeCenter[2] << std::endl;
		mhdFile << "ElementSpacing = ";
		for(int i=0; i<_ndims; i++) {
			mhdFile <<_spacing[i];
			if(i<(_ndims-1))
				mhdFile << " ";
		}
		mhdFile << std::endl;
		mhdFile << "DimSize = ";
		for(int i=0; i<_ndims; i++) {
			mhdFile << _data->getDimension(i);
			if(i<(_ndims-1))
				mhdFile << " ";
		}
		mhdFile << std::endl;
		if(_ndims>2) {
			mhdFile << "AnatomicalOrientation = ???"<< std::endl;
		}
		switch(_data->getDatatype()) {
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
			fsize*=long(_data->getDimensions()[i]);
		T* mdata = _data->getData();
		for(long i=0; i<fsize; i++) {
			rawFile.write(reinterpret_cast<char*>(&mdata[i]), sizeof(T));
		}
		rawFile.close();
	} else {
		std::cout << "Cannot open dat file for writing" << std::endl;
	}
}

template class MhdImage3D<char>;
template class MhdImage3D<unsigned char>;
template class MhdImage3D<short>;
template class MhdImage3D<unsigned short>;
template class MhdImage3D<int>;
template class MhdImage3D<unsigned int>;
template class MhdImage3D<long>;
template class MhdImage3D<unsigned long>;
template class MhdImage3D<float>;
template class MhdImage3D<double>;





template<typename T>
MhdImage4D<T>::MhdImage4D() : _data(NULL), _spacing(NULL), _volumeCenter(NULL), _filename(""), _externalData(false) {

}

template<typename T>
MhdImage4D<T>::~MhdImage4D() {
	if(_data!=NULL) {
		if(!_externalData) {
			_data->clean();
			delete _data;
		}
		_data=NULL;
	}
	if(_spacing!=NULL) {
		delete [] _spacing;
	}
	if(_volumeCenter!=NULL) {
		delete [] _volumeCenter;
	}
}

template<typename T>
void MhdImage4D<T>::read() {
	if(_data==NULL) {
		std::cerr << "empty image in MhdImage2D reader" << std::endl;
		exit(1);
	}

	if(_filename.empty()) {
		std::cerr << "No filename provided to MhdImage2D reader." << std::endl;
		exit(1);
	}
	const int _ndims=4;

	int slashpos = _filename.find_last_of("/");
	int dotpos = _filename.find_last_of(".");
	std::string inputExtension = _filename.substr(dotpos+1, _filename.npos);
	if(inputExtension.find("mhd")==inputExtension.npos && inputExtension.find("raw")==inputExtension.npos) {
		std::cerr << "no volume .mhd file as input - EXITING ..." << std::endl;
		return;
	}
	std::string basepath = _filename.substr(0,slashpos);
	std::string basenameWend = _filename.substr(slashpos+1, _filename.npos);
	dotpos = basenameWend.find_last_of(".");
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
				assert(ndims == _ndims);
				uint dim[_ndims];
				for(uint i=0; i<ndims; i++) {
					dim[i] = uint(atoi(strs.at(i+2).c_str()));
				}
				_data->setDimensions(dim, true);
			} else if(strs.at(0).find("ElementType")!=strs.at(0).npos) {
				std::string indicator = strs.at(2);
				dtype typeRead;
				if(indicator.find("_MET_CHAR")!=indicator.npos) {
					typeRead = CHAR;
				} else if(indicator.find("MET_UCHAR")!=indicator.npos) {
					typeRead = UCHAR;
				} else if(indicator.find("MET_SHORT")!=indicator.npos) {
					typeRead = SHORT;
				} else if(indicator.find("MET_USHORT")!=indicator.npos) {
					typeRead = USHORT;
				} else if(indicator.find("MET_INT")!=indicator.npos) {
					typeRead = INT;
				} else if(indicator.find("MET_UINT")!=indicator.npos) {
					typeRead = UINT;
				} else if(indicator.find("MET_LONG")!=indicator.npos) {
					typeRead = LONG;
				} else if(indicator.find("MET_ULONG")!=indicator.npos) {
					typeRead = ULONG;
				} else if(indicator.find("MET_FLOAT")!=indicator.npos) {
					typeRead = FLOAT;
				} else if(indicator.find("MET_DOUBLE")!=indicator.npos) {
					typeRead = DOUBLE;
				} else {
					typeRead = NONE;
				}

				if((_data->getDatatype()!=NONE) && (_data->getDatatype()!=typeRead)) {
					std::cerr << "Expected datatype and provided datatype do not match - please provide data in " << getDTypeString(_data->getDatatype()) << " format." << std::endl;
					exit(1);
				}
				if(_data->getDatatype()==NONE)
					_data->setDatatype(typeRead);
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

	// now pre-defined in image: dataorder; read: datatype, dimensions; All ready to allocate image
	if(_data->getData()==NULL)
		_data->createImage();
	else {
		_data->clean();
	}

	std::ifstream rawFile;
	rawFile.open(rawFilename.c_str(), std::ios::in|std::ios::binary);
	if(rawFile.is_open()) {
		long fsize = 1;
		for(int i=0; i<_ndims; i++)
			fsize*=long(_data->getDimensions()[i]);

		T dataRead;
		for(long i=0;i<fsize;i++) {
			rawFile.read(reinterpret_cast<char*>(&dataRead), sizeof(T));
			_data->getData()[i]=dataRead;
		}

		rawFile.close();
	} else {
#ifdef DEBUG
		std::cout << "Cannot open raw file for writing" << std::endl;
#endif
	}
}

template<typename T>
void MhdImage4D<T>::write() {
	if(_data==NULL) {
		std::cerr << "No data provided." << std::endl;
		exit(1);
	}
	if(_spacing==NULL) {
		std::cerr << "No spacing information provided." << std::endl;
		exit(1);
	}
	if(_filename.empty()) {
		std::cerr << "empty filename." << std::endl;
		exit(1);
	}
	// check if image is valid
	if(isValid(_data)==false) {
		std::cerr << "Invalid image data to write." << std::endl;
		exit(1);
	}

	const int _ndims=4;

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
	std::string basename = basenameWend.substr(0,dotpos);
	std::string mhdFilename = basepath+"/"+basename+"."+"mhd";
	std::string rawFilename = basepath+"/"+basename+"."+"raw";

	if(_volumeCenter==NULL) {
		_volumeCenter = new float[3];
		//bzero(_volumeCenter, 3*sizeof(float));
		memset(_volumeCenter, 0, 3*sizeof(float));
	}

	std::ofstream mhdFile;
	mhdFile.open(mhdFilename.c_str(), std::ios::out|std::ios::trunc);
	if(mhdFile.is_open())
	{
		mhdFile << "ObjectType = Image" << std::endl;
		mhdFile << "NDims = " << _ndims << std::endl;
		mhdFile << "BinaryData = True"<< std::endl;
		mhdFile << "BinaryDataByteOrderMSB = False"<< std::endl;
		mhdFile << "CompressedData = False"<< std::endl;
		mhdFile << "TransformMatrix = 1 0 0 0 1 0 0 0 1"<< std::endl;
		mhdFile << "Offset = 0 0 0"<< std::endl;
		mhdFile << "CenterOfRotation = " << _volumeCenter[0] << " " << _volumeCenter[1] << " " << _volumeCenter[2] << std::endl;
		mhdFile << "ElementSpacing = ";
		for(int i=0; i<_ndims; i++) {
			mhdFile <<_spacing[i];
			if(i<(_ndims-1))
				mhdFile << " ";
		}
		mhdFile << std::endl;
		mhdFile << "DimSize = ";
		for(int i=0; i<_ndims; i++) {
			mhdFile << _data->getDimension(i);
			if(i<(_ndims-1))
				mhdFile << " ";
		}
		mhdFile << std::endl;
		if(_ndims>2) {
			mhdFile << "AnatomicalOrientation = ???"<< std::endl;
		}
		switch(_data->getDatatype()) {
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
			fsize*=long(_data->getDimensions()[i]);
		T* mdata = _data->getData();
		for(long i=0; i<fsize; i++) {
			rawFile.write(reinterpret_cast<char*>(&mdata[i]), sizeof(T));
		}
		rawFile.close();
	} else {
		std::cout << "Cannot open dat file for writing" << std::endl;
	}
}

template class MhdImage4D<char>;
template class MhdImage4D<unsigned char>;
template class MhdImage4D<short>;
template class MhdImage4D<unsigned short>;
template class MhdImage4D<int>;
template class MhdImage4D<unsigned int>;
template class MhdImage4D<long>;
template class MhdImage4D<unsigned long>;
template class MhdImage4D<float>;
template class MhdImage4D<double>;





