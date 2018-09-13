/**
 * Created on 06 Apr 2018
 * author: christian
 */

#ifndef IMAGE4D_UTILS_HPP_
#define IMAGE4D_UTILS_HPP_

#include "image.hpp"
#ifndef _CXX_COMPILE_
#include "../../common/src/std_typedefs.h"
#else
#include "Utilities/common/src/std_typedefs.h"
#endif
#include <iostream>
#include <vector>
#include <string.h>
#include <cmath>
#include <omp.h>

typedef enum _NEIGHBOURHOOD4D {
	NBR_9s  = 0,
	NBR_81  = 1,
	NBR_81W = 2,
} NEIGHBOURHOOD4D;

// dissimilarity measure between pixels
template<typename Scalar>
inline Scalar SSD_imageStack(Image4D<Scalar>& data, int x1, int y1, int z1, int x2, int y2, int z2) {
	float result = 0;
	for(unsigned int c = 0; c < data.getDimension(3); c++) {
		result += (data.getImageValue(x1, y1, z1, c)-data.getImageValue(x2, y2, z2, c))*(data.getImageValue(x1, y1, z1, c)-data.getImageValue(x2, y2, z2, c));
	}
	return std::sqrt(result);
}

template<typename Scalar>
inline void createNonUniformAverageImage(Image4D<Scalar>& input, std::vector< std::pair<uint, uint> >& bounds, Image4D<Scalar>& output) {
	uint dim0=input.getDimension(0), dim1=input.getDimension(1), dim2=input.getDimension(2);
	for(uint boundaryIdx=0; boundaryIdx<bounds.size(); boundaryIdx++) {
		uint minChannel = bounds[boundaryIdx].first;
		uint maxChannel = bounds[boundaryIdx].second;
		#pragma omp parallel for shared(input,output,minChannel,maxChannel,boundaryIdx,dim0,dim1,dim2)
		for(ulong i=0; i<(dim0*dim1*dim2); i++) {
			uint x=0, y=0, z=0;
			getIndicesFromAddress_field3D(i,x,dim0,y,dim1,z,dim2,input.getDataorder());
			double avgVal = 0;
			for(uint c=minChannel; c<maxChannel; c++) {
				avgVal+=input.getImageValue(x,y,z,c);
			}
			avgVal /= double(maxChannel-minChannel);
			output.setImageValue(avgVal, x,y,z,boundaryIdx);
		}
	}
}

template<typename Scalar>
std::vector< std::vector<Scalar> > computeResponseCurves_uniformSample(Image4D<Scalar>& input, uint stride) {
	std::vector< std::vector<Scalar> > result;
	for(uint x=0; x<input.getDimension(0); x+=stride) {
		for(uint y=0; y<input.getDimension(1); y+=stride) {
			for(uint z=0; z<input.getDimension(2); z+=stride) {
				std::vector<Scalar> sample;
				for(uint channel=0; channel<input.getDimension(3); channel++) {
					sample.push_back(input.getImageValue(x,y,z,channel));
				}
				result.push_back(sample);
			}
		}
	}
	return result;
}

template<typename Scalar>
void computeLocalMean(Image4D<Scalar>& input, Image4D<Scalar>& result, uint kernelSize) {
	int K2 = int(std::floor(kernelSize/2));
	int K2P1 = K2+1;

	ulong numElem = result.getDimension(0)*result.getDimension(1)*result.getDimension(2)*result.getDimension(3);
	uint dim0 = result.getDimension(0); uint dim1 = result.getDimension(1); uint dim2 = result.getDimension(2); uint dim3 = result.getDimension(3);
	Scalar* aux_mem;
	aux_mem = new Scalar[numElem];

	#pragma omp parallel for shared(input,aux_mem)
	for(ulong eNum=0; eNum<numElem; eNum++) {
		uint x=0, y=0, z=0, c=0;
		getIndicesFromAddress_field4D(eNum,x,dim0,y,dim1,z,dim2,c,dim3,input.getDataorder());
		Scalar lmean = 0;
		for(int k=(int(x)-K2); k<(int(x)+K2P1); k++) {
			for(int l=(int(y)-K2); l<(int(y)+K2P1); l++) {
				for(int m=(int(z)-K2); m<(int(z)+K2P1); m++) {
					if((k<0)||(uint(k)>(dim0-1)) || (l<0)||(uint(l)>(dim1-1)) || (m<0)||(uint(m)>(dim2-1))) {
						lmean+=0;
					} else {
						lmean += input.GetImageValue(k,l,m,c);
					}
				}
			}
		}
		lmean /= (kernelSize*kernelSize*kernelSize);
		aux_mem[eNum] = lmean;
	}

	#pragma omp parallel for shared(aux_mem,result)
	for(ulong eNum=0; eNum<numElem; eNum++) {
		uint x=0, y=0, z=0, c=0;
		getIndicesFromAddress_field4D(eNum,x,dim0,y,dim1,z,dim2,c,dim3,result.getDataorder());
		result.setImageValue(aux_mem[eNum],x,y,z,c);
	}

	delete [] aux_mem;
}

template<typename Scalar>
void computeLocalStdDev(Image4D<Scalar>& input, Image4D<Scalar>& local_mean, Image4D<Scalar>& result, uint kernelSize) {
	int K2 = int(std::floor(kernelSize/2));
	int K2P1 = K2+1;
	ulong numElem = result.getDimension(0)*result.getDimension(1)*result.getDimension(2)*result.getDimension(3);
	uint dim0 = result.getDimension(0); uint dim1 = result.getDimension(1); uint dim2 = result.getDimension(2); uint dim3 = result.getDimension(3);
	Scalar* aux_mem;
	aux_mem = new Scalar[numElem];


	#pragma omp parallel for shared(aux_mem,local_mean,input)
	for(ulong eNum=0; eNum<numElem; eNum++) {
		uint x=0, y=0, z=0, c=0;
		getIndicesFromAddress_field4D(eNum,x,dim0,y,dim1,z,dim2,c,dim3,input.getDataorder());
		Scalar var = 0;
		for(int k=(int(x)-K2); k<(int(x)+K2P1); k++) {
			for(int l=(int(y)-K2); l<(int(y)+K2P1); l++) {
				for(int m=(int(z)-K2); m<(int(z)+K2P1); m++) {
					if((k<0)||(uint(k)>(dim0-1)) || (l<0)||(uint(l)>(dim1-1)) || (m<0)||(uint(m)>(dim2-1))) {
						var+=0;
					} else {
						var += pow(double(input.GetImageValue(k,l,m,c) - local_mean.GetImageValue(k,l,m,c)),2);
					}
				}
			}
		}
		var /= (kernelSize*kernelSize*kernelSize);
		aux_mem[eNum] = std::sqrt(var);
	}

	#pragma omp parallel for shared(aux_mem,result)
	for(ulong eNum=0; eNum<numElem; eNum++) {
		uint x=0, y=0, z=0, c=0;
		getIndicesFromAddress_field4D(eNum,x,dim0,y,dim1,z,dim2,c,dim3,result.getDataorder());
		result.setImageValue(aux_mem[eNum],x,y,z,c);
	}

	delete [] aux_mem;
}

template<typename Scalar>
void computeDeviation(Image4D<Scalar>& input, Image4D<Scalar>& local_mean, Image4D<Scalar>& result) {
	ulong numElem = result.getDimension(0)*result.getDimension(1)*result.getDimension(2)*result.getDimension(3);
	uint dim0 = result.getDimension(0); uint dim1 = result.getDimension(1); uint dim2 = result.getDimension(2); uint dim3 = result.getDimension(3);
	Scalar* aux_mem;
	aux_mem = new Scalar[numElem];


	#pragma omp parallel for shared(aux_mem,local_mean,input)
	for(ulong eNum=0; eNum<numElem; eNum++) {
		uint x=0, y=0, z=0, c=0;
		getIndicesFromAddress_field4D(eNum,x,dim0,y,dim1,z,dim2,c,dim3,input.getDataorder());
		Scalar var = pow(double(input.GetImageValue(x,y,z,c) - local_mean.GetImageValue(x,y,z,c)),2);
		aux_mem[eNum] = std::sqrt(var);
	}

	#pragma omp parallel for shared(aux_mem,result)
	for(ulong eNum=0; eNum<numElem; eNum++) {
		uint x=0, y=0, z=0, c=0;
		getIndicesFromAddress_field4D(eNum,x,dim0,y,dim1,z,dim2,c,dim3,result.getDataorder());
		result.setImageValue(aux_mem[eNum],x,y,z,c);
	}

	delete [] aux_mem;
}

#endif



