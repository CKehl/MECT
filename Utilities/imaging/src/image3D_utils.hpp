/**
 * Created on 20 Mar 2018
 * author: christian
 */

#ifndef IMAGE3D_UTILS_HPP_
#define IMAGE3D_UTILS_HPP_

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

typedef enum _NEIGHBOURHOOD3D {
	NBR_7   = 0,
	NBR_27  = 1,
	NBR_27W = 2,
} NEIGHBOURHOOD3D;

// ======== ======== ======== ======== ======== ======== ======== //
//                   Image3D Support Functions                    //
// ======== ======== ======== ======== ======== ======== ======== //

template<typename Scalar>
inline void computeGradient(Image3D<Scalar>& input, Image3D<Scalar>& result) {
	//Image3D(uint* dimensions, dtype datatype, dorder dataorder, Scalar* originalData=NULL);
	//Image3D<Scalar> result(input.getDimensions(), input.getDatatype(), input.getDataorder());
	long numElem = result.getDimension(0)*result.getDimension(1)*result.getDimension(2);
	uint dim0 = result.getDimension(0); uint dim1 = result.getDimension(1); uint dim2 = result.getDimension(2);
	Scalar* aux_mem;
	aux_mem = new Scalar[numElem];


	#pragma omp parallel for
	for(long eNum=0; eNum<numElem; eNum++) {
		uint z = uint(std::floor(eNum/(dim0*dim1)));
		long rest = eNum % (dim0*dim1);
		uint y = uint(std::floor(rest/dim0));
		uint x = uint(rest%dim0);
		Scalar kernel[3][3][3];
		uint kernel_k=0, kernel_l=0, kernel_m=0;
		for(int k=(int(x)-1); k<(int(x)+2); k++) {
			for(int l=(int(y)-1); l<(int(y)+2); l++) {
				for(int m=(int(z)-1); m<(int(z)+2); m++) {
					//zero-based addresses
					kernel_k=k-(x-1);
					kernel_l=l-(y-1);
					kernel_m=m-(z-1);
					if((k<0)||(uint(k)>(dim0-1)) || (l<0)||(uint(l)>(dim1-1)) || (m<0)||(uint(m)>(dim2-1))) {
						kernel[kernel_k][kernel_l][kernel_m]=0;
					} else {
						kernel[kernel_k][kernel_l][kernel_m]=input.GetImageValue(k,l,m);
					}
				}
			}
		}
		Scalar grad = 0;
		grad += (kernel[0][0][2] - kernel[0][0][0]);
		grad += (kernel[0][1][2] - kernel[0][1][0]);
		grad += (kernel[0][2][2] - kernel[0][2][0]);

		grad += (kernel[1][0][2] - kernel[1][0][0]);
		grad += (kernel[1][1][2] - kernel[1][1][0]);
		grad += (kernel[1][2][2] - kernel[1][2][0]);

		grad += (kernel[2][0][2] - kernel[2][0][0]);
		grad += (kernel[2][1][2] - kernel[2][1][0]);
		grad += (kernel[2][2][2] - kernel[2][2][0]);

		grad += (kernel[0][2][0] - kernel[0][0][0]);
		grad += (kernel[0][2][1] - kernel[0][0][1]);
		grad += (kernel[0][2][2] - kernel[0][0][2]);

		grad += (kernel[1][2][0] - kernel[1][0][0]);
		grad += (kernel[1][2][1] - kernel[1][0][1]);
		grad += (kernel[1][2][2] - kernel[1][0][2]);

		grad += (kernel[2][2][0] - kernel[2][0][0]);
		grad += (kernel[2][2][1] - kernel[2][0][1]);
		grad += (kernel[2][2][2] - kernel[2][0][2]);

		grad += (kernel[2][0][0] - kernel[0][0][0]);
		grad += (kernel[2][0][1] - kernel[0][0][1]);
		grad += (kernel[2][0][2] - kernel[0][0][2]);

		grad += (kernel[2][1][0] - kernel[0][1][0]);
		grad += (kernel[2][1][1] - kernel[0][1][1]);
		grad += (kernel[2][1][2] - kernel[0][1][2]);

		grad += (kernel[2][2][0] - kernel[0][2][0]);
		grad += (kernel[2][2][1] - kernel[0][2][1]);
		grad += (kernel[2][2][2] - kernel[0][2][2]);

		//grad += (kernel[0][0][2] - kernel[2][2][0]);
		//grad += (kernel[0][1][2] - kernel[2][1][0]);
		//grad += (kernel[0][2][2] - kernel[2][0][0]);

		//grad += (kernel[2][2][2] - kernel[0][0][0]);
		//grad += (kernel[2][1][2] - kernel[0][1][0]);
		//grad += (kernel[2][0][2] - kernel[0][2][0]);

		grad=grad/27;
		aux_mem[eNum] = grad;
	}

	for(uint x=0; x<dim0; x++) {
		for(uint y=0; y<dim1; y++) {
			for(uint z=0; z<dim2; z++) {
				long eNum = (z*(dim0*dim1))+(y*dim0)+x;
				result.setImageValue(aux_mem[eNum],x,y,z);
			}
		}
	}
}

template<typename Scalar>
inline void computeGaussianSmoothing(Image3D<Scalar>& input, Image3D<Scalar>& result, float sigma) {

}

template<typename Scalar>
std::vector< std::vector<Scalar> > computeResponseCurves_uniformSample(Image3D<Scalar>& input, uint stride) {
	std::vector< std::vector<Scalar> > result;
	for(uint x=0; x<input.getDimension(0); x+=stride) {
		for(uint y=0; y<input.getDimension(1); y+=stride) {
			std::vector<Scalar> sample;
			for(uint channel=0; channel<input.getDimension(2); channel++) {
				sample.push_back(input.getImageValue(x,y,channel));
			}
			result.push_back(sample);
		}
	}
	return result;
}

template<typename Scalar>
void computeLocalMean(Image3D<Scalar>& input, Image3D<Scalar>& result) {
	//Image3D(uint* dimensions, dtype datatype, dorder dataorder, Scalar* originalData=NULL);
	//Image3D<Scalar> result(input.getDimensions(), input.getDatatype(), input.getDataorder());
	long numElem = result.getDimension(0)*result.getDimension(1)*result.getDimension(2);
	uint dim0 = result.getDimension(0); uint dim1 = result.getDimension(1); uint dim2 = result.getDimension(2);
	Scalar* aux_mem;
	aux_mem = new Scalar[numElem];


	#pragma omp parallel for
	for(long eNum=0; eNum<numElem; eNum++) {
		uint z = uint(std::floor(eNum/(dim0*dim1)));
		long rest = eNum % (dim0*dim1);
		uint y = uint(std::floor(rest/dim0));
		uint x = uint(rest%dim0);
		//Scalar kernel[7][7][7];
		//uint kernel_k=0, kernel_l=0, kernel_m=0;
		Scalar lmean = 0;
		//for(int k=(int(x)-3); k<(int(x)+4); k++) {
		//	for(int l=(int(y)-3); l<(int(y)+4); l++) {
		//		for(int m=(int(z)-3); m<(int(z)+4); m++) {
		for(int k=(int(x)-1); k<(int(x)+2); k++) {
			for(int l=(int(y)-1); l<(int(y)+2); l++) {
				for(int m=(int(z)-1); m<(int(z)+2); m++) {
					//kernel_k=k-(x-1);
					//kernel_l=l-(y-1);
					//kernel_m=m-(z-1);
					if((k<0)||(uint(k)>(dim0-1)) || (l<0)||(uint(l)>(dim1-1)) || (m<0)||(uint(m)>(dim2-1))) {
					//	kernel[kernel_k][kernel_l][kernel_m]=0;
						lmean+=0;
					} else {
					//	kernel[kernel_k][kernel_l][kernel_m]=input.GetImageValue(k,l,m);
						lmean += input.GetImageValue(k,l,m);
					}
				}
			}
		}

		//lmean /= (7*7*7);
		lmean /= (3*3*3);
		aux_mem[eNum] = lmean;
	}

	for(uint x=0; x<dim0; x++) {
		for(uint y=0; y<dim1; y++) {
			for(uint z=0; z<dim2; z++) {
				long eNum = (z*(dim0*dim1))+(y*dim0)+x;
				result.setImageValue(aux_mem[eNum],x,y,z);
			}
		}
	}

	delete [] aux_mem;
}

template<typename Scalar>
void computeLocalMean(Image3D<Scalar>& input, Image3D<Scalar>& result, uint kernelSize) {
	int K2 = int(std::floor(kernelSize/2));
	int K2P1 = K2+1;

	ulong numElem = result.getDimension(0)*result.getDimension(1)*result.getDimension(2);
	uint dim0 = result.getDimension(0); uint dim1 = result.getDimension(1); uint dim2 = result.getDimension(2);
	Scalar* aux_mem;
	aux_mem = new Scalar[numElem];


	#pragma omp parallel for shared(input,aux_mem)
	for(ulong eNum=0; eNum<numElem; eNum++) {
		uint x=0, y=0, z=0;
		getIndicesFromAddress_field3D(eNum,x,dim0,y,dim1,z,dim2,input.getDataorder());
		Scalar lmean = 0;
		for(int k=(int(x)-K2); k<(int(x)+K2P1); k++) {
			for(int l=(int(y)-K2); l<(int(y)+K2P1); l++) {
				//for(int m=(int(z)-K2); m<(int(z)+K2P1); m++) {
				// || (m<0)||(uint(m)>(dim2-1))
					if((k<0)||(uint(k)>(dim0-1)) || (l<0)||(uint(l)>(dim1-1))) {
						lmean+=0;
					} else {
						//lmean += input.GetImageValue(k,l,m);
						lmean += input.GetImageValue(k,l,z);
					}
				//}
			}
		}
		lmean /= (kernelSize*kernelSize*kernelSize);
		aux_mem[eNum] = lmean;
	}

	#pragma omp parallel for shared(aux_mem,result)
	for(ulong eNum=0; eNum<numElem; eNum++) {
		uint x=0, y=0, z=0;
		getIndicesFromAddress_field3D(eNum,x,dim0,y,dim1,z,dim2,result.getDataorder());
		result.setImageValue(aux_mem[eNum],x,y,z);
	}

	delete [] aux_mem;
}

template<typename Scalar>
void computeLocalStdDev(Image3D<Scalar>& input, Image3D<Scalar>& local_mean, Image3D<Scalar>& result) {
	//Image3D<Scalar> result(input.getDimensions(), input.getDatatype(), input.getDataorder());
	long numElem = result.getDimension(0)*result.getDimension(1)*result.getDimension(2);
	uint dim0 = result.getDimension(0); uint dim1 = result.getDimension(1); uint dim2 = result.getDimension(2);
	Scalar* aux_mem;
	aux_mem = new Scalar[numElem];


	#pragma omp parallel for
	for(long eNum=0; eNum<numElem; eNum++) {
		uint z = uint(std::floor(eNum/(dim0*dim1)));
		long rest = eNum % (dim0*dim1);
		uint y = uint(std::floor(rest/dim0));
		uint x = uint(rest%dim0);
		Scalar var = 0;
		//for(int k=(int(x)-1); k<(int(x)+2); k++) {
		//	for(int l=(int(y)-1); l<(int(y)+2); l++) {
		//		for(int m=(int(z)-1); m<(int(z)+2); m++) {
		for(int k=(int(x)-2); k<(int(x)+3); k++) {
			for(int l=(int(y)-2); l<(int(y)+3); l++) {
				for(int m=(int(z)-2); m<(int(z)+3); m++) {
		//for(int k=(int(x)-3); k<(int(x)+4); k++) {
		//	for(int l=(int(y)-3); l<(int(y)+4); l++) {
		//		for(int m=(int(z)-3); m<(int(z)+4); m++) {
					if((k<0)||(uint(k)>(dim0-1)) || (l<0)||(uint(l)>(dim1-1)) || (m<0)||(uint(m)>(dim2-1))) {
					//	kernel[kernel_k][kernel_l][kernel_m]=0;
						var+=0;
					} else {
						var += pow(double(input.GetImageValue(k,l,m) - local_mean.GetImageValue(k,l,m)),2);
					}
				}
			}
		}

		//var /= (7*7*7);
		var /= (5*5*5);
		//var /= (3*3*3);
		aux_mem[eNum] = std::sqrt(var);
	}

	for(uint x=0; x<dim0; x++) {
		for(uint y=0; y<dim1; y++) {
			for(uint z=0; z<dim2; z++) {
				long eNum = (z*(dim0*dim1))+(y*dim0)+x;
				result.setImageValue(aux_mem[eNum],x,y,z);
			}
		}
	}

	delete [] aux_mem;
}

template<typename Scalar>
void computeLocalStdDev(Image3D<Scalar>& input, Image3D<Scalar>& local_mean, Image3D<Scalar>& result, uint kernelSize) {
	int K2 = int(std::floor(kernelSize/2));
	int K2P1 = K2+1;
	ulong numElem = result.getDimension(0)*result.getDimension(1)*result.getDimension(2);
	uint dim0 = result.getDimension(0); uint dim1 = result.getDimension(1); uint dim2 = result.getDimension(2);
	Scalar* aux_mem;
	aux_mem = new Scalar[numElem];


	#pragma omp parallel for shared(aux_mem,local_mean,input)
	for(ulong eNum=0; eNum<numElem; eNum++) {
		uint x=0, y=0, z=0;
		getIndicesFromAddress_field3D(eNum,x,dim0,y,dim1,z,dim2,input.getDataorder());
		Scalar var = 0;
		for(int k=(int(x)-K2); k<(int(x)+K2P1); k++) {
			for(int l=(int(y)-K2); l<(int(y)+K2P1); l++) {
				//for(int m=(int(z)-K2); m<(int(z)+K2P1); m++) {
					//if((k<0)||(uint(k)>(dim0-1)) || (l<0)||(uint(l)>(dim1-1)) || (m<0)||(uint(m)>(dim2-1))) {
					if((k<0)||(uint(k)>(dim0-1)) || (l<0)||(uint(l)>(dim1-1))) {
						var+=0;
					} else {
						//var += pow(double(input.GetImageValue(k,l,m) - local_mean.GetImageValue(k,l,m)),2);
						var += pow(double(input.GetImageValue(k,l,z) - local_mean.GetImageValue(k,l,z)),2);
					}
				//}
			}
		}
		var /= (kernelSize*kernelSize*kernelSize);
		aux_mem[eNum] = std::sqrt(var);
	}

	#pragma omp parallel for shared(aux_mem,result)
	for(ulong eNum=0; eNum<numElem; eNum++) {
		uint x=0, y=0, z=0;
		getIndicesFromAddress_field3D(eNum,x,dim0,y,dim1,z,dim2,result.getDataorder());
		result.setImageValue(aux_mem[eNum],x,y,z);
	}

	delete [] aux_mem;
}

template<typename Scalar>
void computeDeviation(Image3D<Scalar>& input, Image3D<Scalar>& local_mean, Image3D<Scalar>& result) {
	ulong numElem = result.getDimension(0)*result.getDimension(1)*result.getDimension(2);
	uint dim0 = result.getDimension(0); uint dim1 = result.getDimension(1); uint dim2 = result.getDimension(2);
	Scalar* aux_mem;
	aux_mem = new Scalar[numElem];


	#pragma omp parallel for shared(aux_mem,local_mean,input)
	for(ulong eNum=0; eNum<numElem; eNum++) {
		uint x=0, y=0, c=0;
		getIndicesFromAddress_field3D(eNum,x,dim0,y,dim1,c,dim2,input.getDataorder());
		Scalar var = pow(double(input.GetImageValue(x,y,c) - local_mean.GetImageValue(x,y,c)),2);
		aux_mem[eNum] = std::sqrt(var);
	}

	#pragma omp parallel for shared(aux_mem,result)
	for(ulong eNum=0; eNum<numElem; eNum++) {
		uint x=0, y=0,c=0;
		getIndicesFromAddress_field3D(eNum,x,dim0,y,dim1,c,dim2,result.getDataorder());
		result.setImageValue(aux_mem[eNum],x,y,c);
	}

	delete [] aux_mem;
}

template<typename Scalar>
inline void computeLocalStdDev(Image3D<Scalar>& input, Image3D<Scalar>& result) {
	Image3D<Scalar>* lmean = new Image3D<Scalar>(input.getDimensions(), input.getDatatype(), input.getDataorder());
	lmean->createImage();
	computeLocalMean(input, *lmean);
	computeLocalStdDev(input, *lmean, result);
	//delete lmean;
}

template<typename Scalar>
inline MinMaxUI locateOuterBoundary(Image3D<Scalar>& input, Scalar bVal) {
	MinMaxUI result;

	MinMaxUI lowStackBounds, highStackBounds;
	uint dim0 = input.getDimension(0); uint dim1 = input.getDimension(1); uint dim2 = input.getDimension(2);
	lowStackBounds.min_x=dim0; lowStackBounds.max_x=0; lowStackBounds.min_y=dim1; lowStackBounds.max_y=0; lowStackBounds.min_z=dim2; lowStackBounds.max_z=0;
	highStackBounds.min_x=dim0; highStackBounds.max_x=0; highStackBounds.min_y=dim1; highStackBounds.max_y=0; highStackBounds.min_z=dim2; highStackBounds.max_z=0;
	result.min_x=dim0; result.max_x=0; result.min_y=dim1; result.max_y=0; result.min_z=dim2; result.max_z=0;
	bool change=false; Scalar avg=0;
	for(uint z=0; z<dim2; z++) {
		change=false;
		avg=0;
		for(uint y=0; y<dim1; y++) {
			for(uint x=0; x<dim0; x++) {
				Scalar cVal = input.GetImageValue(x,y,z);
				avg+=cVal;
				if(cVal<=bVal) {
					if((x*x+y*y) < (lowStackBounds.min_x*lowStackBounds.min_x+lowStackBounds.min_y*lowStackBounds.min_y)) {
						//optVal = cVal;
						lowStackBounds.min_x=x;
						lowStackBounds.min_y=y;
						lowStackBounds.min_z=z;
						change=true;
					}
					if((x*x+y*y) > (lowStackBounds.max_x*lowStackBounds.max_x+lowStackBounds.max_y*lowStackBounds.max_y)) {
						//optVal = cVal;
						lowStackBounds.max_x=x;
						lowStackBounds.max_y=y;
						lowStackBounds.max_z=z;
						change=true;
					}
				}
			}
		}
		avg/=(dim1*dim0);
		if((change==false) && (avg<=(bVal*1.05))) {
			lowStackBounds.min_z=z;
			lowStackBounds.max_z=z;
			break;
		}
	}
	for(uint z=dim2-1; z>0; --z) {
		change=false;
		avg=0;
		for(uint y=0; y<dim1; y++) {
			for(uint x=0; x<dim0; x++) {
				Scalar cVal = input.GetImageValue(x,y,z);
				avg+=cVal;
				if(cVal<=bVal) {
					if((x*x+y*y) < (highStackBounds.min_x*highStackBounds.min_x+highStackBounds.min_y*highStackBounds.min_y)) {
						//optVal = cVal;
						highStackBounds.min_x=x;
						highStackBounds.min_y=y;
						highStackBounds.min_z=z;
						change=true;
					}
					if((x*x+y*y) > (highStackBounds.max_x*highStackBounds.max_x+highStackBounds.max_y*highStackBounds.max_y)) {
						//optVal = cVal;
						highStackBounds.max_x=x;
						highStackBounds.max_y=y;
						highStackBounds.max_z=z;
						change=true;
					}
				}
			}
		}
		avg/=(dim1*dim0);
		if((change==false) && (avg<=(bVal*1.05))) {
			highStackBounds.min_z=z;
			highStackBounds.max_z=z;
			break;
		}
	}

	result.min_x = std::max(lowStackBounds.min_x, highStackBounds.min_x);
	result.min_y = std::max(lowStackBounds.min_y, highStackBounds.min_y);
	result.max_x = std::min(lowStackBounds.max_x, highStackBounds.max_x);
	result.max_y = std::min(lowStackBounds.max_y, highStackBounds.max_y);
	result.max_z=highStackBounds.max_z;
	result.min_z=lowStackBounds.min_z;
	return result;
}

// dissimilarity measure between pixels
template<typename Scalar>
inline Scalar SSD_imageStack(Image3D<Scalar>& data, int x1, int y1, int x2, int y2) {
	float result = 0;
	for(unsigned int c = 0; c < data.getDimension(2); c++) {
		result += (data.getImageValue(x1, y1, c)-data.getImageValue(x2, y2, c))*(data.getImageValue(x1, y1, c)-data.getImageValue(x2, y2, c));
	}
	return std::sqrt(result);
}

template<typename Scalar>
inline void createNonUniformAverageImage(Image3D<Scalar>& input, std::vector< std::pair<uint, uint> >& bounds, Image3D<Scalar>& output) {
	uint dim0=input.getDimension(0), dim1=input.getDimension(1);
	for(uint boundaryIdx=0; boundaryIdx<bounds.size(); boundaryIdx++) {
		uint minChannel = bounds[boundaryIdx].first;
		uint maxChannel = bounds[boundaryIdx].second;
		#pragma omp parallel for shared(input,output,minChannel,maxChannel,boundaryIdx,dim0,dim1)
		for(ulong i=0; i<(dim0*dim1); i++) {
			uint x=0, y=0;
			getIndicesFromAddress_field2D(i,x,dim0,y,dim1,input.getDataorder());
			double avgVal = 0;
			for(uint c=minChannel; c<maxChannel; c++) {
				avgVal+=input.getImageValue(x,y,c);
			}
			avgVal /= double(maxChannel-minChannel);
			output.setImageValue(avgVal, x,y,boundaryIdx);
		}
	}
}


#endif /* IMAGE3D_UTILS_HPP_ */
