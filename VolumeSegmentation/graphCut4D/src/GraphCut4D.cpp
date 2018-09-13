/*
 * GraphCut4D.cpp
 *
 *  Created on: Apr 5, 2018
 *      Author: christian
 */

#include "GraphCut4D.hpp"
#include <map>
#include <algorithm>

#define WIDTH 4.0

/* normalize mask so it integrates to one */
void normalizeMask(std::vector<float> &mask) {
  int len = mask.size();
  float sum = 0;
  for (int i = 1; i < len; i++) {
    sum += fabs(mask[i]);
  }
  sum = 2*sum + fabs(mask[0]);
  for (int i = 0; i < len; i++) {
    mask[i] /= sum;
  }
}

MAKE_FILTER(fgauss, exp(-0.5*((i/sigma)*(i/sigma))));

/* convolve src with mask in x-direction */
template<typename Scalar>
void GraphCut4D<Scalar>::convolve_x(Image4D<Scalar> *src, Image4D<Scalar> *dst, std::vector<float> &mask) {
	int width = src->getDimension(0);
	int height = src->getDimension(1);
	int depth = src->getDimension(2);
	int channels = src->getDimension(3);
	int len = mask.size();

	for (int c = 0; c < channels; c++) {
		for (int z = 0; z < depth; z++) {
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					Scalar sum = (Scalar)(mask[0] * src->getImageValue(x,y,z,c));
					for (int i = 1; i < len; i++) {
						sum += (Scalar)(mask[i] * (src->getImageValue(std::max(x-i,0), y, z, c) + src->getImageValue(std::min(x+i, width-1), y, z, c)));
					}
					dst->setImageValue(sum,x,y,z,c);
				}
			}
		}
	}
}

/* convolve src with mask in x-direction */
template<typename Scalar>
void GraphCut4D<Scalar>::convolve_y(Image4D<Scalar> *src, Image4D<Scalar> *dst, std::vector<float> &mask) {
	int width = src->getDimension(0);
	int height = src->getDimension(1);
	int depth = src->getDimension(2);
	int channels = src->getDimension(3);
	int len = mask.size();

	for (int c = 0; c < channels; c++) {
		for (int z = 0; z < depth; z++) {
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					Scalar sum = (Scalar)(mask[0] * src->getImageValue(x,y,z,c));
					for (int i = 1; i < len; i++) {
						sum += (Scalar)(mask[i] * (src->getImageValue(x, std::max(y-i,0), z, c) + src->getImageValue(x, std::min(y+i, height-1), z, c)));
					}
					dst->setImageValue(sum,x,y,z,c);
				}
			}
		}
	}
}

/* convolve src with mask in x-direction */
template<typename Scalar>
void GraphCut4D<Scalar>::convolve_z(Image4D<Scalar> *src, Image4D<Scalar> *dst, std::vector<float> &mask) {
	int width = src->getDimension(0);
	int height = src->getDimension(1);
	int depth = src->getDimension(2);
	int channels = src->getDimension(3);
	int len = mask.size();

	for (int c = 0; c < channels; c++) {
		for (int z = 0; z < depth; z++) {
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					Scalar sum = (Scalar)(mask[0] * src->getImageValue(x,y,z,c));
					for (int i = 1; i < len; i++) {
						sum += (Scalar)(mask[i] * (src->getImageValue(x, y, std::max(z-i,0), c) + src->getImageValue(x, y, std::min(z+i, depth-1), c)));
					}
					dst->setImageValue(sum,x,y,z,c);
				}
			}
		}
	}
}

template<typename Scalar>
void GraphCut4D<Scalar>::smooth(Image4D<Scalar>& input, Image4D<Scalar>& result, float sigma) {
	  std::vector<float> mask = make_fgauss(sigma);
#ifdef DEBUG
	  std::cout << "gaussian mask constructed" << std::endl;
#endif
	  normalizeMask(mask);
#ifdef DEBUG
	  std::cout << "gaussian mask normalised" << std::endl;
#endif

	  Image4D<Scalar> *tmp_x = new Image4D<Scalar>(input.getDimension(0), input.getDimension(1), input.getDimension(2), input.getDimension(3), input.getDatatype(), input.getDataorder());
	  Image4D<Scalar> *tmp_y = new Image4D<Scalar>(input.getDimension(0), input.getDimension(1), input.getDimension(2), input.getDimension(3), input.getDatatype(), input.getDataorder());
	  //Image3D<Scalar> *dst = new Image3D<Scalar>(input.getDimension(0), input.getDimension(1), input.getDimension(2), input.getDatatype(), input.getDataorder());
#ifdef DEBUG
	  std::cout << "containers build" << std::endl;
#endif
	  convolve_x(&input, tmp_x, mask);
	  convolve_y(tmp_x, tmp_y, mask);
	  convolve_z(tmp_y, &result, mask);
#ifdef DEBUG
	  std::cout << "gaussian convolved with image" << std::endl;
#endif

	  tmp_x->clean(); tmp_y->clean();
	  delete tmp_x; delete tmp_y;
}

template<typename Scalar>
Image4D<Scalar>* GraphCut4D<Scalar>::smooth(Image4D<Scalar>& input, float sigma) {
	  std::vector<float> mask = make_fgauss(sigma);
#ifdef DEBUG
	  std::cout << "gaussian mask constructed" << std::endl;
#endif
	  normalizeMask(mask);
#ifdef DEBUG
	  std::cout << "gaussian mask normalised" << std::endl;
#endif

	  Image4D<Scalar> *tmp_x = new Image4D<Scalar>(input.getDimension(0), input.getDimension(1), input.getDimension(2), input.getDimension(3), input.getDatatype(), input.getDataorder());
	  Image4D<Scalar> *tmp_y = new Image4D<Scalar>(input.getDimension(0), input.getDimension(1), input.getDimension(2), input.getDimension(3), input.getDatatype(), input.getDataorder());
	  Image4D<Scalar> *tmp_z = new Image4D<Scalar>(input.getDimension(0), input.getDimension(1), input.getDimension(2), input.getDimension(3), input.getDatatype(), input.getDataorder());
#ifdef DEBUG
	  std::cout << "containers build" << std::endl;
#endif
	  convolve_x(&input, tmp_x, mask);
	  convolve_y(tmp_x, tmp_y, mask);
	  convolve_z(tmp_y, tmp_z, mask);
#ifdef DEBUG
	  std::cout << "gaussian convolved with image" << std::endl;
#endif

	  tmp_x->clean(); tmp_y->clean();
	  delete tmp_x; delete tmp_y;
	  return tmp_z;
}

template<typename Scalar>
GraphCut4D<Scalar>::GraphCut4D() : _weightMethod(SSD) {

}

template<typename Scalar>
GraphCut4D<Scalar>::~GraphCut4D() {

}

//ascending - index pairs
struct asc_value_compare {
  bool operator() (const std::pair<unsigned long, int>& i,const std::pair<unsigned long, int>& j)
  {
	  bool result = false;
	  if(i.second < j.second)
		  result = true;
	  else if(i.second == j.second)
	  {
		  if(i.first < j.first)
			  result = true;
		  else
			  result = false;
	  }
	  else
	  {
		  result = false;
	  }
	  return result;
  }
};

template<typename Scalar>
int GraphCut4D<Scalar>::segmentImageStack_7(Image4D<Scalar>& image, Image3D<short>& output, float sigma, float constantTreshold, int min_size) {
	uint width = image.getDimension(0);
	uint height = image.getDimension(1);
	uint depth = image.getDimension(2);

#ifdef DEBUG
	uint channels = image.getDimension(3);
	std::cout << "w: " << width << ", h:" << height <<", d: " << depth <<", c: " << channels << std::endl;
#endif


	Image4D<Scalar>* smoothVar = smooth(image, sigma);
	//Image4D<Scalar>* smoothVar = new Image4D<Scalar>(image.getDimension(0),image.getDimension(1),image.getDimension(2),image.getDimension(3), image.getDatatype(), image.getDataorder());
	//smoothVar->setData(image.getData(), true);
#ifdef DEBUG
	std::cout << "performed smoothing" << std::endl;
#endif

	// build graph
	unsigned long numVerts = width*height*depth;
	//graphEdge *edges = new graphEdge[numVerts*(9+4+6)];
	std::vector<graphEdge> edges = std::vector<graphEdge>(numVerts*6);
	// FIND OUT THE COMBINATION OF UNIQUE EDGES!! => paper Felzenszwalb2004: undirected edges
	// alternative: add all, then remove duplicates automatically (empty cells at borders), then remove mutual duplicates.
	unsigned long numEdges = 0;
	for (uint y = 0; y < height; y++) {
		for (uint x = 0; x < width; x++) {
			for (uint z = 0; z < depth; z++) {
				if (z > 0) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y,height, z-1,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y,z-1));
					numEdges++;
				}

				if (y < height-1) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y+1,height, z,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y+1,z));
					numEdges++;
				}

				if (x > 0) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y,height, z,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y,z));
					numEdges++;
				}

				if (x < width-1) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y,height, z,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y,z));
					numEdges++;
				}

				if (y > 0) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y-1,height, z,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y-1,z));
					numEdges++;
				}

				if (z < depth-1) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y,z+1));
					numEdges++;
				}
			}
		}
	}
	std::vector<graphEdge> uEdges;
	std::unique_copy(edges.begin(), edges.end(), std::back_inserter(uEdges));
	//std::unique(edges.begin(), edges.end());
	edges.clear();
	numEdges = uEdges.size();
	delete smoothVar;
#ifdef DEBUG
	std::cout << "constructed graph - # edges: " << numEdges << std::endl;
#endif

	// segment
	disjointForest_implicitGrid* u = segmentGraph_implicitGrid(numVerts, numEdges, uEdges, constantTreshold);
#ifdef DEBUG
	std::cout << "segmented graph" << std::endl;
#endif

	// post process small components
	for (ulong i = 0; i < numEdges; i++) {
		ulong a = u->find(uEdges[i].a);
		ulong b = u->find(uEdges[i].b);
		if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
			u->join(a, b);
	}
	//delete [] edges;
	uEdges.clear();
	int num_ccs = (int)(u->num_sets());
#ifdef DEBUG
	std::cout << "fused smaller segments" << std::endl;
#endif

	std::map<ulong, int> componentMap;
	//std::map<ulong, int> sizeMap;
	std::vector< std::pair<ulong, int> > sizeMap;
	std::vector<ulong> uniqueSet = u->unique_set(width*height*depth);
	std::vector<int> uniqueSizeSet = u->unique_sizeSet(width*height*depth);
	for(size_t i=0; i<uniqueSet.size(); i++) {
		componentMap[uniqueSet[i]] = i;
		//sizeMap[uniqueSet[i]] = uniqueSizeSet[i];
		sizeMap.push_back(std::pair<ulong, int>(uniqueSet[i],uniqueSizeSet[i]));
	}
	asc_value_compare _comparator;
	std::sort(sizeMap.begin(), sizeMap.end(), _comparator);
#ifdef DEBUG
	std::cout << "constructed component map" << std::endl;
	//int i=0;
	//for(std::map<ulong, int>::iterator sizeMapIterator = sizeMap.begin(); sizeMapIterator!=sizeMap.end(); sizeMapIterator++) {
	for(int i=0; i<5; i++) {
		std::cout << sizeMap[i].second << std::endl;
	}
#endif
	//image_multichannel<short> *output = new image_multichannel<short>(width, height, 1, true);

	for (uint z = 0; z < depth; z++) {
		for (uint y = 0; y < height; y++) {
			for (uint x = 0; x < width; x++) {
				//ulong comp = u->find(y * width + x);
				ulong address = getDataAddress_field3D(x,width, y,height, z,depth, image.getDataorder());
//#ifdef DEBUG
//				std::cout << "search for connected component addess " << address << " ("<< x <<","<< y <<","<< z <<") ..." << std::endl;
//#endif
				ulong comp = u->find(address);
//#ifdef DEBUG
//				std::cout << "setting segment " << componentMap[comp] << " for ("<< x <<","<< y <<","<< z <<") ..." << std::endl;
//#endif
				output.setImageValue(short(componentMap[comp]), x,y,z);
			}
		}
	}
	delete u;

	return num_ccs;
}

template<typename Scalar>
int GraphCut4D<Scalar>::segmentImageStack_27(Image4D<Scalar>& image, Image3D<short>& output, float sigma, float constantTreshold, int min_size) {
	uint width = image.getDimension(0);
	uint height = image.getDimension(1);
	uint depth = image.getDimension(2);

#ifdef DEBUG
	uint channels = image.getDimension(3);
	std::cout << "w: " << width << ", h:" << height <<", d: " << depth <<", c: " << channels << std::endl;
#endif


	Image4D<Scalar>* smoothVar = smooth(image, sigma);
	//Image4D<Scalar>* smoothVar = new Image4D<Scalar>(image.getDimension(0),image.getDimension(1),image.getDimension(2),image.getDimension(3), image.getDatatype(), image.getDataorder());
	//smoothVar->setData(image.getData(), true);
#ifdef DEBUG
	std::cout << "performed smoothing" << std::endl;
#endif

	// build graph
	unsigned long numVerts = width*height*depth;
	//graphEdge *edges = new graphEdge[numVerts*(9+4+6)];
	std::vector<graphEdge> edges = std::vector<graphEdge>(numVerts*26);
	// FIND OUT THE COMBINATION OF UNIQUE EDGES!! => paper Felzenszwalb2004: undirected edges
	// WEIGHT THE WEIGHTS DEPENDING ON THE KERNEL DISTANCE!!! (1-distance) rather -> higher weights mean higher difference
	// alternative: add all, then remove duplicates automatically (empty cells at borders), then remove mutual duplicates.
	unsigned long numEdges = 0;
	for (uint y = 0; y < height; y++) {
		for (uint x = 0; x < width; x++) {
			for (uint z = 0; z < depth; z++) {
				// front-stack
				if ((x > 0) && (y < height-1) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y+1,height, z-1,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y+1,z-1));
					numEdges++;
				}

				if ((y < height-1) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y+1,height, z-1,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y+1,z-1));
					numEdges++;
				}

				if ((x < width-1) && (y < height-1) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y+1,height, z-1,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y+1,z-1));
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if ((x > 0) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y,height, z-1,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y,z-1));
					numEdges++;
				}

				if (z > 0) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y,height, z-1,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y,z-1));
					numEdges++;
				}

				if ((x < width-1) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y,height, z-1,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y,z-1));
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if ((x > 0) && (y > 0) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y-1,height, z-1,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y-1,z-1));
					numEdges++;
				}

				if ((y > 0) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y-1,height, z-1,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y-1,z-1));
					numEdges++;
				}

				if ((x < width-1) && (y > 0) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y-1,height, z-1,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y-1,z-1));
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



				// mid-stack
				if ((x > 0) && (y < height-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y+1,height, z,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y+1,z));
					numEdges++;
				}

				if (y < height-1) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y+1,height, z,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y+1,z));
					numEdges++;
				}

				if ((x < width-1) && (y < height-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y+1,height, z,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y+1,z));
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if (x > 0) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y,height, z,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y,z));
					numEdges++;
				}

				if (x < width-1) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y,height, z,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y,z));
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if ((x > 0) && (y > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y-1,height, z,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y-1,z));
					numEdges++;
				}

				if (y > 0) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y-1,height, z,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y-1,z));
					numEdges++;
				}

				if ((x < width-1) && (y > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y-1,height, z,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y-1,z));
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



				//back-stack
				if ((x > 0) && (y < height-1) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y+1,height, z+1,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y+1,z+1));
					numEdges++;
				}

				if ((y < height-1) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y+1,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y+1,z+1));
					numEdges++;
				}

				if ((x < width-1) && (y < height-1) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y+1,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y+1,z+1));
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if ((x > 0) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y,height, z+1,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y,z+1));
					numEdges++;
				}

				if (z < depth-1) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y,z+1));
					numEdges++;
				}

				if ((x < width-1) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y,z+1));
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if ((x > 0) && (y > 0) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y-1,height, z+1,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y-1,z+1));
					numEdges++;
				}

				if ((y > 0) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y-1,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y-1,z+1));
					numEdges++;
				}

				if ((x < width-1) && (y > 0) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y-1,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y-1,z+1));
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			}
		}
	}
	std::vector<graphEdge> uEdges;
	std::unique_copy(edges.begin(), edges.end(), std::back_inserter(uEdges));
	//std::unique(edges.begin(), edges.end());
	edges.clear();
	numEdges = uEdges.size();
	delete smoothVar;
#ifdef DEBUG
	std::cout << "constructed graph - # edges: " << numEdges << std::endl;
#endif

	// segment
	disjointForest_implicitGrid* u = segmentGraph_implicitGrid(numVerts, numEdges, uEdges, constantTreshold);
#ifdef DEBUG
	std::cout << "segmented graph" << std::endl;
#endif

	// post process small components
	for (ulong i = 0; i < numEdges; i++) {
		ulong a = u->find(uEdges[i].a);
		ulong b = u->find(uEdges[i].b);
		if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
			u->join(a, b);
	}
	//delete [] edges;
	uEdges.clear();
	int num_ccs = (int)(u->num_sets());
#ifdef DEBUG
	std::cout << "fused smaller segments" << std::endl;
#endif

	std::map<ulong, int> componentMap;
	std::vector<ulong> uniqueSet = u->unique_set(width*height*depth);
	for(size_t i=0; i<uniqueSet.size(); i++) {
		componentMap[uniqueSet[i]] = i;
	}
#ifdef DEBUG
	std::cout << "constructed component map - " << uniqueSet.size() << " unique components." << std::endl;
#endif
	//image_multichannel<short> *output = new image_multichannel<short>(width, height, 1, true);

	for (uint z = 0; z < depth; z++) {
		for (uint y = 0; y < height; y++) {
			for (uint x = 0; x < width; x++) {
				//ulong comp = u->find(y * width + x);
				ulong address = getDataAddress_field3D(x,width, y,height, z,depth, image.getDataorder());
//#ifdef DEBUG
//				std::cout << "search for connected component addess " << address << " ("<< x <<","<< y <<","<< z <<") ..." << std::endl;
//#endif
				ulong comp = u->find(address);
//#ifdef DEBUG
//				std::cout << "setting segment " << componentMap[comp] << " for ("<< x <<","<< y <<","<< z <<") ..." << std::endl;
//#endif
				output.setImageValue(short(componentMap[comp]), x,y,z);
			}
		}
	}
	delete u;

	return num_ccs;
}

template<typename Scalar>
int GraphCut4D<Scalar>::segmentImageStack_27weighted(Image4D<Scalar>& image, Image3D<short>& output, float sigma, float constantTreshold, int min_size) {
	uint width = image.getDimension(0);
	uint height = image.getDimension(1);
	uint depth = image.getDimension(2);

#ifdef DEBUG
	uint channels = image.getDimension(3);
	std::cout << "w: " << width << ", h:" << height <<", d: " << depth <<", c: " << channels << std::endl;
#endif


	Image4D<Scalar>* smoothVar = smooth(image, sigma);
	//Image4D<Scalar>* smoothVar = new Image4D<Scalar>(image.getDimension(0),image.getDimension(1),image.getDimension(2),image.getDimension(3), image.getDatatype(), image.getDataorder());
	//smoothVar->setData(image.getData(), true);
#ifdef DEBUG
	std::cout << "performed smoothing" << std::endl;
#endif

	// build graph
	unsigned long numVerts = width*height*depth;
	//graphEdge *edges = new graphEdge[numVerts*(9+4+6)];
	std::vector<graphEdge> edges = std::vector<graphEdge>(numVerts*26);
	// FIND OUT THE COMBINATION OF UNIQUE EDGES!! => paper Felzenszwalb2004: undirected edges
	// WEIGHT THE WEIGHTS DEPENDING ON THE KERNEL DISTANCE!!! (1-distance) rather -> higher weights mean higher difference
	// alternative: add all, then remove duplicates automatically (empty cells at borders), then remove mutual duplicates.
	unsigned long numEdges = 0;
	for (uint y = 0; y < height; y++) {
		for (uint x = 0; x < width; x++) {
			for (uint z = 0; z < depth; z++) {
				// front-stack
				if ((x > 0) && (y < height-1) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y+1,height, z-1,depth, smoothVar->getDataorder());
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y+1,z-1))*(1.0f-0.005f);
					numEdges++;
				}

				if ((y < height-1) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());
					edges[numEdges].b = getDataAddress_field3D(x,width, y+1,height, z-1,depth, smoothVar->getDataorder());
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y+1,z-1))*(1.0f-0.18f);
					numEdges++;
				}

				if ((x < width-1) && (y < height-1) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y+1,height, z-1,depth, smoothVar->getDataorder());
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y+1,z-1))*(1.0f-0.005f);
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if ((x > 0) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y,height, z-1,depth, smoothVar->getDataorder());
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y,z-1))*(1.0f-0.18f);
					numEdges++;
				}

				if (z > 0) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());
					edges[numEdges].b = getDataAddress_field3D(x,width, y,height, z-1,depth, smoothVar->getDataorder());
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y,z-1))*(1.0f-0.92f);
					numEdges++;
				}

				if ((x < width-1) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y,height, z-1,depth, smoothVar->getDataorder());
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y,z-1))*(1.0f-0.18f);
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if ((x > 0) && (y > 0) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y-1,height, z-1,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y-1,z-1))*(1.0f-0.005f);
					numEdges++;
				}

				if ((y > 0) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y-1,height, z-1,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y-1,z-1))*(1.0f-0.18f);
					numEdges++;
				}

				if ((x < width-1) && (y > 0) && (z > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y-1,height, z-1,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y-1,z-1))*(1.0f-0.005f);
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



				// mid-stack
				if ((x > 0) && (y < height-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y+1,height, z,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y+1,z))*(1.0f-0.18f);
					numEdges++;
				}

				if (y < height-1) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y+1,height, z,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y+1,z))*(1.0f-0.92f);
					numEdges++;
				}

				if ((x < width-1) && (y < height-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y+1,height, z,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y+1,z))*(1.0f-0.18f);
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if (x > 0) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y,height, z,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y,z))*(1.0f-0.92f);
					numEdges++;
				}

				if (x < width-1) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y,height, z,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y,z))*(1.0f-0.92f);
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if ((x > 0) && (y > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y-1,height, z,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y-1,z))*(1.0f-0.18f);
					numEdges++;
				}

				if (y > 0) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y-1,height, z,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y-1,z))*(1.0f-0.92f);
					numEdges++;
				}

				if ((x < width-1) && (y > 0)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y-1,height, z,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y-1,z))*(1.0f-0.18f);
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



				//back-stack
				if ((x > 0) && (y < height-1) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y+1,height, z+1,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y+1,z+1))*(1.0f-0.005f);
					numEdges++;
				}

				if ((y < height-1) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y+1,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y+1,z+1))*(1.0f-0.18f);
					numEdges++;
				}

				if ((x < width-1) && (y < height-1) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y+1,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y+1,z+1))*(1.0f-0.005f);
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if ((x > 0) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y,height, z+1,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y,z+1))*(1.0f-0.18f);
					numEdges++;
				}

				if (z < depth-1) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y,z+1))*(1.0f-0.92f);
					numEdges++;
				}

				if ((x < width-1) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y,z+1))*(1.0f-0.18f);
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if ((x > 0) && (y > 0) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x-1,width, y-1,height, z+1,depth, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x-1,y-1,z+1))*(1.0f-0.005f);
					numEdges++;
				}

				if ((y > 0) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x,width, y-1,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x,y-1,z+1))*(1.0f-0.18f);
					numEdges++;
				}

				if ((x < width-1) && (y > 0) && (z < depth-1)) {
					edges[numEdges].a = getDataAddress_field3D(x,width, y,height, z,depth, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field3D(x+1,width, y-1,height, z+1,depth, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y,z, x+1,y-1,z+1))*(1.0f-0.005f);
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			}
		}
	}
	std::vector<graphEdge> uEdges;
	std::unique_copy(edges.begin(), edges.end(), std::back_inserter(uEdges));
	//std::unique(edges.begin(), edges.end());
	edges.clear();
	numEdges = uEdges.size();
	delete smoothVar;
#ifdef DEBUG
	std::cout << "constructed graph - # edges: " << numEdges << std::endl;
#endif

	// segment
	disjointForest_implicitGrid* u = segmentGraph_implicitGrid(numVerts, numEdges, uEdges, constantTreshold);
#ifdef DEBUG
	std::cout << "segmented graph" << std::endl;
#endif

	// post process small components
	for (ulong i = 0; i < numEdges; i++) {
		ulong a = u->find(uEdges[i].a);
		ulong b = u->find(uEdges[i].b);
		if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
			u->join(a, b);
	}
	//delete [] edges;
	uEdges.clear();
	int num_ccs = (int)(u->num_sets());
#ifdef DEBUG
	std::cout << "fused smaller segments" << std::endl;
#endif

	std::map<ulong, int> componentMap;
	std::map<ulong, int> sizeMap;
	std::vector<ulong> uniqueSet = u->unique_set(width*height*depth);
	std::vector<int> uniqueSizeSet = u->unique_sizeSet(width*height*depth);
	for(size_t i=0; i<uniqueSet.size(); i++) {
		componentMap[uniqueSet[i]] = i;
		sizeMap[uniqueSet[i]] = uniqueSizeSet[i];
	}

#ifdef DEBUG
	std::cout << "constructed component map - " << uniqueSet.size() << " unique components." << std::endl;
#endif
	//image_multichannel<short> *output = new image_multichannel<short>(width, height, 1, true);

	for (uint z = 0; z < depth; z++) {
		for (uint y = 0; y < height; y++) {
			for (uint x = 0; x < width; x++) {
				//ulong comp = u->find(y * width + x);
				ulong address = getDataAddress_field3D(x,width, y,height, z,depth, image.getDataorder());
//#ifdef DEBUG
//				std::cout << "search for connected component addess " << address << " ("<< x <<","<< y <<","<< z <<") ..." << std::endl;
//#endif
				ulong comp = u->find(address);
//#ifdef DEBUG
//				std::cout << "setting segment " << componentMap[comp] << " for ("<< x <<","<< y <<","<< z <<") ..." << std::endl;
//#endif
				output.setImageValue(short(componentMap[comp]), x,y,z);
			}
		}
	}
	delete u;

	return num_ccs;
}

template class GraphCut4D<char>;
template class GraphCut4D<unsigned char>;
template class GraphCut4D<short>;
template class GraphCut4D<unsigned short>;
template class GraphCut4D<int>;
template class GraphCut4D<unsigned int>;
template class GraphCut4D<long>;
template class GraphCut4D<unsigned long>;
template class GraphCut4D<float>;
template class GraphCut4D<double>;




