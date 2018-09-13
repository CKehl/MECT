/*
 * GraphCut3D.cpp
 *
 *  Created on: Apr 5, 2018
 *      Author: christian
 */

#include "GraphCut3D.hpp"
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
void GraphCut3D<Scalar>::convolve_x(Image3D<Scalar> *src, Image3D<Scalar> *dst, std::vector<float> &mask) {
	int width = src->getDimension(0);
	int height = src->getDimension(1);
	int channels = src->getDimension(2);
	int len = mask.size();

	for (int c = 0; c < channels; c++) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				Scalar sum = (Scalar)(mask[0] * src->getImageValue(x,y,c));
				for (int i = 1; i < len; i++) {
					sum += (Scalar)(mask[i] * (src->getImageValue(std::max(x-i,0), y, c) + src->getImageValue(std::min(x+i, width-1), y, c)));
				}
				dst->setImageValue(sum,x,y,c);
			}
		}
	}
}

/* convolve src with mask in x-direction */
template<typename Scalar>
void GraphCut3D<Scalar>::convolve_y(Image3D<Scalar> *src, Image3D<Scalar> *dst, std::vector<float> &mask) {
	int width = src->getDimension(0);
	int height = src->getDimension(1);
	int channels = src->getDimension(2);
	int len = mask.size();

	for (int c = 0; c < channels; c++) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				Scalar sum = (Scalar)(mask[0] * src->getImageValue(x,y,c));
				for (int i = 1; i < len; i++) {
					sum += (Scalar)(mask[i] * (src->getImageValue(x, std::max(y-i,0), c) + src->getImageValue(x, std::min(y+i, height-1), c)));
				}
				dst->setImageValue(sum,x,y,c);
			}
		}
	}
}

template<typename Scalar>
void GraphCut3D<Scalar>::smooth(Image3D<Scalar>& input, Image3D<Scalar>& result, float sigma) {
	  std::vector<float> mask = make_fgauss(sigma);
#ifdef DEBUG
	  std::cout << "gaussian mask constructed" << std::endl;
#endif
	  normalizeMask(mask);
#ifdef DEBUG
	  std::cout << "gaussian mask normalised" << std::endl;
#endif

	  Image3D<Scalar> *tmp_x = new Image3D<Scalar>(input.getDimension(0), input.getDimension(1), input.getDimension(2), input.getDatatype(), input.getDataorder());
	  //Image3D<Scalar> *tmp_y = new Image3D<Scalar>(input.getDimension(0), input.getDimension(1), input.getDimension(2), input.getDatatype(), input.getDataorder());
	  //Image3D<Scalar> *dst = new Image3D<Scalar>(input.getDimension(0), input.getDimension(1), input.getDimension(2), input.getDatatype(), input.getDataorder());
#ifdef DEBUG
	  std::cout << "containers build" << std::endl;
#endif
	  convolve_x(&input, tmp_x, mask);
	  convolve_y(tmp_x, &result, mask);
#ifdef DEBUG
	  std::cout << "gaussian convolved with image" << std::endl;
#endif

	  tmp_x->clean();
	  delete tmp_x;
}

template<typename Scalar>
Image3D<Scalar>* GraphCut3D<Scalar>::smooth(Image3D<Scalar>& input, float sigma) {
	  std::vector<float> mask = make_fgauss(sigma);
#ifdef DEBUG
	  std::cout << "gaussian mask constructed" << std::endl;
#endif
	  normalizeMask(mask);
#ifdef DEBUG
	  std::cout << "gaussian mask normalised" << std::endl;
#endif

	  Image3D<Scalar> *tmp_x = new Image3D<Scalar>(input.getDimension(0), input.getDimension(1), input.getDimension(2), input.getDatatype(), input.getDataorder());
	  Image3D<Scalar> *tmp_y = new Image3D<Scalar>(input.getDimension(0), input.getDimension(1), input.getDimension(2), input.getDatatype(), input.getDataorder());
	  //Image3D<Scalar> *tmp_z = new Image3D<Scalar>(input.getDimension(0), input.getDimension(1), input.getDimension(2), input.getDatatype(), input.getDataorder());
#ifdef DEBUG
	  std::cout << "containers build" << std::endl;
#endif
	  convolve_x(&input, tmp_x, mask);
	  convolve_y(tmp_x, tmp_y, mask);
#ifdef DEBUG
	  std::cout << "gaussian convolved with image" << std::endl;
#endif

	  tmp_x->clean();
	  delete tmp_x;
	  return tmp_y;
}

template<typename Scalar>
GraphCut3D<Scalar>::GraphCut3D() : _weightMethod(SSD) {

}

template<typename Scalar>
GraphCut3D<Scalar>::~GraphCut3D() {

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
int GraphCut3D<Scalar>::segmentImageStack_5(Image3D<Scalar>& image, Image2D<short>& output, float sigma, float constantTreshold, int min_size) {
	uint width = image.getDimension(0);
	uint height = image.getDimension(1);

#ifdef DEBUG
	uint channels = image.getDimension(2);
	std::cout << "w: " << width << ", h:" << height <<", c: " << channels << std::endl;
#endif


	Image3D<Scalar>* smoothVar = smooth(image, sigma);
#ifdef DEBUG
	std::cout << "performed smoothing" << std::endl;
#endif

	// build graph
	unsigned long numVerts = width*height;
	//graphEdge *edges = new graphEdge[numVerts*(9+4+6)];
	std::vector<graphEdge> edges = std::vector<graphEdge>(numVerts*4);
	// FIND OUT THE COMBINATION OF UNIQUE EDGES!! => paper Felzenszwalb2004: undirected edges
	// alternative: add all, then remove duplicates automatically (empty cells at borders), then remove mutual duplicates.
	unsigned long numEdges = 0;
	for (uint y = 0; y < height; y++) {
		for (uint x = 0; x < width; x++) {

				if (y < height-1) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x,width, y+1,height, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x,y+1));
					numEdges++;
				}

				if (x > 0) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x-1,width, y,height, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x-1,y));
					numEdges++;
				}

				if (x < width-1) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x+1,width, y,height, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x+1,y));
					numEdges++;
				}

				if (y > 0) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x,width, y-1,height, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x,y-1));
					numEdges++;
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
	std::vector<ulong> uniqueSet = u->unique_set(width*height);
	std::vector<int> uniqueSizeSet = u->unique_sizeSet(width*height);
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


	for (uint y = 0; y < height; y++) {
		for (uint x = 0; x < width; x++) {
			//ulong comp = u->find(y * width + x);
			ulong address = getDataAddress_field2D(x,width, y,height, image.getDataorder());
//#ifdef DEBUG
//			std::cout << "search for connected component addess " << address << " ("<< x <<","<< y <<","<< z <<") ..." << std::endl;
//#endif
			ulong comp = u->find(address);
//#ifdef DEBUG
//			std::cout << "setting segment " << componentMap[comp] << " for ("<< x <<","<< y <<","<< z <<") ..." << std::endl;
//#endif
			output.setImageValue(short(componentMap[comp]), x,y);
		}
	}

	delete u;

	return num_ccs;
}

template<typename Scalar>
int GraphCut3D<Scalar>::segmentImageStack_9(Image3D<Scalar>& image, Image2D<short>& output, float sigma, float constantTreshold, int min_size) {
	uint width = image.getDimension(0);
	uint height = image.getDimension(1);

#ifdef DEBUG
	uint channels = image.getDimension(2);
	std::cout << "w: " << width << ", h:" << height <<", c: " << channels << std::endl;
#endif


	Image3D<Scalar>* smoothVar = smooth(image, sigma);
#ifdef DEBUG
	std::cout << "performed smoothing" << std::endl;
#endif

	// build graph
	unsigned long numVerts = width*height;
	//graphEdge *edges = new graphEdge[numVerts*(9+4+6)];
	std::vector<graphEdge> edges = std::vector<graphEdge>(numVerts*8);
	// FIND OUT THE COMBINATION OF UNIQUE EDGES!! => paper Felzenszwalb2004: undirected edges
	// WEIGHT THE WEIGHTS DEPENDING ON THE KERNEL DISTANCE!!! (1-distance) rather -> higher weights mean higher difference
	// alternative: add all, then remove duplicates automatically (empty cells at borders), then remove mutual duplicates.
	unsigned long numEdges = 0;
	for (uint y = 0; y < height; y++) {
		for (uint x = 0; x < width; x++) {

				if ((x > 0) && (y < height-1)) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x-1,width, y+1,height, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x-1,y+1));
					numEdges++;
				}

				if (y < height-1) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x,width, y+1,height, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x,y+1));
					numEdges++;
				}

				if ((x < width-1) && (y < height-1)) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x+1,width, y+1,height, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x+1,y+1));
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if (x > 0) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x-1,width, y,height, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x-1,y));
					numEdges++;
				}

				if (x < width-1) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x+1,width, y,height, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x+1,y));
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if ((x > 0) && (y > 0)) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x-1,width, y-1,height, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x-1,y-1));
					numEdges++;
				}

				if (y > 0) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x,width, y-1,height, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x,y-1));
					numEdges++;
				}

				if ((x < width-1) && (y > 0)) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x+1,width, y-1,height, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x+1,y-1));
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


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
	std::vector<ulong> uniqueSet = u->unique_set(width*height);
	for(size_t i=0; i<uniqueSet.size(); i++) {
		componentMap[uniqueSet[i]] = i;
	}
#ifdef DEBUG
	std::cout << "constructed component map - " << uniqueSet.size() << " unique components." << std::endl;
#endif
	//image_multichannel<short> *output = new image_multichannel<short>(width, height, 1, true);

	for (uint y = 0; y < height; y++) {
		for (uint x = 0; x < width; x++) {
			//ulong comp = u->find(y * width + x);
			ulong address = getDataAddress_field2D(x,width, y,height, image.getDataorder());
//#ifdef DEBUG
//			std::cout << "search for connected component addess " << address << " ("<< x <<","<< y <<","<< z <<") ..." << std::endl;
//#endif
			ulong comp = u->find(address);
//#ifdef DEBUG
//			std::cout << "setting segment " << componentMap[comp] << " for ("<< x <<","<< y <<","<< z <<") ..." << std::endl;
//#endif
			output.setImageValue(short(componentMap[comp]), x,y);
		}
	}

	delete u;

	return num_ccs;
}

template<typename Scalar>
int GraphCut3D<Scalar>::segmentImageStack_9weighted(Image3D<Scalar>& image, Image2D<short>& output, float sigma, float constantTreshold, int min_size) {
	uint width = image.getDimension(0);
	uint height = image.getDimension(1);

#ifdef DEBUG
	uint channels = image.getDimension(2);
	std::cout << "w: " << width << ", h:" << height <<", c: " << channels << std::endl;
#endif


	Image3D<Scalar>* smoothVar = smooth(image, sigma);
	//Image3D<Scalar>* smoothVar = new Image3D<Scalar>(image.getDimension(0),image.getDimension(1),image.getDimension(2),image.getDimension(3), image.getDatatype(), image.getDataorder());
	//smoothVar->setData(image.getData(), true);
#ifdef DEBUG
	std::cout << "performed smoothing" << std::endl;
#endif

	// build graph
	unsigned long numVerts = width*height;
	//graphEdge *edges = new graphEdge[numVerts*(9+4+6)];
	std::vector<graphEdge> edges = std::vector<graphEdge>(numVerts*8);
	// FIND OUT THE COMBINATION OF UNIQUE EDGES!! => paper Felzenszwalb2004: undirected edges
	// WEIGHT THE WEIGHTS DEPENDING ON THE KERNEL DISTANCE!!! (1-distance) rather -> higher weights mean higher difference
	// alternative: add all, then remove duplicates automatically (empty cells at borders), then remove mutual duplicates.
	unsigned long numEdges = 0;
	for (uint y = 0; y < height; y++) {
		for (uint x = 0; x < width; x++) {

				if ((x > 0) && (y < height-1)) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x-1,width, y+1,height, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x-1,y+1))*(1.0f-0.18f);
					numEdges++;
				}

				if (y < height-1) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x,width, y+1,height, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x,y+1))*(1.0f-0.92f);
					numEdges++;
				}

				if ((x < width-1) && (y < height-1)) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x+1,width, y+1,height, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x+1,y+1))*(1.0f-0.18f);
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if (x > 0) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x-1,width, y,height, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x-1,y))*(1.0f-0.92f);
					numEdges++;
				}

				if (x < width-1) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x+1,width, y,height, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x+1,y))*(1.0f-0.92f);
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
				if ((x > 0) && (y > 0)) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());// y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x-1,width, y-1,height, smoothVar->getDataorder());//y * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x-1,y-1))*(1.0f-0.18f);
					numEdges++;
				}

				if (y > 0) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x,width, y-1,height, smoothVar->getDataorder());//(y+1) * width + x;
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x,y-1))*(1.0f-0.92f);
					numEdges++;
				}

				if ((x < width-1) && (y > 0)) {
					edges[numEdges].a = getDataAddress_field2D(x,width, y,height, smoothVar->getDataorder());//y * width + x;
					edges[numEdges].b = getDataAddress_field2D(x+1,width, y-1,height, smoothVar->getDataorder());//(y+1) * width + (x+1);
					edges[numEdges].w = (float)(SSD_imageStack(*smoothVar, x,y, x+1,y-1))*(1.0f-0.18f);
					numEdges++;
				}
				// ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

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
	std::vector<ulong> uniqueSet = u->unique_set(width*height);
	std::vector<int> uniqueSizeSet = u->unique_sizeSet(width*height);
	for(size_t i=0; i<uniqueSet.size(); i++) {
		componentMap[uniqueSet[i]] = i;
		sizeMap[uniqueSet[i]] = uniqueSizeSet[i];
	}

#ifdef DEBUG
	std::cout << "constructed component map - " << uniqueSet.size() << " unique components." << std::endl;
#endif
	//image_multichannel<short> *output = new image_multichannel<short>(width, height, 1, true);


	for (uint y = 0; y < height; y++) {
		for (uint x = 0; x < width; x++) {
			ulong address = getDataAddress_field2D(x,width, y,height, image.getDataorder());
//#ifdef DEBUG
//			std::cout << "search for connected component addess " << address << " ("<< x <<","<< y <<","<< z <<") ..." << std::endl;
//#endif
			ulong comp = u->find(address);
//#ifdef DEBUG
//			std::cout << "setting segment " << componentMap[comp] << " for ("<< x <<","<< y <<","<< z <<") ..." << std::endl;
//#endif
				output.setImageValue(short(componentMap[comp]), x,y);
		}
	}

	delete u;

	return num_ccs;
}

template class GraphCut3D<char>;
template class GraphCut3D<unsigned char>;
template class GraphCut3D<short>;
template class GraphCut3D<unsigned short>;
template class GraphCut3D<int>;
template class GraphCut3D<unsigned int>;
template class GraphCut3D<long>;
template class GraphCut3D<unsigned long>;
template class GraphCut3D<float>;
template class GraphCut3D<double>;




