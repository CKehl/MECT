/*
 * graphDisjointForest.h
 *
 *  Created on: Apr 5, 2018
 *      Author: christian
 */

#ifndef GRAPHDISJOINTFOREST_H_
#define GRAPHDISJOINTFOREST_H_

#include "graphEdge.hpp"
#include <vector>
#include <algorithm>

// disjoint-set forests using union-by-rank and path compression (sort of).
typedef struct {
  int rank;
  unsigned long p;
  unsigned long size;
} uni_elt;

class disjointForest_implicitGrid {
public:
	disjointForest_implicitGrid();
	disjointForest_implicitGrid(unsigned long elements);
	virtual ~disjointForest_implicitGrid();

	void setElementNumber(unsigned long elements);

	unsigned long find(unsigned long x);
	void join(unsigned long x, unsigned long y);
	int size(unsigned long x) const { return elts[x].size; }
	unsigned long num_sets() const { return num; }

	inline std::vector<unsigned long> unique_set(unsigned long inputRange) {
		std::vector<unsigned long> result;
		for(unsigned long i=0; i<inputRange; i++) {
			unsigned long val = elts[i].p;
			if(std::find(result.begin(), result.end(), val)==result.end())
				result.push_back(val);
		}
		return result;
	}
	inline std::vector<int> unique_sizeSet(unsigned long inputRange) {
		std::vector<int> result;
		std::vector<unsigned long> uniqueIds;
		for(unsigned long i=0; i<inputRange; i++) {
			unsigned long val = elts[i].p;
			if(std::find(uniqueIds.begin(), uniqueIds.end(), val)==uniqueIds.end()) {
				result.push_back(elts[i].size);
				uniqueIds.push_back(val);
			}
		}
		return result;
	}

private:
	uni_elt* elts;
	unsigned long num;
};

disjointForest_implicitGrid* segmentGraph_implicitGrid(unsigned long num_vertices, unsigned long num_edges, std::vector<graphEdge>& edges, float constantTreshold);
disjointForest_implicitGrid* segmentGraph_implicitGrid(unsigned long num_vertices, unsigned long num_edges, graphEdge *edges, float constantTreshold);


#endif /* GRAPHDISJOINTFOREST_H_ */
