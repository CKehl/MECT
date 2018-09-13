/*
 * graphDisjointForest.cpp
 *
 *  Created on: Apr 5, 2018
 *      Author: christian
 */

#include "graphDisjointForest.hpp"
#include <iostream>

// threshold function
#define THRESHOLD(size, c) (c/size)

disjointForest_implicitGrid::disjointForest_implicitGrid() : elts(NULL), num(0) {
	// TODO Auto-generated constructor stub

}

disjointForest_implicitGrid::disjointForest_implicitGrid(unsigned long elements) {
	// TODO Auto-generated constructor stub
	elts = new uni_elt[elements];
	num = elements;
	for (ulong i = 0; i < elements; i++) {
		elts[i].rank = 0;
		elts[i].size = 1;
		elts[i].p = i;
	}
}

disjointForest_implicitGrid::~disjointForest_implicitGrid() {
	// TODO Auto-generated destructor stub
	delete [] elts;
}

void disjointForest_implicitGrid::setElementNumber(unsigned long elements) {
	elts = new uni_elt[elements];
	num = elements;
	for (ulong i = 0; i < elements; i++) {
		elts[i].rank = 0;
		elts[i].size = 1;
		elts[i].p = i;
	}
}

unsigned long disjointForest_implicitGrid::find(unsigned long x) {
	unsigned long y = x;
	while (y != elts[y].p)
		y = elts[y].p;
	elts[x].p = y;
	return y;
}

void disjointForest_implicitGrid::join(unsigned long x, unsigned long y) {
	if (elts[x].rank > elts[y].rank) {
		elts[y].p = x;
		elts[x].size += elts[y].size;
	} else {
		elts[x].p = y;
		elts[y].size += elts[x].size;
		if (elts[x].rank == elts[y].rank)
			elts[y].rank++;
	}
	num--;
}

disjointForest_implicitGrid* segmentGraph_implicitGrid(unsigned long num_vertices, unsigned long num_edges, std::vector<graphEdge>& edges, float constantTreshold) {
	  // sort edges by weight
	  std::sort(edges.begin(), edges.end());

	  // make a disjoint-set forest
	  disjointForest_implicitGrid *u = new disjointForest_implicitGrid(num_vertices);
//#ifdef DEBUG
//	  std::cout << "created disjoint set." << std::endl;
//#endif

	  // init thresholds
	  float *threshold = new float[num_vertices];
//#ifdef DEBUG
//	  std::cout << "created threshold array." << std::endl;
//#endif
	  for (ulong i = 0; i < num_vertices; i++)
	    threshold[i] = THRESHOLD(1,constantTreshold);

	  // for each edge, in non-decreasing weight order...
	  for (ulong i = 0; i < edges.size(); i++) {
	    graphEdge *pedge = &(edges[i]);

	    // components connected by this edge
//#ifdef DEBUG
//	    std::cout << "search for a ("<< pedge->a <<") in edge " << i << " ... ";
//#endif
	    int a = u->find(pedge->a);
//#ifdef DEBUG
//	    std::cout << "found " << a << "." << std::endl;
//	    std::cout << "search for b ("<< pedge->b <<") in edge " << i << " ... ";
//#endif
	    int b = u->find(pedge->b);
//#ifdef DEBUG
//	    std::cout << "found " << b << "." << std::endl;
//#endif
	    if (a != b) {
	      if ((pedge->w <= threshold[a]) && (pedge->w <= threshold[b])) {
	    	  u->join(a, b);
	    	  a = u->find(a);
	    	  threshold[a] = pedge->w + THRESHOLD(u->size(a), constantTreshold);
	      }
	    }
	  }

	  // free up
	  delete threshold;
	  return u;
}

disjointForest_implicitGrid* segmentGraph_implicitGrid(unsigned long num_vertices, unsigned long num_edges, graphEdge *edges, float constantTreshold) {
	  // sort edges by weight
	  std::sort(edges, edges + num_edges);

	  // make a disjoint-set forest
	  disjointForest_implicitGrid *u = new disjointForest_implicitGrid(num_vertices);
#ifdef DEBUG
	  std::cout << "created disjoint set." << std::endl;
#endif

	  // init thresholds
	  float *threshold = new float[num_vertices];
#ifdef DEBUG
	  std::cout << "created threshold array." << std::endl;
#endif
	  for (ulong i = 0; i < num_vertices; i++)
	    threshold[i] = THRESHOLD(1,constantTreshold);

	  // for each edge, in non-decreasing weight order...
	  for (ulong i = 0; i < num_edges; i++) {
	    graphEdge *pedge = &edges[i];

	    // components conected by this edge
	    int a = u->find(pedge->a);
	    int b = u->find(pedge->b);
	    if (a != b) {
	      if ((pedge->w <= threshold[a]) &&
		  (pedge->w <= threshold[b])) {
		u->join(a, b);
		a = u->find(a);
		threshold[a] = pedge->w + THRESHOLD(u->size(a), constantTreshold);
	      }
	    }
	  }

	  // free up
	  delete threshold;
	  return u;
}
