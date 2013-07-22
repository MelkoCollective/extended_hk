/* Tobin Fricke's original code with an extended version of the HK algorithm
 added. The extension borrows some ideas from:

 http://gaia.pge.utexas.edu/papers/AFTWPPhysicaA.pdf

 This program has been updated to work in C++11 ("-std=c++0x" in a compatible
 compiler).

 ----------------------------------------------------------------------------
 Fricke's original file header:
 ----------------------------------------------------------------------------

 Tobin Fricke's implementation of the
 Hoshen-Kopelman algorithm for
 cluster labeling.

 Copyright (c) September 9, 2000, by Tobin Fricke <tobin@pas.rochester.edu>

 Modified 2002-03-09 Tobin Fricke
 Modified substantially 2004-04-21 by Tobin Fricke

 This program is written in the 1999 standard of the C language (C99).  Older C
 compilers will refuse to compile it.   You can use a C++ compiler, a C99 compiler,
 or you can modify this code to comply with a previous version of the C standard.
 The GCC compiler supports C99 as of version 3.0.  Compile this program with:

 gcc-3.0 -Wall -std=c99 hk.c -o hk

 http://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
 */

//#include "hk.h"
#include <boost/multi_array.hpp>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <array>
#include <new>
#include <iostream>

using namespace std;

/* Implementation of Union-Find Algorithm */

/* The 'labels' array has the meaning that labels[x] is an alias for the label x; by
 following this chain until x == labels[x], you can find the canonical name of an
 equivalence class.  The labels start at one; labels[0] is a special value indicating
 the highest label already used. */

int *labels; //TODO: consider removing from the top-of-file scope.
int n_labels = 0; /* length of the labels array */

/*  uf_find returns the canonical label for the equivalence class containing x */

int uf_find(int x) {
  int y = x;
  while (labels[y] != y)
    y = labels[y];

  while (labels[x] != x) {
    int z = labels[x];
    labels[x] = y;
    x = z;
  }
  return y;
}

/*  uf_union joins two equivalence classes and returns the canonical label of the resulting class. */

int uf_union(int x, int y) {
  return labels[uf_find(x)] = uf_find(y);
}

/*  uf_union_multi is a generalized version of uf_union for more than 2 nodes. */

int uf_union_multi(int* nodes, int N) {
  //TODO: pairwise for now. perhaps there's a quicker way? Think about this.
  //TODO: soln: just pair all up to 1 guy (equivalence is transitive duuuhh...)
  int out;
  //double-for loop across upper-triangle of cartesian product matrix
  for (int i = 0; i < N; i++)
    for (int j = i; j < N; j++)
      out = uf_union(nodes[i], nodes[j]); //TODO: redundant assgt here. prob not costly.

  return out;
}

/*  uf_make_set creates a new equivalence class and returns its label */

int uf_make_set(void) {
  labels[0]++;
  assert(labels[0] < n_labels);
  labels[labels[0]] = labels[0];
  return labels[0];
}

/*  uf_intitialize sets up the data structures needed by the union-find implementation. */

void uf_initialize(int max_labels) {
  n_labels = max_labels;
//  labels = calloc(sizeof(int), n_labels);
  labels = new int[n_labels];
  labels[0] = 0;
}

/*  uf_done frees the memory used by the union-find data structures */

void uf_done(void) {
  n_labels = 0;
//  free(labels);
  delete[] labels;
  labels = 0;
}

/* End Union-Find implementation */

#define max(a,b) (a>b?a:b)
#define min(a,b) (a>b?b:a)

/* print_matrix prints out a matrix that is set up in the "pointer to pointers" scheme
 (aka, an array of arrays); this is incompatible with C's usual representation of 2D
 arrays, but allows for 2D arrays with dimensions determined at run-time */

void print_matrix(int **matrix, int m, int n) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      printf("%3d ", matrix[i][j]);
    printf("\n");
  }
}

/* A generalized version of the HK algorithm for arbitrary networks of nodes.
 *
 * INPUT:
 * -nbs: 2d matrix. ith row is the neighbours of node i.
 * -occupancy: vector with the occupation number (0 or 1) of the nodes.
 * OUTPUT:
 * -node_labels: the labels of the nodes.
 */
/* TODO: conventionalize the above docstring.
 * TODO: Add a C matrix flavour where you have the number of neighbours
 * TODO: add a different flavour where we can assume no-look-ahead. In this
 * case, we would be able to assume labelling, and would not need the
 * node_labels entity!
 * TODO: is the redeclaration of node_nbs inefficient or does compiler handle this?
 * TODO: benefit of the type initializer way in loop, or declare then assign?
 * TODO: add check to see if already allocated space for size.
 * TODO: polish the min element part with some iteratory stuff to derive new
 * iterators over the subarrays.
 * TODO: currently using an extreme memory bound for label initialization..can
 *  wedo better?
 * TODO: perhaps switch the N+1 placeholder for 0, and make more computation
 * occur when finding min element...readability vs significant optimization.
 */
void extended_hoshen_kopelman(boost::multi_array<int, 1>& node_labels,
                              const boost::multi_array<int, 2>& nbs,
                              int* occupancy) {
  // Typedefs
  typedef boost::multi_array<int, 2> array_2t;
  typedef boost::multi_array<int, 1> array_1t;
  typedef boost::multi_array_types::index_range range_t;
  typedef boost::multi_array<int, 1>::const_iterator const_iter;

  const int N = nbs.shape()[0]; //number of nodes.

  // Initialize node_labels
  array_1t::extent_gen extents;
  node_labels.resize(extents[N]);
  fill(node_labels.begin(),node_labels.end(),N+1); // N+1 is unlabelled node

  // Initialize memory for binary forest of labels.
  uf_initialize(N);

  // Iterate over nodes and perform clustering. TODO: replace with iteration
  for (int i = 0; i < N; ++i) {
    if (occupancy[i]) {
      // Get neighbours of node i ('i'th row of neighbours)
      array_2t::index_gen indices;
      array_2t::const_array_view<1>::type node_nbs = nbs[indices[i][range_t()]];

      // Check if there exist any neighbouring labels.
      if (all_of(node_nbs.begin(), node_nbs.end(),
                 [&](int j) {return node_labels[j]==N+1;}))
        node_labels[i] = uf_make_set();
      else {
        // Get subset of labels using node_nbs as indices.
        array_1t node_nbs_labels(extents[node_nbs.shape()[0]]);
        for (unsigned int j = 0; j < node_nbs_labels.shape()[0]; ++j){
          node_nbs_labels[j] = node_labels[node_nbs[j]];
        }

        // Find smallest label of the neighbours.
        int min_label = *min_element(node_nbs_labels.begin(),
                                      node_nbs_labels.end());

        // Apply the minimum label to all labelled neighbours + current node.
        node_labels[i] = min_label;
        for (const_iter nb = node_nbs.begin(); nb != node_nbs.end(); ++nb)
          if (node_labels[*nb] != N+1)
            uf_union(min_label, node_labels[*nb]);
      } //else
    } //occupancy
  } //node

  // Recursive relabelling
  int *new_labels = new int[n_labels]; // allocate array, initialized to zero
  for (int i = 0; i < N; i++)
    if (occupancy[i]) {
      int x = uf_find(node_labels[i]);
      if (new_labels[x] == 0) {
        new_labels[0]++;
        new_labels[x] = new_labels[0];
      }
      node_labels[i] = new_labels[x];
    }
    else {
      node_labels[i] = 0; // Replace placeholders with 0.
    }

  // Cleanup
  delete[] new_labels;
  uf_done();
}

/* Label the clusters in "matrix". Return the total number of clusters found. */

int hoshen_kopelman(int **matrix, int m, int n) {

  uf_initialize(m * n / 2);

  /* scan the matrix */

  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (matrix[i][j]) {                        // if occupied ...

        int up = (i == 0 ? 0 : matrix[i - 1][j]);    //  look up
        int left = (j == 0 ? 0 : matrix[i][j - 1]);  //  look left

        switch (!!up + !!left) {

        case 0:
          matrix[i][j] = uf_make_set();      // a new cluster
          break;

        case 1:                              // part of an existing cluster
          matrix[i][j] = max(up,left);       // whichever is nonzero is labelled
          break;

        case 2:                              // this site binds two clusters
          matrix[i][j] = uf_union(up, left);
          break;
        }

      }

  /* apply the relabeling to the matrix */

  /* This is a little bit sneaky.. we create a mapping from the canonical labels
   determined by union/find into a new set of canonical labels, which are
   guaranteed to be sequential. */

  int *new_labels = new int[n_labels]; // allocate array, initialized to zero

  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (matrix[i][j]) {
        int x = uf_find(matrix[i][j]);
        if (new_labels[x] == 0) {
          new_labels[0]++;
          new_labels[x] = new_labels[0];
        }
        matrix[i][j] = new_labels[x];
      }

  int total_clusters = new_labels[0];

  delete [] new_labels;
  uf_done();

  return total_clusters;
}

/* This procedure checks to see that any occupied neighbors of an occupied site
 have the same label. */

void check_labelling(int **matrix, int m, int n) {
  int N, S, E, W;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (matrix[i][j]) {
        N = (i == 0 ? 0 : matrix[i - 1][j]);
        S = (i == m - 1 ? 0 : matrix[i + 1][j]);
        E = (j == n - 1 ? 0 : matrix[i][j + 1]);
        W = (j == 0 ? 0 : matrix[i][j - 1]);

        assert(N == 0 || matrix[i][j] == N);
        assert(S == 0 || matrix[i][j] == S);
        assert(E == 0 || matrix[i][j] == E);
        assert(W == 0 || matrix[i][j] == W);
      }
}
