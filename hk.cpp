/* Tobin Fricke's original code with an extended version of the HK algorithm
 added. The extension borrows some ideas from:
 http://gaia.pge.utexas.edu/papers/AFTWPPhysicaA.pdf

 Requirements:
 -Boost
 -compiler ? - ?

 Copyright (c) September 9, 2000, by Tobin Fricke <tobin@pas.rochester.edu>

 Extended 2013-08-14 Grant Watson
 Modified 2002-03-09 Tobin Fricke
 Modified substantially 2004-04-21 by Tobin Fricke

 Fricke's code available at:
 http://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
 */

#include <boost/multi_array.hpp>
#include <boost/phoenix.hpp>
#include <algorithm>
#include <cassert>
#include <new>

using namespace std;

/* Implementation of Union-Find Algorithm */

/* The 'labels' array has the meaning that labels[x] is an alias for the label x; by
 following this chain until x == labels[x], you can find the canonical name of an
 equivalence class.  The labels start at one; labels[0] is a special value indicating
 the highest label already used. */

int *labels;
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
  labels = new int[n_labels];
  labels[0] = 0;
}

/*  uf_done frees the memory used by the union-find data structures */

void uf_done(void) {
  n_labels = 0;
  delete[] labels;
  labels = 0;
}

/* End Union-Find implementation */

#define max(a,b) (a>b?a:b)
#define min(a,b) (a>b?b:a)

/* A generalized version of the HK algorithm for arbitrary networks of nodes.
 *
 * INPUT:
 * -nbs: 2d matrix. ith row is the neighbours of node i.
 * -occupancy: vector with the occupation number (0 or 1) of the nodes.
 * OUTPUT:
 * -node_labels: the labels of the nodes.
 */
/*
 * TODO: add a different flavour where we can assume no-look-ahead. In this
 * case, we would be able to assume labelling, and would not need the
 * node_labels entity!
 * TODO: is the redeclaration of node_nbs inefficient or does compiler handle this?
 * TODO: benefit of the type initializer way in loop, or declare then assign?
 * TODO: polish the min element part with some iteratory stuff to derive new
 * iterators over the subarrays.
 * TODO: Without C++11 we lost the syntactically sleek lambda function...is
 * there still a sleek way to do the all_of?  This would be useful for a
 * couple locations in code (ex: finding is_alone)
 */
void extended_hoshen_kopelman(boost::multi_array<int, 1>& node_labels,
                              const boost::multi_array<int, 2>& nbs,
                              const boost::multi_array<int, 1>& occupancy) {
  // Typedefs
  typedef boost::multi_array<int, 2> array_2t;
  typedef boost::multi_array<int, 1> array_1t;
  typedef boost::multi_array_types::index_range range_t;
  typedef boost::multi_array<int, 1>::const_iterator c_iter;

  const int N = nbs.shape()[0]; //number of nodes.

  // Initialize node_labels
  array_1t::extent_gen extents;
  node_labels.resize(extents[N]);
  // labels live in [1,N]. Use N+1 instead of 0 to reduce some computation.
  int unlabelled = N+1;
  fill(node_labels.begin(),node_labels.end(),N+1);

  // Initialize memory for binary forest of labels.
  uf_initialize(N);

  // Iterate over nodes and perform clustering. TODO: replace with iteration
  for (int i = 0; i < N; ++i) {
    if (occupancy[i]) {

      // Get neighbours of node i ('i'th row of neighbours)
      array_2t::index_gen indices;
      array_2t::const_array_view<1>::type node_nbs =
                                              nbs[ indices[i][range_t()] ];

      // Get subset of labels using node_nbs as indices.
      array_1t node_nbs_labels(extents[node_nbs.shape()[0]]);
      for (unsigned int j = 0; j < node_nbs_labels.shape()[0]; ++j) {
        node_nbs_labels[j] = node_labels[node_nbs[j]];
      }
      c_iter begin = node_nbs_labels.begin();
      c_iter end = node_nbs_labels.end();

      // Check if node has no labeled neighbours.
      bool is_alone = 1;
      for (c_iter it = begin; is_alone && (it != end); ++it){
        is_alone *= (*it == unlabelled);
      }

      // Labelling + merging
      if(is_alone)
        node_labels[i] = uf_make_set();
      else {
        // Find smallest label of the neighbours.
        int min_label = *min_element(begin, end);

        // Apply the minimum label to all labelled neighbours + current node.
        node_labels[i] = min_label;
        for (c_iter nb = node_nbs.begin(); nb != node_nbs.end(); ++nb)
          if (node_labels[*nb] != unlabelled)
            uf_union(min_label, node_labels[*nb]);
      }

    } //occupancy
  } //node

  /* This is a little bit sneaky.. we create a mapping from the canonical labels
   determined by union/find into a new set of canonical labels, which are
   guaranteed to be sequential. */

  int *new_labels = new int[n_labels]; // allocate array, initialized to zero
  for (int i = 0; i < n_labels; i++) new_labels[i] = 0;

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

/* A flavour of extended_hoshen_kopelman that uses standard arrays instead of
 * boost. Note that it requires a constant number of neighbours per site.
 *
 * INPUT:
 * -nbs: 2d matrix (N x m). ith row is the neighbours of node i.
 * -occupancy: vector with the occupation number (0 or 1) of the nodes.
 * -N: the number of nodes.
 * -m: the number of neighbours per node
 * OUTPUT:
 * -node_labels: the labels of the nodes. Assumed that space already allocated
 */

void extended_hk_no_boost(int* node_labels, int const* const* nbs,
                          const int* occupancy, int N, int m) {

  // Initialize node_labels with N+1 since labels live in [1,N].
  int unlabelled = N+1;
  for (int i = 0; i < N; ++i)
    node_labels[i] = unlabelled;

  // Initialize memory for binary forest of labels.
  uf_initialize(N);

  int node_nbs_labels [m];

  // Iterate over nodes and perform clustering.
  for (int i = 0; i < N; ++i) {
    if (occupancy[i]) {

      // Get neighbours of node i ('i'th row of neighbours)
      const int* node_nbs = nbs[i];

      // Get subset of labels using node_nbs as indices (ie node_labels[node_nbs])
      for (int j = 0; j < m; ++j) {
        node_nbs_labels[j] = node_labels[node_nbs[j]];
      }

      // Check if node has no labeled neighbours.
      bool is_alone = 1;
      for (int j = 0; is_alone && (j < m); ++j){
        is_alone *= (node_nbs_labels[j] == unlabelled);
      }

      // Labelling + merging
      if(is_alone)
        node_labels[i] = uf_make_set();
      else {
        // Find smallest label of the neighbours.
        int min_label = *min_element(&node_nbs_labels[0], &node_nbs_labels[m]);

        // Apply the minimum label to all labelled neighbours + current node.
        node_labels[i] = min_label;
        for (int j = 0; j < m; ++j)
//          if (node_labels[node_nbs[j]] != unlabelled)
//            uf_union(min_label, node_labels[node_nbs[j]]);
          if (node_nbs_labels[j] != unlabelled)
            uf_union(min_label,node_nbs_labels[j]);
      }

    } //occupancy
  } //node

  /* This is a little bit sneaky.. we create a mapping from the canonical labels
   determined by union/find into a new set of canonical labels, which are
   guaranteed to be sequential. */

  int *new_labels = new int[n_labels]; // allocate array, initialized to zero
  for (int i = 0; i < n_labels; i++) new_labels[i] = 0;

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
  for (int i = 0; i < n_labels; i++) new_labels[i] = 0;

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
