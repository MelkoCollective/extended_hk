//TODO: Replace with unit-testing structure later.

#include "hk.h"
#include "MersenneTwister.h"
#include <cstdio>

using namespace std;

/* Script for the boost array case -----------------------*/

int main(int argc, char *argv[]) {
  int L = atoi(argv[1]);   //The length of the lattice (ex: L =1 gives 4 nodes)
  double p = atof(argv[2]); //The probability of populating a node.
  bool verbose = !!atoi(argv[3]);
  const int N = L*L;

  // Typedefs
  typedef boost::multi_array<int, 1> array_1t;
  typedef boost::multi_array<int, 2> array_2t;

  // Build a square lattice network with periodic lattice conditions.
  boost::multi_array<int, 2> nbs(boost::extents[N][4]);
  for (int k = 0; k < L*L; ++k) {
    nbs[k][0] = (k/L)*L + (k%L+1)%L; //top
    nbs[k][1] = ((k/L+1)%L)*L + k%L; //right
    nbs[k][2] = (k/L)*L + (k%L+L-1)%L; //bottom
    nbs[k][3] = ((k/L+L-1)%L)*L+k%L; //left
  }

  // Build the occupancies of the sites.
  MTRand mrand = MTRand();
  //int* occupancy = new int[N];
  array_1t occupancy(boost::extents[N]);
  for (int i=0; i<N; ++i) {
    double number = mrand();
    if (number < p)
      occupancy[i] = 1;
    else
      occupancy[i] = 0;
  }

  // Print the matrix to be input into HK.
  if (verbose) {
    cout << "---INPUT---" << endl << endl;
    for (int i = 0; i < N; i++){
      if (i%L == 0)
        cout << endl;
      printf("%3d ", occupancy[i]);
    }
  }

  boost::multi_array<int, 1> node_labels;
  extended_hoshen_kopelman(node_labels,nbs,occupancy);

  // Print the output matrix.
  if (verbose) {
    cout << "---OUTPUT---" << endl << endl;
    for (int i = 0; i < N; i++){
      if (i%L == 0)
        cout << std::endl;
      printf("%3d ", node_labels[i]);
    }
  }

  return 0;
}

/* Script for the no boost function --------------------------*/

//int main(int argc, char *argv[]) {
//  int L = atoi(argv[1]);   //The length of the lattice (ex: L =1 gives 4 nodes)
//  double p = atof(argv[2]); //The probability of populating a node.
//  bool verbose = !!atoi(argv[3]);
//  const int N = L*L;
//
//  // Build a square lattice network with periodic lattice conditions.
//  int ** nbs = new int*[N];
//  for (int i = 0; i < N; ++i) {
//    nbs[i] = new int[4];
//  }
//  for (int k = 0; k < L*L; ++k) {
//    nbs[k][0] = (k/L)*L + (k%L+1)%L; //top
//    nbs[k][1] = ((k/L+1)%L)*L + k%L; //right
//    nbs[k][2] = (k/L)*L + (k%L+L-1)%L; //bottom
//    nbs[k][3] = ((k/L+L-1)%L)*L+k%L; //left
//  }
//
//  // Build the occupancies of the sites.
//  MTRand mrand = MTRand();
//  int* occupancy = new int[N];
//  for (int i=0; i<N; ++i) {
//    double number = mrand();
//    if (number < p)
//      occupancy[i] = 1;
//    else
//      occupancy[i] = 0;
//  }
////  int myarr[] = {1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0};
////  occupancy = myarr;
//
//  // Print the matrix to be input into HK.
//  if (verbose) {
//    cout << "---INPUT---" << endl << endl;
//    for (int i = 0; i < N; i++){
//      if (i%L == 0)
//        cout << endl;
//      printf("%3d ", occupancy[i]);
//    }
//    cout << endl << endl << endl;
//  }
//
//  int node_labels [N];
//  extended_hk_no_boost(node_labels,nbs,occupancy,N,4);
//
//  // Print the output matrix.
//  if (verbose) {
//    cout << "---OUTPUT---" << endl << endl;
//    for (int i = 0; i < N; i++){
//      if (i%L == 0)
//        cout << std::endl;
//      printf("%3d ", node_labels[i]);
//    }
//  }
//
//  delete [] occupancy;
//  for (int i = 0; i < N; ++i) {
//    delete [] nbs[i];
//  }
//  delete [] nbs;
//  return 0;
//}
