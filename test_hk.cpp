//TODO: Replace with unit-testing structure later.

#include "hk.h"
#include <random>

using namespace std;

int main(int argc, char *argv[]) {
  int L = atoi(argv[1]);   //The length of the lattice (ex: L =1 gives 4 nodes)
  double p = atof(argv[2]); //The probability of populating a node.
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
  int* occupancy = new int[N];
  default_random_engine generator;
  uniform_real_distribution<double> distribution(0.0,1.0);
  for (int i=0; i<N; ++i) {
    double number = distribution(generator);
    if (number < p)
      occupancy[i] = 1;
    else
      occupancy[i] = 0;
  }

  // Print the matrix to be input into HK.
  cout << "---INPUT---" << endl << endl;
  for (int i = 0; i < N; i++){
    if (i%L == 0)
      cout << endl;
    printf("%3d ", occupancy[i]);
  }

  boost::multi_array<int, 1> node_labels;
  extended_hoshen_kopelman(node_labels,nbs,occupancy);

  // Print the output matrix.
  cout << "---OUTPUT---" << endl << endl;
  for (int i = 0; i < N; i++){
    if (i%L == 0)
      cout << std::endl;
    printf("%3d ", node_labels[i]);
  }

  return 0;
}
