#ifndef HK_H_
#define HK_H_

#include <boost/multi_array.hpp>

void extended_hoshen_kopelman(boost::multi_array<int, 1>& node_labels,
                              const boost::multi_array<int, 2>& nbs,
                               const boost::multi_array<int, 1>& occupancy);

void extended_hk_no_boost(int* node_labels, int const* const* nbs,
                          const int* occupancy, int N, int m);

int hoshen_kopelman(int **matrix, int m, int n);

#endif /* HK_H_ */
