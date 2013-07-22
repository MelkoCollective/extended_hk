#ifndef HK_H_
#define HK_H_

#include <boost/multi_array.hpp>

void extended_hoshen_kopelman(boost::multi_array<int, 1>& node_labels,
                              const boost::multi_array<int, 2>& nbs,
                              int* occupancy);
int hoshen_kopelman(int **matrix, int m, int n);

#endif /* HK_H_ */
