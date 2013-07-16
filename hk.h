#ifndef HK_H_
#define HK_H_

int* extended_hoshen_kopelman(const boost::multi_array<int, 2>& nbs,
															int* occupancy, int N);
int hoshen_kopelman(int **matrix, int m, int n);

#endif /* HK_H_ */
