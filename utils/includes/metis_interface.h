//
// Created by kazem on 10/13/19.
//

#ifndef FUSION_METIS_INTERFACE_H
#define FUSION_METIS_INTERFACE_H

#include "sparse_utilities.h"

namespace sym_lib{

 /// Returns the partition done by METIS
 /// \param A matrix A to partition
 /// \param npar number of partitions to make
 /// \param partition output array of node to partition mapping
 /// \return the permuted matrix if any
 CSC *metis_partitioning(CSC *A, int npar, int *&partition);


 /// Computes the Nested-dissection permutation of the input matrix A
 /// \param A symmetric matrix with only half stored
 /// \param perm the output permutation array
 /// \return operation status set by METIS
 int metis_perm_symmetric(CSC *A, int *&perm);


 /// Computes the Nested-dissection permutation of the input matrix A
 /// \param A general matrix
 /// \param perm the output permutation array
 /// \return operation status set by METIS
 int metis_perm_general(CSC *A, int *&perm);


 /// k-part partitioning of symmetric matrix A
 /// \param A symmetric matrix with only half stored
 /// \param part the output partitioning array
 /// \param k number of partitions to do
 /// \return operation status set by METIS
 int metis_partition_symmetric(CSC *A, int *&part, int k);


 /// k-part partitioning of general matrix A
 /// \param A general matrix
 /// \param part the output partitioning array
 /// \param k number of partitions to do
 /// \return operation status set by METIS
 int metis_partition_general(CSC *A, int *&part, int k);


 /// First coarsens every coarsen rows/cols and then partition it with metis
 /// \param A in can be symmetric or general
 /// \param part  output
 /// \param k in
 /// \param coarsen  in
 /// \return
 int metis_partition_coarsened(CSC *A, int *&part, int k, int coarsen);
}


#endif //FUSION_METIS_INTERFACE_H
