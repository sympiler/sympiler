//
// Created by Kazem on 10/14/19.
//

#ifndef PROJECT_LBC_H
#define PROJECT_LBC_H

#include <cstddef>

namespace sym_lib{

 int get_coarse_levelSet_DAG_CSC(size_t n,
                                 int *lC,
                                 int *lR,
                                 int &finaLevelNo,
                                 int *&finaLevelPtr,
                                 int &partNo,
                                 int *&finalPartPtr,
                                 int *&finalNodePtr,
                                 int innerParts,
                                 int minLevelDist,
                                 int divRate,
                                 double *nodeCost);

 /// Comnputes coarsened level-set for CSC DAG using etree
 /// \param n
 /// \param lC
 /// \param lR
 /// \param finaLevelNo
 /// \param finaLevelPtr
 /// \param partNo
 /// \param finalPartPtr
 /// \param finalNodePtr
 /// \param innerParts
 /// \param minLevelDist
 /// \param divRate
 /// \param nodeCost
 /// \return
 int get_coarse_levelSet_DAG_CSC_tree(size_t n,
                                 int *lC,
                                 int *lR,
                                 int stype,
                                 int &finaLevelNo,
                                 int *&finaLevelPtr,
                                 int &partNo,
                                 int *&finalPartPtr,
                                 int *&finalNodePtr,
                                 int innerParts,
                                 int minLevelDist,
                                 int divRate,
                                 double *nodeCost);

 int getCoarseLevelSet_DAG_BCSC02(size_t n,
                                  size_t *lC,
                                  size_t *Li_ptr,
                                  int* lR,
                                  const int *blk2col,
                                  const int *col2blk,
                                  int &finaLevelNo,
                                  int* &finaLevelPtr,
                                  int* &parLevelSet,
                                  int &partNo,
                                  int* &finalPartPtr,
                                  int* &finalNodePtr,
                                  int innerParts,
                                  int minLevelDist,
                                  int divRate,
                                  double *nodeCost);

 int get_coarse_levelSet_tree(size_t n,
                              const int *eTree,
                              int &finaLevelNo,
                              int* &finaLevelPtr,
                              int &partNo,
                              int* &finalPartPtr,
                              int* &finalNodePtr,
                              int innerParts,
                              int minLevelDist,
                              int divRate,
                              double *nodeCost);
}

#endif //PROJECT_LBC_H
