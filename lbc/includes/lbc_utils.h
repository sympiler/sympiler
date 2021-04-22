//
// Created by Kazem on 10/14/19.
//

/// Expert routines that are used in LBC algorithm.

#ifndef PROJECT_LBC_UTILS_H
#define PROJECT_LBC_UTILS_H

#include <vector>
namespace sym_lib{

 struct subTree {
  double cost;
  std::vector<int> nodeList;
 };

 ///
 /// \param node2Par
 /// \param list
 /// \param n
 /// \param ws
 /// \return
 int make_unique(int *node2Par, std::vector<int> &list, int n, bool *ws);


 /// depth-first-search of the graph of a matrix, starting at node j
 /// modified to report the intersections
 /// \param n
 /// \param j
 /// \param Gp
 /// \param Gi
 /// \param marked marked[i]=0; not visited marked[i]=1; visited now
 ///                 marked[i]=-1; visited in previous CC
 /// \param top
 /// \param xi
 /// \param pstack
 /// \param clashedNodes has the list of nodes clashed with visited nodes
 ///                       in previous CCs
 /// \param pinv
 /// \return
 int dfs_CSC_CC (size_t  n, int j, int *Gp, int *Gi, int *marked, int top,
                 int *xi, int *pstack, std::vector<int> &clashedNodes,
                 const int *pinv);


 /// depth-first-search of the graph of a matrix, starting at node j
 /// modified to report the intersections
 /// \param n
 /// \param j
 /// \param Gp
 /// \param Gi_ptr
 /// \param Gi
 /// \param blk2Col
 /// \param col2Blk
 /// \param marked marked[i]=0; not visited, marked[i]=1; visited now,
 ///         marked[i]=-1; visited in previous CC
 /// \param top
 /// \param xi
 /// \param pstack
 /// \param clashedNodes
 /// \param pinv
 /// \return
 int dfs_BCSC_CC (size_t  n, int j, size_t *Gp, size_t *Gi_ptr, int *Gi,
                  const int *blk2Col,const int *col2Blk, int *marked,
                  int top, int *xi, int *pstack,
                  std::vector<int> &clashedNodes, const int *pinv);


 ///
 /// \param newLeveledParList
 /// \param inCost
 /// \param mergedLeveledParList
 /// \param outCost
 /// \param costThreshold
 /// \return
 int merge_inner_part(std::vector<std::vector<int>> newLeveledParList,
                    int *inCost,
                    std::vector<std::vector<int>> &mergedLeveledParList,
                    int *outCost, int costThreshold);

 ///
 /// \param a
 /// \param b
 /// \return
 bool cmp_cost(subTree a, subTree b);


 ///
 /// \param cost
 /// \param size
 /// \return
 int find_min(double *cost, int size);


 ///
 /// \param newLeveledParList
 /// \param inCost
 /// \param mergedLeveledParList
 /// \param outCost
 /// \param costThreshold
 /// \param numOfBins
 /// \return
 int worst_fit_bin_pack(const std::vector<std::vector<int>> &newLeveledParList,
                     double *inCost,
                     std::vector<std::vector<int>> &mergedLeveledParList,
                     double *outCost, int costThreshold, int numOfBins);


 ///
 /// \param levelNo
 /// \param levelPtr
 /// \param node2Level
 /// \param originalHeight
 /// \param innerParts
 /// \param minLevelDist
 /// \param divRate
 /// \param innerPartsSize
 /// \param slackGroups
 /// \param subTreeCost
 /// \param partition2Level
 /// \param sw
 /// \return
 int height_partitioning_DAG_trng(int levelNo, int *levelPtr, int *node2Level,
                                 int originalHeight, int innerParts,
                                 int minLevelDist, int divRate,
                                 std::vector<int> &innerPartsSize,
                                 std::vector<std::vector<int>> &slackGroups,
                                 double *subTreeCost, int *partition2Level,
                                 bool sw = false);


 ///
 /// \param levelNo
 /// \param levelPtr
 /// \param node2Level
 /// \param originalHeight
 /// \param innerParts
 /// \param minLevelDist
 /// \param divRate
 /// \param innerPartsSize
 /// \param slackGroups
 /// \param subTreeCost
 /// \param partition2Level
 /// \param sw
 /// \return
 int height_partitioning(int levelNo, int *levelPtr, int *node2Level,
                        int originalHeight, int innerParts, int minLevelDist,
                        int divRate, std::vector<int> &innerPartsSize,
                        std::vector<std::vector<int>> &slackGroups,
                        double *subTreeCost, int *partition2Level,
                        bool sw = false);
}

#endif //PROJECT_LBC_UTILS_H
