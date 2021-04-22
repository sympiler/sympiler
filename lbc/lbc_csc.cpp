//
// Created by kazem on 9/18/19.
//

#ifndef CROSSKERNEL_LBC_CSC_H
#define CROSSKERNEL_LBC_CSC_H

#include <sparse_utilities.h>
#include <sparse_inspector.h>
#include <lbc.h>
#include <sparse_io.h>
#include "includes/lbc_utils.h"
double elap, elap_trans;
namespace sym_lib {
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
                                 double *nodeCost) {
  auto *node2partition = new int[n];
  auto *outCost = new double[n];
  auto *newOutCost = new double[n];
  auto *isChanged = new bool[n]();
  auto *visited = new bool[n]();
  auto *isMarked = new int[n]();
  auto *xi = new int[2 * n];
  int curNumOfPart = 0;
  std::vector<int> remainingNodes, remainingTmp;
  std::vector<std::vector<int>> newLeveledParList, mergedLeveledParList;
  int originalHeight = 0;
  finaLevelPtr = new int[n + 1];
  finalPartPtr = new int[n]();
  finalNodePtr = new int[n];

  finaLevelPtr[0] = 0;
  for (int i = 0; i < n; ++i) {
   node2partition[i] = -1;
   outCost[i] = 0.0;
   newOutCost[i] = 0.0;
  }
  int averageCC = 0;
  //making levelset
  int *levelPtr, *levelSet;;
  auto *node2Level = new int[n];
//  int levelNo = build_levelSet_CSC(n, lC, lR, levelPtr, levelSet);
  int levelNo = level_set_multi_graphs(n, lC, lR, levelPtr, levelSet, 1);

  //COMPUTING NODE2lEVEL
  for (int i = 0; i < levelNo; ++i) {
   for (int j = levelPtr[i]; j < levelPtr[i + 1]; ++j) {
    int node = levelSet[j];
    node2Level[node] = i;
   }
  }

#if 0
  for (int i = 0; i < levelNo; ++i) {
   std::cout<<i<<"::";
   for (int j = levelPtr[i]; j < levelPtr[i+1]; ++j) {
    std::cout<<levelSet[j]<<";";
   }
   std::cout<<"\n";
  }
#endif
  //Filling indegree array with incoming edges
  auto *inDegree = new int[n]();
  for (int s = 0; s < n; ++s) {
   for (int r = lC[s]; r < lC[s + 1]; ++r) {
    int cn = lR[r];
    inDegree[cn]++;
   }
  }
  /** l-partitioning of height **/
  int *partition2Level = new int[levelNo + 1]();
  std::vector<int> innerPartsSize;
  originalHeight = levelNo;
  std::vector<std::vector<int>> slackGroups(originalHeight + 1);
  std::vector<std::vector<int>> slackedLevelSet(originalHeight + 1);
  int lClusterCnt = height_partitioning_DAG_trng(levelNo, levelPtr,
                                                 nullptr, originalHeight, innerParts,
                                                 minLevelDist, divRate,
                                                 innerPartsSize, slackGroups,
                                                 nullptr, partition2Level);

  /** w-partitioning of each chunk of level **/
  for (int l = 0; l < lClusterCnt; ++l) {//for each leveled partition
   int lbLevel = partition2Level[l] - 1;
   int ubLevel = partition2Level[l + 1];
   int dfsLevel = partition2Level[l];
   int curLeveledParCost = 0;
   //Marking lower bound

   // FIXME: we might need to do all levels below for general DAG
   for (int j = levelPtr[lbLevel > 0 ? lbLevel : 0];
        j < levelPtr[lbLevel + 1]; ++j) {
    int curNode = levelSet[j];
    isMarked[curNode] = true;
   }
   //Marking upper bound
   for (int ii = ubLevel; ii < originalHeight; ++ii) {
    for (int j = levelPtr[ii]; j < levelPtr[ii + 1]; ++j) {
     int curNode = levelSet[j];
     isMarked[curNode] = true;
    }
   }
   //Iterating over non-visited leaf nodes_ to compute CCs
   //CC:connected component
   int stackStart = 0, cc = 0, max_cc=0;
   std::vector<int> needAliased;
   bool *isUniq = new bool[n]();
   for (int k = levelPtr[dfsLevel]; k < levelPtr[dfsLevel + 1]; ++k) {
    int curLeaf = levelSet[k];
    bool isCC = true;
    int minAliasedPar = INT32_MAX;
    if (!isMarked[curLeaf]) {
     /**
      * find clusters within the vertical section. Label overlapped nodes and put them into partition
      * with minimum number of nodes
      **/
     stackStart = dfs_CSC_CC(n, curLeaf, lC, lR, isMarked, n, xi, xi + n, needAliased, nullptr);
     //Finding unique clusters from needAliased
     minAliasedPar = make_unique(node2partition, needAliased, n, isUniq);

     isCC = needAliased.empty();
     if (!isCC) {//There are some intersection between found CCs.
      for (int j : needAliased) {//the first is min
       int tn = node2partition[j];
       if (tn != minAliasedPar) {
       // cc--;
        for (int i = 0; i < n; ++i) {
         // Replace all needAliased node with their min part number.
         if (node2partition[i] == tn) {
          node2partition[i] = minAliasedPar;
         }
        }
       }
      }
      needAliased.clear();

      //Set the nodes_ in the current cluster
      for (int i = stackStart; i < n; ++i) {
       int node = xi[i];
       node2partition[node] = minAliasedPar;
       max_cc = std::max(max_cc, minAliasedPar);
       outCost[minAliasedPar] += nodeCost[node];
       curLeveledParCost += nodeCost[node];
       if (node2Level[node] != dfsLevel)
        isMarked[node] = -1;
      }
     } else {
      // reset all nodes_ but leaves to -1
      for (int i = stackStart; i < n; ++i) {
       int node = xi[i];
       node2partition[node] = cc;
       outCost[cc] += nodeCost[node];
       curLeveledParCost += nodeCost[node];
       if (node2Level[node] != dfsLevel)
        isMarked[node] = -1;
      }
      cc++;// one more CC.
     }
    }
   }
   //Reset all marked node in the DAG
   //std::cout<<cc<<"\n";
   int lb = lbLevel > 0 ? lbLevel : 0;
   for (int j = levelPtr[lb]; j < levelPtr[lb + 1]; ++j) {
    int curNode = levelSet[j];
    isMarked[curNode] = false;
   }
   delete[]isUniq;
   //Marking upper bound
   for (int ii = ubLevel; ii < originalHeight; ++ii) {
    for (int j = levelPtr[ii]; j < levelPtr[ii + 1]; ++j) {
     int curNode = levelSet[j];
     isMarked[curNode] = false;
     visited[curNode] = true;//Make it ready for mod-BFS
    }
   }
   //Topological sort of each cc, the fastest way, FIXME: make it more
   //local
   std::vector<int> extraDim;
   for (int i = 0; i < cc; ++i) {
    newLeveledParList.push_back(extraDim);
   }

   modified_BFS_CSC(n, lC, lR, inDegree, visited,
                    node2partition, levelPtr, levelSet, dfsLevel,
                    newLeveledParList);
   //Marking upper bound
   for (int ii = ubLevel; ii < originalHeight; ++ii) {
    for (int j = levelPtr[ii]; j < levelPtr[ii + 1]; ++j) {
     int curNode = levelSet[j];
     visited[curNode] = false;//Make it ready for mod-BFS
    }
   }
   //Bin packing and form W-partitions
   int levelParCostThresh = curLeveledParCost / innerParts;
   levelParCostThresh += int(0.1 * levelParCostThresh);
   int outinnerParts = 0;
   averageCC += newLeveledParList.size();
   mergedLeveledParList.resize(innerPartsSize[l]);//FIXME;
   if (newLeveledParList.size() > innerPartsSize[l]) {
    outinnerParts = worst_fit_bin_pack(newLeveledParList, outCost,
                                       mergedLeveledParList, newOutCost,
                                       levelParCostThresh, innerPartsSize[l]);
    //assert(outinnerParts<=innerParts);
   } else {
    mergedLeveledParList.clear();
    mergedLeveledParList = newLeveledParList;
    outinnerParts = newLeveledParList.size();
#if 0
    if(outinnerParts>1) {
     for (int ii = 0; ii < outinnerParts; ++ii) {
      std::cout << outCost[ii] << ";";
     }
     for (int ii = outinnerParts; ii < innerParts; ++ii) {
      std::cout << "0;";
     }
    }
#endif
   }

   double curPartCost = 0;
   finaLevelPtr[l + 1] = finaLevelPtr[l] + outinnerParts;

   for (int i = 0; i < outinnerParts; ++i) {
    int curPartElem = 0;
    curPartCost = 0;
    for (int ii=0; ii < mergedLeveledParList[i].size(); ii++) {
     int j = mergedLeveledParList[i][ii];
     curPartCost += nodeCost[j];
     finalNodePtr[finalPartPtr[curNumOfPart] + curPartElem] = j;
     node2partition[j] = curNumOfPart;
     curPartElem++;
    }
#if 0
    //std::cout<<"parts: "<<newLeveledParList.size()<<","<<curPartCost<<", ";
    std::cout<<curPartCost<<", ";
#endif
    finalPartPtr[curNumOfPart + 1] = finalPartPtr[curNumOfPart] + curPartElem;
    curNumOfPart++;
   }

   //Cleaning the current sets.
   for (auto &i : mergedLeveledParList)
    i.clear();
   mergedLeveledParList.clear();
   for (auto &i : newLeveledParList)
    i.clear();
   newLeveledParList.clear();

  }

  finaLevelNo = lClusterCnt;
#ifdef CHECK
  if (true) {//Verification of the set.
      bool *checkExist = new bool[n];
      for (int i = 0; i < n; ++i) checkExist[i] = false;
      for (int i = 0; i < lClusterCnt; ++i) {
          for (int k = finaLevelPtr[i]; k < finaLevelPtr[i + 1]; ++k) {
              for (int j = finalPartPtr[k]; j < finalPartPtr[k + 1]; ++j) {
                  assert(checkExist[finalNodePtr[j]] == false);
                  checkExist[finalNodePtr[j]] = true;
              }
          }
      }
      for (int i = 0; i < n; ++i) {
          assert(checkExist[i] == true);
      }
      delete[] checkExist;
  }
#endif
  delete[] levelPtr;
  delete[] levelSet;
  delete[] outCost;
  delete[] newOutCost;
  delete[] partition2Level;
  delete[] node2partition;
  delete[] node2Level;
  delete[] isChanged;
  delete[] isMarked;
  delete[]xi;
  delete[]visited;
  delete[]inDegree;

  return averageCC / lClusterCnt;
 }


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
                                 double *nodeCost){
  int ret = 0;
  CSC *A = new CSC(n,n,lC[n],lC,lR, stype);

  // make it symmetric and upper triangular
  timing_measurement t1, t2;
  t1.start_timer();
  int *etree = new int[n]();
  if(stype != 1){
   CSC* A_sym = make_symmetric(A, false);
   t1.start_timer();
   elap_trans = t1.measure_elapsed_time();
   t2.start_timer();
   ret = compute_etree(A_sym, etree);
   t2.start_timer();
  elap = t2.measure_elapsed_time();
   delete A_sym;
  } else{
   ret = compute_etree(A, etree);
  }
  //std::cout<<t.elapsed_time<<"\n";
#if 0
  print_csc(1, "merged \n", A);
  print_vec("ETREE\n",0,n,etree);
  CSC *tmp = tree_to_csc(n, etree);
  print_csc(1,"Tree\n",n,tmp->p, tmp->i, NULLPNTR);
  delete tmp;
#endif
  //timing_measurement t1;
  //t1.start_timer();
  ret = get_coarse_levelSet_tree(n,etree,finaLevelNo,finaLevelPtr,partNo,
    finalPartPtr, finalNodePtr, innerParts, minLevelDist, divRate,
     nodeCost);
  //t1.measure_elapsed_time();
  //std::cout<<"--> "<<t1.elapsed_time<<"\n";
  delete []etree;
  delete A;
  return ret;
 }
}
#endif //CROSSKERNEL_LBC_CSC_H
