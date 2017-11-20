
#include <cstdlib>
#include <chrono>
#include <string>
#include <fstream>
#include <iostream>
#include "Eigen/Sparse"
#include "../../util/Util.h"

using namespace std;


template <typename valType>
void load( int row, int *colP, int *rowIdx, double *val,
           std::vector<Eigen::Triplet<valType>>& triplets){
 for (int i = 0; i < row; ++i) {
  for (int j = colP[i]; j < colP[i+1]; ++j) {
   int curRow, curCol;
   valType curVal;
   curRow=rowIdx[j]; curCol=i; curVal=val[j];
   Eigen::Triplet<valType> tmp(curRow,curCol,curVal);
   triplets.push_back(tmp);
  }
 }
}


int main(int argc, char *argv[]){
#ifdef COLDCACHE
 if(argc<3){
  cout<<"Two arguments are needed, input matrix path and a "
    "large matrix path\n";
  return -1;
 }
#else
 if(argc<2){
        printf("Please enter a path for the input matrix");
        return -1;
    }
#endif
 std::chrono::time_point<std::chrono::system_clock> start, end;
 std::chrono::duration<double> elapsed_seconds;

 std::string f1 = argv[1];
#ifdef COLDCACHE
 std::string waste=argv[2];
 ifstream wasteFile;
 wasteFile.open(waste);
 if(wasteFile.fail()){
  cout<<"The input file for simulating cold-cache conditions does "
    "not exist\n";
  return -1;
 }
#endif

 int *col, *row;
 double  *val;
 int n, nnz;

 if (!readMatrix(f1,n,nnz,col,row,val))
  return -1;

 int rown=n;
 double durationEigen=0,durationSym=0;
 std::vector<Eigen::Triplet<double>> triplets;
 double *x=new double[n];
 //Eigen
 load(rown,col,row,val,triplets);
 Eigen::SparseMatrix<double> sm1(rown,rown);
 sm1.setFromTriplets(triplets.begin(),triplets.end());//col major

 Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Lower> spdSolver;
 start = std::chrono::system_clock::now();
 spdSolver.analyzePattern(sm1);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 durationSym=elapsed_seconds.count();

 enableColdCache(1200,wasteFile);//1.2k*1.2k*8 ~ 15MB
 start = std::chrono::system_clock::now();
 spdSolver.factorize(sm1);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 durationEigen=elapsed_seconds.count();
 cout<<"symbolic time: "<<durationSym<<" Numeric time: "<<durationEigen<<";";

#if 0
 for(int i=0; i<20; ++i){
        cout<<x[i]<<",";
    }
    std::cout<<"\n";
#endif
 delete []x;
 return 1;
}
