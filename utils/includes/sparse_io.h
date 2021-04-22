//
// Created by Kazem on 10/10/19.
//

#ifndef PROJECT_SPARSE_IO_H
#define PROJECT_SPARSE_IO_H

#include <iostream>
#include <string>
#include <iomanip>
#include "def.h"
namespace sym_lib{


///
/// \param fName
/// \return
 CSC* read_mtx(std::string fName);

 ///
 /// \param fname
 /// \param A
 void CSC_to_mtx(std::string fname, CSC *A);

 ///
 /// \param fname
 /// \param A
 void BCSC_to_mtx(std::string fname, BCSC *A);

 /// Converts dense to CSC
 /// \param rows in
 /// \param cols in
 /// \param val in
 /// \return
 CSC *dense_to_csc(int rows, int cols, double **val);

 /// Converts indices to one-based
 /// \param A input CSC
 /// \return
 CSC * convert_to_one_based(const CSC *A);

 ///
/// \tparam type
/// \param header the string at the beginning of output
/// \param beg the beginning index to print
/// \param end the last index to print
/// \param vec the vector of values
 template<class type> void print_vec(std::string header,
                                     int beg,int end, type *vec){

  std::cout<<header;
  for (int i = beg; i < end; ++i) {
   std::cout<<std::setprecision(15)<<vec[i]<<", ";
  }
  std::cout<<"\n";
 }


 /// copying vec_in to vec_out from beg to end
 /// \tparam type
 /// \param beg
 /// \param end
 /// \param vec_in
 /// \param vec_out
 template<class type> void copy_vector(int beg, int end, type *vec_in, type *vec_out){
  for (int i = beg; i < end; ++i) {
   vec_out[i] = vec_in[i];
  }
 }

/// Print CSC matrix into output
/// \param beg
/// \param n
/// \param Ap
/// \param Ai
/// \param Ax
 void print_csc(int fd, std::string beg, size_t n, int *Ap, int *Ai, double *Ax);
 void print_csc(int fd, std::string beg, CSC *A);

 ///
 /// \param fd
 /// \param A
 void print_dense(int fd, Dense *A);

 /// Printing a level-set
 /// \param beg
 /// \param n
 /// \param level_ptr
 /// \param level_set
 void print_level_set(std::string beg, int n, int *level_ptr, int *level_set);

/// Prints an H-level set
/// \param beg
/// \param n
/// \param level_ptr
/// \param level_part_ptr
/// \param level_set
 void print_hlevel_set(std::string beg, int n,
                       const int *level_ptr, const int *level_part_ptr,
                       int *level_set);
}


#endif //PROJECT_SPARSE_IO_H
