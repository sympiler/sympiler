//
// Created by kazem on 6/13/20.
//
namespace sym_lib{

/*
 * Scale a vector by another vector vec2 = vec2/vec1
 */
 void scale_vec_vec(int n, double *vec1, double *vec2){
  for (int i = 0; i < n; ++i) {
   double tmp = *vec2 / *(vec1++);
   *(vec2++) = tmp;
   //*(vec2++) = *vec2 / *(vec1++);
  }
 }


}
