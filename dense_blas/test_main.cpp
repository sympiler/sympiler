//
// Created by George Huang on 2019-11-21.
//

#include <cmath>
#include <iostream>
#include <chrono>
#include "testBLAS.cpp"
#include "BLAS.cpp"

std::chrono::time_point<std::chrono::system_clock> start, end;

int main() {
 /// Benchmark for dot product
 int N = 1000000;
 auto x0 = new double[N]();
 auto x1 = new double[N]();
 for(int i = 0; i < N; i++) {
  x0[i] = ((double) random() / (RAND_MAX));
  x1[i] = ((double) random() / (RAND_MAX));
 }
 int num_test = 100;

 // test basic dot_prod
 double result0;
 start = std::chrono::system_clock::now();
 for(int i = 0; i < num_test; i++)
  result0 = dot_prod(N, x0, x1);
 end = std::chrono::system_clock::now();
 auto elapsed_seconds = end - start;
 std::cout << elapsed_seconds.count() << std::endl;

 // test dotprod01
 double result1;
 start = std::chrono::system_clock::now();
 for(int i = 0; i < num_test; i++)
  result1 = dotprod01(N, x0, x1);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end - start;
 std::cout << elapsed_seconds.count() << std::endl;

 double result2;
 start = std::chrono::system_clock::now();
 for(int i = 0; i < num_test; i++)
  result2 = dotprod02(N, x0, x1);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end - start;
 std::cout << elapsed_seconds.count() << std::endl;

 std::cout << std::fabs(result0 - result1) << std::endl;
 std::cout << std::fabs(result1 - result2) << std::endl;
}