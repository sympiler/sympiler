//
// Created by kazem on 10/12/19.
//

#define DBG_LOG
#define CSV_LOG

#include <iostream>
#include <sparse_io.h>
#include <test_utils.h>
#include <omp.h>
#include <metis_interface.h>

#include "sptrsv_demo_utils.h"

using namespace sym_lib;

/// Evaluate spmv-sptrsv based on random matrices
/// \return
int sptrsv_csc_demo02(int argc, char *argv[]);

int main(int argc, char *argv[]){
 int ret_val;
 ret_val = sptrsv_csc_demo02(argc,argv);
 return ret_val;
}



int sptrsv_csc_demo02(int argc, char *argv[]){
 CSC *L1_csc, *A = NULLPNTR;
 CSR *L2_csr;
 size_t n;
 int num_threads = 6;
 int p2 = -1, p3 = 4000; // LBC params
 int header = 0;
 int *perm;
 std::string matrix_name;
 std::vector<timing_measurement> time_array;
 if (argc < 2) {
  PRINT_LOG("Not enough input args, switching to random mode.\n");
  n = 16;
  double density = 0.2;
  matrix_name = "Random_" + std::to_string(n);
  A = random_square_sparse(n, density);
  if (A == NULLPNTR)
   return -1;
  L1_csc = make_half(A->n, A->p, A->i, A->x);
 } else {
  std::string f1 = argv[1];
  matrix_name = f1;
  L1_csc = read_mtx(f1);
  if (L1_csc == NULLPNTR)
   return -1;
  n = L1_csc->n;
 }
 if(argc >= 3)
  p2 = atoi(argv[2]);
 omp_set_num_threads(num_threads);
 if(argc >= 4)
  p3 = atoi(argv[3]);
 /// Re-ordering L matrix
#ifdef METIS
 //We only reorder L since dependency matters more in l-solve.
 //perm = new int[n]();
 CSC *L1_csc_full = make_full(L1_csc);
 delete L1_csc;
 metis_perm_general(L1_csc_full, perm);
 L1_csc = make_half(L1_csc_full->n, L1_csc_full->p, L1_csc_full->i,
                    L1_csc_full->x);
 CSC *Lt = transpose_symmetric(L1_csc, perm);
 CSC *L1_ord = transpose_symmetric(Lt, NULLPNTR);
 delete L1_csc;
 L1_csc = L1_ord;
 delete Lt;
 delete L1_csc_full;
 delete[]perm;
#endif

 L2_csr = csc_to_csr(L1_csc);

 double *y_serial, *y_correct = new double[n];

 timing_measurement t_ser, t_par, t_par2, t_blocked, t_blocked_mkl,
 t_blocked_levelset, t_levelset;

 SptrsvSerial *ss = new SptrsvSerial(L2_csr, L1_csc, NULLPNTR, "serial");
 t_ser = ss->evaluate();
 y_serial = ss->solution();
 copy_vector(0,n,y_serial,y_correct);
 //print_vec("x:\n", 0, n, y_correct);

 auto *sls = new SptrsvLevelSet(L2_csr, L1_csc, y_correct, "levelset csc");
 t_levelset = sls->evaluate();

 auto *sl = new SptrsvLBC(L2_csr, L1_csc, y_serial, "lbc",num_threads, p2, p3);
 t_par = sl->evaluate();


 if(header)
  std::cout<<"Matrix Name,Metis Enabled,"
             "Number of Threads,"
             "Serial Non-fused,Parallel Levelset CSC,Parallel LBC CSR,";

 #ifdef METIS
 PRINT_CSV("METIS");
#else
 PRINT_CSV("No");
#endif
 PRINT_CSV(num_threads);
 PRINT_CSV(p2);
 PRINT_CSV(p3);
 PRINT_CSV(t_ser.elapsed_time);
 PRINT_CSV(t_levelset.elapsed_time);
 PRINT_CSV(t_par.elapsed_time);

 delete []y_correct;
 delete A;
 delete L1_csc;
 delete L2_csr;

 delete ss;
 delete sl;
 delete sls;


 return 0;
}
