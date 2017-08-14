//
// Created by kazem on 8/8/17.
//

#ifndef SYMPILER_PROJ_TRIANGULAR_H
#define SYMPILER_PROJ_TRIANGULAR_H

#include "Kernel.h"

namespace Sympiler {
    namespace Internal {
        class Triangular : public Kernel {
            Matrix *L, *rhs, *rhsDense;
            Expr rhsCol;

            Stmt fwdSolve();

            Stmt copy2Dense();

        public:
            Triangular();

            Triangular(Matrix &L, Matrix &rhs);

            ~Triangular();

            virtual Stmt baseCode();

            virtual Stmt
            VSBlockIG(SymbolicObject *sym); //Inspector-guided transformations
            virtual Stmt VIPruneIG(SymbolicObject *sym);

            virtual void sympile_to_c(std::string fName, Target t);

            virtual void sympile_to_c(std::string fName);
        };
    }
}


#endif //SYMPILER_PROJ_TRIANGULAR_H
