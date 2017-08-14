//
// Created by kazem on 6/13/17.
//

#ifndef SYMPILER_PROJ_FACTORIZATIONINSPECTOR_H
#define SYMPILER_PROJ_FACTORIZATIONINSPECTOR_H

#include "Inspector.h"

namespace Sympiler{
namespace Internal{

class FactorizationInspector: public Inspector{
protected:
    Matrix *A;
    int *eTree;
    int *postETree;
    int *colCount;

public:
    FactorizationInspector(Matrix *A);
    ~FactorizationInspector();
    virtual void inspectionGraph();
    virtual void pruneCalculation(){};
    virtual void blockDetection(){};
    virtual SymbolicObject* strategy();

};

class CholeskyInspector: public FactorizationInspector{
protected:

public:
    CholeskyInspector(Matrix *A);
    ~CholeskyInspector();
    virtual void inspectionGraph();
    virtual void pruneCalculation();
    virtual void blockDetection();
    virtual SymbolicObject* strategy(Tuning );

};

}
}

#endif //SYMPILER_PROJ_FACTORIZATIONINSPECTOR_H
