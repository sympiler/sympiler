//
// Created by kazem on 4/20/17.
//

#include <fstream>
#include "Matrix.h"
#include "Util.h"

namespace Sympiler {
namespace Internal {

MatrixPattern::MatrixPattern() { }


MatrixPattern::~MatrixPattern() {
    delete []p;
    delete []i;
}

bool MatrixPattern::read() {//TODO allocate p an i
    return false;
}

int Matrix::mNo=0;

Matrix::Matrix() { }

Matrix::Matrix(Type t, int order_num, int dim,
               std::string path, std::string name)
        :t(t),name(name),dim(dim),path(path) {
    incNo();
    order="n"+std::to_string(mNo);
}

Matrix::Matrix(const Matrix &m) {
    t=m.t;
    name=m.name;
    dim=m.dim;
    path=m.path;
    order = m.order;
}


Expr Matrix::Order() {
    return Variable::make(t,order);
}

void Matrix::setOrder_Nym(int oNo) {
    if(oNo>0)
        mPattern->n=oNo;
    else
        mPattern->n=0;
}

Expr Matrix::diagonal(Expr e){//TODO: should make it virtual=0 or?
    return e;
}

Expr Matrix::accessRowIdx(Expr e) {
    return e;
}

Expr Matrix::accessRowIdxPntr(Expr e){
    return e;
}

Expr Matrix::accessCol(Expr e) {
    return e;
}

Expr Matrix::accessNNZ(Expr e) {
    return e;
}

Stmt Matrix::allocateCol(){
    Stmt s;
    return s;
}

Stmt Matrix::allocateRow(){
    Stmt s;
    return s;
}

Stmt Matrix::allocateNNZ(){
    Stmt s;
    return s;
}

/*bool Matrix::readPattern() {

    return false;
}*/

MatrixPattern* Matrix::readPattern() {
    return mPattern;
}

Sparse::Sparse(int nnz, Type t, int order, int dim, std::string path, std::string name):
        nnz(nnz),Matrix(t,order,dim,path,name),isBlock(true) {
    std::string mnumber = std::to_string(mNo);
    row = Matrix::name+"i"+mnumber;
    col = Matrix::name+"p"+mnumber;
    nz = Matrix::name+"x"+mnumber;
 //   if(isBlock)
        rowP = Matrix::name+"ip"+mnumber;
}

Sparse::Sparse(Type t, std::string path):
        Matrix(t,0,2,path,"M"){
    std::string mnumber = std::to_string(mNo);
    row = Matrix::name+"i"+mnumber;
    col = Matrix::name+"p"+mnumber;
    nz = Matrix::name+"x"+mnumber;
  //  if(isBlock)
        rowP = Matrix::name+"ip"+mnumber;
    mPattern = new MatrixPattern();
}

Sparse::~Sparse() {
    delete mPattern;
}

Expr Sparse::diagonal(Expr e) {
    Expr idx = Pointer::make(halide_type_t(halide_type_int,32),col,e);
    return Pointer::make(halide_type_t(halide_type_float,64),nz, idx);
}

Expr Sparse::accessRowIdx(Expr e) {
    return Pointer::make(halide_type_t(halide_type_int,32),row,e);
}

Expr Sparse::accessRowIdxPntr(Expr e) {
    return Pointer::make(halide_type_t(halide_type_int,32),rowP,e);
}

Expr Sparse::accessNNZ(Expr e) {
    return Pointer::make(halide_type_t(halide_type_float,64),nz,e);
}

Expr Sparse::accessCol(Expr e) {
    return Pointer::make(halide_type_t(halide_type_int,32),col,e);
}

Stmt Sparse::allocateCol(){//TODO will be needed in future
    Stmt s;
    return s;
}

Stmt Sparse::allocateRow(){
    Stmt s;
    return s;
}

Stmt Sparse::allocateNNZ(){
    Stmt s;
    return s;
}

void Sparse::getDecl(std::vector<Expr>& exprList,
                     std::vector<Argument>& argList){
  //  exprList.push_back(Variable::make(Matrix::order));
    argList.push_back(Argument(order,Argument::Kind::InputScalar,
                               halide_type_t(halide_type_int,32),0));
   // exprList.push_back(Variable::make("*"+col));
    argList.push_back(Argument(col,Argument::Kind::InputBuffer,
                               halide_type_t(halide_type_int,32),1));
   // exprList.push_back(Variable::make("*"+row));
    argList.push_back(Argument(row,Argument::Kind::InputBuffer,
                               halide_type_t(halide_type_int,32),1));
   // exprList.push_back(Variable::make("*"+nnz));
    argList.push_back(Argument(nz,Argument::Kind::InputBuffer,
                               halide_type_t(halide_type_float,64),1));
   // if(isBlock)
        argList.push_back(Argument(rowP,Argument::Kind::InputBuffer,
                                   halide_type_t(halide_type_int,32),1));
}

/*bool Sparse::readPattern(){
    mPattern = new MatrixPattern();
    readCSCMatrixPattern(mPattern->path, order_num, mPattern->nnz,
                         mPattern->p,mPattern->i);
    return true;
}*/

MatrixPattern* Sparse::readPattern() {

    readCSCMatrixPattern(path, mPattern->n, mPattern->m, mPattern->nnz,
                         mPattern->p,mPattern->i);
    return mPattern;
}

Dense::Dense():Matrix() { }

Dense::Dense(Type t, int order, int dim, std::string path, std::string name)
        :Matrix(t, order, dim, path, name ) {
    Matrix::name += std::to_string(mNo);
}

Dense::Dense(Type t, std::string path)
        :Matrix(t,0,1,path,"D") {
    Matrix::name += std::to_string(mNo);
}

Dense::Dense(Sparse sp) // might use colNo to point
    :Matrix(sp.getType(),0,1,sp.Path(),"D"){
    Matrix::name += std::to_string(mNo);

}

Dense::Dense(const Matrix& m) // might use colNo to point
        :Matrix(m){
    Matrix::name = "D" + std::to_string(mNo);
}

Expr Dense::diagonal(Expr e){
    return Pointer::make(halide_type_t(halide_type_float,64),Matrix::name, e);
}

Expr Dense::accessNNZ(Expr e) {
    return Pointer::make(halide_type_t(halide_type_float,64),Matrix::name,e);
}

void Dense::getDecl(std::vector<Expr>& exprList, std::vector<Argument>& argList){
    argList.push_back(Argument(name,Argument::Kind::InputBuffer,
                               halide_type_t(halide_type_float,64),1));
}


}
}
