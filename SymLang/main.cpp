//
// Created by Shujian Qian on 2021-05-30.
//

#include <iostream>
#include "SymLang.h"

int main() {

  SymLang::Func f =
      SYMLANG_FUNC(
          SymLang::Sparse<SYMLANG_FROM_SOURCE("somefile")> a,
          SymLang::Sparse<SYMLANG_FROM_SOURCE("somefile2")> b,
          SymLang::Scalar<SYMLANG_FROM_SOURCE("somefile3")> c,
          SymLang::Scalar<SYMLANG_FROM_SOURCE("somefile4")> d
      )
        auto e = c + d;
        return a;
      SYMLANG_ENDFUNC

  f.compile_to_header("sym_out.h");

  auto v = SymLang::Dense<SYMLANG_FROM_SOURCE("dense_input_file.mtx")>();
  std::cout << v.get_source() << std::endl;

  auto s = SymLang::Sparse<SYMLANG_FROM_SOURCE("sparse_input_file.mtx")>();
  std::cout << v.get_source() << std::endl;
}
