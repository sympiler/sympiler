//
// Created by Shujian Qian on 2021-05-30.
//

#ifndef SYMPILER_FRONT_SCALAR_H
#define SYMPILER_FRONT_SCALAR_H


#include <string>

#include "front_var.h"


namespace SymLang {

  template<typename TSource = void>
  class Scalar : public Var<TSource> {
  private:
    double m_value;
  public:
    Scalar(double value);

    Scalar() = default;

    template <typename _>
    Scalar<> operator+(Scalar<_> rhs);

  };

  template<typename TSource>
  Scalar<TSource>::Scalar(double value) : m_value(value) {}

  template<typename TSource>
  template<typename _>
  Scalar<> Scalar<TSource>::operator+(Scalar<_> rhs) {
    // TODO:: Implement
    return Scalar<>(0);
  }

}


#endif //SYMPILER_FRONT_SCALAR_H
