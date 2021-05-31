//
// Created by Shujian Qian on 2021-05-28.
//

#ifndef SYMPILER_FRONT_FUNC_H
#define SYMPILER_FRONT_FUNC_H


#include <memory>

#include "callable_object.h"


#define SYMLANG_FUNC(func_args...) SymLang::FrontEnd::Internal::FuncBuilder([](func_args){
#define SYMLANG_ENDFUNC }).get_func();


namespace SymLang::FrontEnd::Internal {
  template<typename TLambda, typename TTuple>
  class FuncBuilder;
} // namespace SymLang::FrontEnd::Internal


namespace SymLang {

  /**
   * \brief A class that represents a SymLang function.
   */
  class Func {
    // TODO: Implement
  private:
    Func() = default;

  public:
    /**
     * \brief Output SymLang Func to C header file.
     */
    void compile_to_header(const std::string &fname);

    template<typename TLambda, typename TTuple>
    friend
    class FrontEnd::Internal::FuncBuilder;
  };

} // namespace SymLang


namespace SymLang::FrontEnd::Internal {

  /**
   * \brief Utility class that builds a SymLang function.
   *
   * \tparam TLambda
   * \tparam TTuple
   */
  template<typename TLambda, typename TTuple = typename CallableObject<TLambda>::ArgsTuple>
  class FuncBuilder {
  private:
    TTuple m_args_tuple;
    TLambda m_lambda;

  public:
    /**
     * \brief Constructor from a callable object.
     *
     * \param lambda A callable object (with an operator() defined).
     */
    explicit FuncBuilder(TLambda lambda)
        : m_args_tuple(TTuple()),
          m_lambda(lambda) {}

    /**
     * \brief Create a SymLang::Func from the callable object.
     *
     * \return SymLang::Func object.
     */
    Func get_func();

  };

  template<typename TLambda, typename TTuple>
  Func FuncBuilder<TLambda, TTuple>::get_func() {
    // TODO: Implement
    return Func();
  }

} // namespace SymLang::Internal


#endif //SYMPILER_FRONT_FUNC_H
