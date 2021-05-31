//
// Created by Shujian Qian on 2021-05-30.
//

#ifndef SYMPILER_FRONT_VAR_H
#define SYMPILER_FRONT_VAR_H


#include <string>
#include <utility>

#include "front_expr.h"
#include "callable_object.h"


#define SYMLANG_FROM_SOURCE(fname) decltype(fname ""_ToSymLangStr)


namespace SymLang {
  template <char... chars>
  using TString = std::integer_sequence<char, chars...>;
} // namespace SymLang


template <typename T, T... chars>
constexpr SymLang::TString<chars...> operator""_ToSymLangStr() { return { }; }


namespace SymLang {

  template<typename TSource = void>
  class Var : public FrontEnd::Internal::Expr {
  private:
    static constexpr bool m_is_intermetiate = false;
    std::string m_source;

  protected:
    Var() = default;

  public:
    inline std::string get_source() const { return m_source; }
  };

  template<>
  class Var<void> : public FrontEnd::Internal::Expr {
  private:
    static constexpr bool m_is_intermetiate = true;
    std::string m_source;

  protected:
    Var() = default;

  public:
    inline std::string get_source() const { return m_source; }
  };

  template<char... Elements>
  class Var<TString<Elements...> > {
  private:
    static constexpr char m_source_str[sizeof...(Elements)+1] = {Elements..., '\0'};
    std::string m_source = (m_source_str);

  protected:
    Var() = default;

  public:
    inline std::string get_source() const { return m_source; }
  };

//  auto v = Var<decltype("something"_tostr)>();

}


#endif //SYMPILER_FRONT_VAR_H
