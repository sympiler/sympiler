//
// Created by Shujian Qian on 2021-05-30.
//

#ifndef SYMPILER_CALLABLE_OBJECT_H
#define SYMPILER_CALLABLE_OBJECT_H

#include <tuple>
#include <type_traits>

namespace SymLang::FrontEnd::Internal {

/**
 * \brief Template metaprogram that checks if a type has a call operator.
 *
 * \tparam T
 */
template <typename T>
class HasCallOperator
{
  typedef char Correct;
  struct Incorrect { char _[2]; };

  template <typename C> static Correct test(decltype(&C::operator())) { }
  template <typename C> static Incorrect test(C) { }
public:
  /**
   * \brief Does the type have a call operator?
   */
  static constexpr bool value = sizeof(test<T>(0)) == sizeof(Correct);
};

/**
 * \brief Base of a template metaprogram.
 *
 * \tparam TFuncType
 */
template<typename TFuncType>
struct MemberFunction {
};

/**
 * \brief Template metaprogram that extracts the arity, args types and return
 *        type of a callable type.
 *
 * \tparam TClass
 * \tparam TReturn
 * \tparam TArgs
 */
template<typename TClass, typename TReturn, typename... TArgs>
struct MemberFunction<TReturn (TClass::*)(TArgs...)> {
  static constexpr size_t arity = sizeof...(TArgs); /**< arity of the callable */

  typedef TReturn return_type; /**< return type of the callable */
  typedef TClass declare_type; /**< declared type of the callable */

  /**
   * \brief Get type of the idx-th argument.
   */
  template<size_t idx>
  using get_arg_type = typename std::tuple_element<idx, std::tuple<TArgs...> >::type;

  typedef std::tuple<TArgs...> ArgsTuple; /**< tuple of args of callable */
  typedef std::tuple<typename std::decay<TArgs>::type...> DecayedArgsTuple; /**< tuple of decayed types of args of callable */
};

/**
 * \brief Template metaprogram that extracts the arity, args types and return
 *        type of a callable type.
 *
 * \tparam TClass
 * \tparam TReturn
 */
template<typename TClass, typename TReturn>
struct MemberFunction<TReturn (TClass::*)()>
{
  static constexpr size_t arity = 0; /**< arity of the callable */

  typedef TReturn return_type; /**< return type of the callable */
  typedef TClass declare_type; /**< declared type of the callable */

  /**
   * \brief Get type of the idx-th argument.
   */
  template <size_t idx>
  using Args = void;

  typedef std::tuple<> ArgsTuple; /**< tuple of args of callable */
  typedef std::tuple<> DecayedArgsTuple; /**< tuple of decayed types of args of callable */

};

/**
 * \brief Template metaprogram that extracts the arity, args types and return
 *        type of a callable type.
 *
 * \tparam TClass
 * \tparam TReturn
 * \tparam TArgs
 */
template <typename TClass, typename TReturn, typename... TArgs>
struct MemberFunction<TReturn (TClass::*)(TArgs...) const> :
    MemberFunction<TReturn (TClass::*)(TArgs...)> { };

/**
 * \brief Template metaprogram to detect if a type is callable. Get the arg
 *        types of the call operator if callable.
 *
 * \tparam T
 */
template<typename T, bool = HasCallOperator<T>::value>
struct CallableObject;

template <typename T>
struct CallableObject<T, false> {
  static constexpr bool is_callable = false;
};

template <typename T>
struct CallableObject<T, true> : MemberFunction<decltype(&T::operator())> {
  static constexpr bool is_callable = true;
};

} // namespace SymLang::Internal::FrontEnd

#endif //SYMPILER_CALLABLE_OBJECT_H
