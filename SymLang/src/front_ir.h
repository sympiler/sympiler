//
// Created by Shujian Qian on 2021-05-24.
//

/**
 * \file ir.h
 * \brief Implements high-level intermediate representation of SymLang.
 *
 *
 */

#ifndef SYMPILER_FRONT_IR_H
#define SYMPILER_FRONT_IR_H

#include <memory>

namespace SymLang::FrontEnd {

/**
 * \brief Visitor base class to traverse IR nodes.
 */
class IRVisitor {
  // TODO: Define virtual visit functions.
};

/**
 * \brief Base IR node class.
 */
class IRNode {
public:
  /**
   * \brief Virtual accept method that accepts visitors.
   *
   * \param visitor Visitor to traverse IR nodes.
   */
  virtual void accept(IRVisitor *visitor) = 0;

  virtual ~IRNode() = default;
};

} // namespace SymLang::FrontEnd

#endif //SYMPILER_FRONT_IR_H
