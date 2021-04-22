//
// Created by george on 2020-02-22.
//

#include "IRProfiler.h"
#include "IRPrinter.h"

#include <utility>
#include <iostream>
#include <sstream>


namespace Sympiler {
namespace Internal {


void IRProfiler::use_var(int line, const std::string& name) {
    auto cache = variables.find(name);
    if(cache != variables.end()) {
        cache->second->uses.emplace_back(line);
    }
}



IRProfiler::IRProfiler() {
    cur_lines = 0;
    id = id_ptr = "";
    new_v = "";
    variables.clear();
    expr_map.clear();

    printer = new IRPrinter(ss);

    loops.clear();
    parseloop = false;
    loopscope = -1;
}

// TODO: print the AST using DOT file
void IRProfiler::print_AST(const std::string &fname, Kernel &ker) {
    print_ast = true;
    op_count = 0; // for unique operator output
    stream.open(fname.c_str(), std::ios::out | std::ios::trunc);

    stream << "graph AST {\n";
    visit_stmt(ker.LoweredKer());
    stream << "}\n";
}



/** Comparison functions **/

std::string IRProfiler::visit_expr(Expr e) {
    id = "$$ BAD ID $$";
    e.accept(this);
    return id;
}

void IRProfiler::visit_stmt(Stmt s) {
    s.accept(this);
}

void IRProfiler::start_analysis(Module module) {
    start_analysis(module.kers[0]);
}

void IRProfiler::start_analysis(Kernel &ker) {
    // store the declarations of argument
    for (auto & i : ker.argType)
        new_var(i.name);
    visit_stmt(ker.LoweredKer());
}


void IRProfiler::print_trace() {
//    auto it = variables.begin();
//    while (it != variables.end()) {
//        std::cout << it->first << " "
//            << it->second->def << " "
//            << it->second->dead << " "
//            << it->second->array << "\n\t";
//
//        for(auto &i : it->second->uses)
//            std::cout << i << " ";
//        std::cout << "\n---------------------------\n";
//        it++;
//    }
    std::cout << loopstack.size() << "\n";

    std::cout << "-------------- LOOPS ----------------\n";
    for(auto &l : loops) {
        std::cout << l->start << " " << l->end << " " << l->scope << "\n";
        std::cout << "\t\t";
        auto it = l->reads.begin();
        while (it != l->reads.end()) {
            std::cout << it->c_str() << " ";
            it++;
        }
        std::cout << "\n\t\t";
        it = l->writes.begin();
        while(it != l->writes.end()) {
            std::cout << it->c_str() << " ";
            it++;
        }
        std::cout << "\n===================\n";
    }
}


void IRProfiler::visit_binop(Type t, Expr a, Expr b, const char *op) {
    std::string a_str, b_str;
    loop *parentloop;

    auto cur_var1 = new_v;
    a_str = visit_expr(std::move(a));
    if(new_v != cur_var1) {
        use_var(cur_lines, new_v);
        cur_lines++;
    }

    auto cur_var2 = new_v;
    b_str = visit_expr(std::move(b));
    if(cur_var2 != new_v) {
        use_var(cur_lines, new_v);
        cur_lines++;
    }

    std::string name = "("+a_str+" "+op+" "+b_str+")";
    std::string s = store_name(name);
    if(s[0] == '_') {
        new_var(s);
        cur_lines++;
    }
    id = name;
}

void IRProfiler::visit(const Variable *v) {
    id = v->name;
    use_var(cur_lines, v->name);
}

void IRProfiler::visit(const Pointer *p) {
    id = p->name;
    scope_var.emplace_back(id);
    id_ptr = p->name + "[" + visit_expr(p->idx) + "]";
    id = id_ptr;
    id_ptr = "";
    use_var(cur_lines, p->name);
}

void IRProfiler::visit(const IntImm *i) {
    id = std::to_string(i->value);
}

void IRProfiler::visit(const UIntImm *i) {
    id = std::to_string(i->value);
}

void IRProfiler::visit(const FloatImm *f) {
    id = std::to_string(f->value);
}

void IRProfiler::visit(const StringImm *s) {
    id = s->value;
}

void IRProfiler::visit(const Add *op) {
    visit_binop(op->type, op->a, op->b, "+");
}

void IRProfiler::visit(const Sub *op) {
    visit_binop(op->type, op->a, op->b, "-");
}

void IRProfiler::visit(const Mul *op) {
    visit_binop(op->type, op->a, op->b, "*");
}

void IRProfiler::visit(const Div *op) {
    visit_binop(op->type, op->a, op->b, "/");
}

void IRProfiler::visit(const Mod *op) {
    visit_binop(op->type, op->a, op->b, "%");
}

void IRProfiler::visit(const Min *op) {
    visit_expr(Call::make(op->type, "min", {op->a, op->b}));
}

void IRProfiler::visit(const Max *op) {
    visit_expr(Call::make(op->type, "max", {op->a, op->b}));
}

void IRProfiler::visit(const EQ *op) {
    visit_binop(op->type, op->a, op->b, "==");
}

void IRProfiler::visit(const NE *op) {
    visit_binop(op->type, op->a, op->b, "!=");
}

void IRProfiler::visit(const LT *op) {
    visit_binop(op->type, op->a, op->b, "<");
}

void IRProfiler::visit(const LE *op) {
    visit_binop(op->type, op->a, op->b, "<=");
}

void IRProfiler::visit(const GT *op) {
    visit_binop(op->type, op->a, op->b, ">");
}

void IRProfiler::visit(const GE *op) {
    visit_binop(op->type, op->a, op->b, ">=");
}

void IRProfiler::visit(const And *op) {
    visit_binop(op->type, op->a, op->b, "&&");
}

void IRProfiler::visit(const Or *op) {
    visit_binop(op->type, op->a, op->b, "||");
}

void IRProfiler::visit(const Not *op) {

}

void IRProfiler::visit(const Select *op) {
}

void IRProfiler::visit(const For *op) {
    // set up for loop
    loopscope++;
    parseloop = true;
    loop *l = new loop();
    l->l = op;
    l->start = cur_lines;
    l->scope = loopscope;
    l->innerloops.clear();

    loopstack.push(l);

    cur_lines++;
    auto cur_var1 = new_v;
    visit_expr(op->min);
    if(new_v != cur_var1)
        use_var(cur_lines, new_v);
    new_var(op->name);
    cur_lines++;

    auto cur_var2 = new_v;
    visit_expr(op->extent);
    if(new_v != cur_var2) {
        use_var(cur_lines, new_v);
        cur_lines++;
    }
    visit_stmt(op->body);

    // loop finished
    l->end = cur_lines;
    loopstack.pop();
    loopscope--;
    if(loopscope < 0)
        parseloop = false;
    loops.emplace_back(l);
}

void IRProfiler::visit(const Let *op) {

}

void IRProfiler::visit(const LetStmt *op) {
}

void IRProfiler::visit(const LetStmtPtr *op) {

}

void IRProfiler::visit(const Evaluate *op) {
//    scope_var.clear();
    auto name = visit_expr(op->value);

    cur_lines++;
    auto s = store_name(name);
    new_var(s);

//    if(parseloop) {
//        loop *parentloop = loopstack.top();
//        parentloop->reads.insert(scope_var.begin(), scope_var.end());
//    }
}

void IRProfiler::visit(const Block *op) {
    visit_stmt(op->first);
    visit_stmt(op->rest);
}

void IRProfiler::visit(const Load *op) {

}

// FIXME: does Store have to be in for loops?
void IRProfiler::visit(const Store *op) {
    loop *parentloop = loopstack.top();

    // process buffer to get name
    scope_var.clear();
    auto cur_var1 = new_v;
    visit_expr(op->value);
    if(new_v != cur_var1) {
        use_var(cur_lines, new_v);
    }
    parentloop->reads.insert(scope_var.begin(), scope_var.end());

    scope_var.clear();
    auto cur_var2 = new_v;
    visit_expr(op->ptr);
    if(new_v != cur_var2) {
        use_var(cur_lines, new_v);
        cur_lines++;
    }
    int nvar = scope_var.size();
    parentloop->writes.insert(scope_var[0]);
    parentloop->reads.insert(scope_var.begin()+1, scope_var.end());
}

void IRProfiler::visit(const Call *op) {

}

void IRProfiler::visit(const CallX *op) {

}

void IRProfiler::visit(const Allocate *op) {

}

void IRProfiler::visit(const Free *op) {

}

void IRProfiler::visit(const IfThenElse *op) {

}

void IRProfiler::clear_map() {
    auto it = variables.begin();
    while (it != variables.end()) {
        it->second->uses.clear();
        delete(it->second);
        it++;
    }
    variables.clear();
    expr_map.clear();
}

std::string IRProfiler::store_name(const std::string& s) {
    auto cached = expr_map.find(s);
    if (cached == expr_map.end()) {
        std::string name = unique_name('_');
        expr_map[s] = name;
        return name;
    } else {
        return cached->second;
    }
}

var *IRProfiler::new_var(const std::string &s) {
    var *v = new var();
    v->def = cur_lines;
//    v->array = array;
    v->dead = INT32_MAX;
    v->uses.clear();
    variables[s] = v;

    new_v = s;
    return v;
}

}
}