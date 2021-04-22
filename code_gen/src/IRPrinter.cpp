//
// Created by kazem on 13/02/16.
//

#include "IRPrinter.h"
#include "IROperator.h"

namespace Sympiler {
namespace Internal {

    using std::ostream;
    using std::vector;
    using std::string;
    using std::ostringstream;

    std::ostream &operator<<(std::ostream &out, const Type &type) {
        switch (type.code()) {
            case Type::Int:
                out << "int";
                break;
            case Type::UInt:
                out << "uint";
                break;
            case Type::Float:
                out << "float";
                break;
            case Type::Handle:
                out << "handle";
                break;
        }
        out << type.bits();
        if (type.lanes() > 1) out << 'x' << type.lanes();
        return out;
    }

    ostream &operator<<(ostream &stream, const Expr &ir) {
        if (!ir.defined()) {
            stream << "(undefined)";
        } else {
            IRPrinter p(stream);
            p.print(ir);
        }
        return stream;
    }

/*
    ostream &operator <<(ostream &stream, const Buffer &buffer) {
        return stream << "buffer " << buffer.name() << " = {...}\n";
    }
*/

/*    ostream &operator<<(ostream &stream, const Module &m) {
        stream << "Target = " << m.target().to_string() << "\n";
        for (size_t i = 0; i < m.buffers.size(); i++) {
            stream << m.buffers[i] << "\n";
        }
        for (size_t i = 0; i < m.functions.size(); i++) {
            stream << m.functions[i] << "\n";
        }
        return stream;
    }*/

    ostream &operator<<(ostream &out, const DeviceAPI &api) {
        switch (api) {
            case DeviceAPI::Host:
                break;
            case DeviceAPI::Parent:
                out << "<Parent>";
                break;
            case DeviceAPI::Default_GPU:
                out << "<Default_GPU>";
                break;
            case DeviceAPI::CUDA:
                out << "<CUDA>";
                break;
            case DeviceAPI::OpenCL:
                out << "<OpenCL>";
                break;
            case DeviceAPI::OpenGLCompute:
                out << "<OpenGLCompute>";
                break;
            case DeviceAPI::GLSL:
                out << "<GLSL>";
                break;
            case DeviceAPI::Renderscript:
                out << "<Renderscript>";
                break;
            case DeviceAPI::Metal:
                out << "<Metal>";
                break;
        }
        return out;
    }

IRPrinter::IRPrinter(std::ostream &s) : stream(s), indent(0) {
    s.setf(std::ios::fixed, std::ios::floatfield);
}

void IRPrinter::print(Expr ir) {
    ir.accept(this);
}

void IRPrinter::print(Stmt ir) {
    ir.accept(this);
}


void IRPrinter::do_indent() {
    for (int i = 0; i < indent; i++) stream << ' ';
}

void IRPrinter::visit(const IntImm *op) {
    if (op->type == Int(32)) {
        stream << op->value;
    } else {
        stream << "(" << op->type << ")" << op->value;
    }
}

void IRPrinter::visit(const UIntImm *op) {
    stream << "(" << op->type << ")" << op->value;
}

void IRPrinter::visit(const FloatImm *op) {
    switch (op->type.bits()) {
        case 64:
            stream << op->value;
            break;
        case 32:
            stream << op->value << 'f';
            break;
        case 16:
            stream << op->value << 'h';
            break;
        /*default:
            internal_error << "Bad bit-width for float: " << op->type << "\n";*/
    }
}

void IRPrinter::visit(const Variable *v) {
    stream << v->name;
}
void IRPrinter::visit(const Pointer *p) {
    stream << p->name;// << "[";
//    print(p->idx);
//    stream << "]";
}


void IRPrinter::visit(const StringImm *op) {
    stream << '"';
    for (size_t i = 0; i < op->value.size(); i++) {
        unsigned char c = op->value[i];
        if (c >= ' ' && c <= '~' && c != '\\' && c != '"') {
            stream << c;
        } else {
            stream << '\\';
            switch (c) {
                case '"':
                    stream << '"';
                    break;
                case '\\':
                    stream << '\\';
                    break;
                case '\t':
                    stream << 't';
                    break;
                case '\r':
                    stream << 'r';
                    break;
                case '\n':
                    stream << 'n';
                    break;
                default:
                    string hex_digits = "0123456789ABCDEF";
                    stream << 'x' << hex_digits[c >> 4] << hex_digits[c & 0xf];
            }
        }
    }
    stream << '"';
}

void IRPrinter::visit(const Add *op) {
    stream << '(';
    print(op->a);
    stream << " + ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const Sub *op) {
    stream << '(';
    print(op->a);
    stream << " - ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const Mul *op) {
    stream << '(';
    print(op->a);
    stream << " * ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const Div *op) {
    stream << '(';
    print(op->a);
    stream << " / ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const Mod *op) {
    stream << '(';
    print(op->a);
    stream << " % ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const Min *op) {
    stream << "min(";
    print(op->a);
    stream << ", ";
    print(op->b);
    stream << ")";
}

void IRPrinter::visit(const Max *op) {
    stream << "max(";
    print(op->a);
    stream << ", ";
    print(op->b);
    stream << ")";
}

void IRPrinter::visit(const EQ *op) {
    stream << '(';
    print(op->a);
    stream << " == ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const NE *op) {
    stream << '(';
    print(op->a);
    stream << " != ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const LT *op) {
    stream << '(';
    print(op->a);
    stream << " < ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const LE *op) {
    stream << '(';
    print(op->a);
    stream << " <= ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const GT *op) {
    stream << '(';
    print(op->a);
    stream << " > ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const GE *op) {
    stream << '(';
    print(op->a);
    stream << " >= ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const And *op) {
    stream << '(';
    print(op->a);
    stream << " && ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const Or *op) {
    stream << '(';
    print(op->a);
    stream << " || ";
    print(op->b);
    stream << ')';
}

void IRPrinter::visit(const Not *op) {
    stream << '!';
    print(op->a);
}

void IRPrinter::visit(const Select *op) {
    stream << "select(";
    print(op->condition);
    stream << ", ";
    print(op->true_value);
    stream << ", ";
    print(op->false_value);
    stream << ")";
}

void IRPrinter::visit(const Evaluate *op){
    print(op->value);
}

void IRPrinter::visit(const Block *op) {
    print(op->first);
    if (op->rest.defined()) print(op->rest);
}

void IRPrinter::visit(const Let *op) {
    stream << "(let " << op->name << " = ";
    print(op->value);
    stream << " in ";
    print(op->body);
    stream << ")";
}

void IRPrinter::visit(const LetStmt *op) {
    do_indent();
    stream << "let " << op->name << " = ";
    print(op->value);
    stream << '\n';

    print(op->body);
}

void IRPrinter::visit(const Load *op){
    //print(op->ptr);
}

void IRPrinter::visit(const Store *op){
    //print(op->ptr);
}

void IRPrinter::visit(const Call *op) {
    // Special-case some intrinsics for readability
/*    if (op->call_type == Call::Intrinsic) {
        if (op->name == Call::extract_buffer_host) {
            print(op->args[0]);
            stream << ".host";
            return;
        } else if (op->name == Call::extract_buffer_min) {
            print(op->args[0]);
            stream << ".min[" << op->args[1] << "]";
            return;
        } else if (op->name == Call::extract_buffer_max) {
            print(op->args[0]);
            stream << ".max[" << op->args[1] << "]";
            return;
        }
    }*/

    stream << op->name << "(";
    for (size_t i = 0; i < op->args.size(); i++) {
        print(op->args[i]);
        if (i < op->args.size() - 1) {
            stream << ", ";
        }
    }
    stream << ")";
}

void IRPrinter::visit(const CallX *op) {
    // Special-case some intrinsics for readability
/*    if (op->call_type == Call::Intrinsic) {
if (op->name == Call::extract_buffer_host) {
    print(op->args[0]);
    stream << ".host";
    return;
} else if (op->name == Call::extract_buffer_min) {
    print(op->args[0]);
    stream << ".min[" << op->args[1] << "]";
    return;
} else if (op->name == Call::extract_buffer_max) {
    print(op->args[0]);
    stream << ".max[" << op->args[1] << "]";
    return;
}
}*/

    stream << op->name << "(";
    for (size_t i = 0; i < op->args.size(); i++) {
        print(op->args[i]);
        if (i < op->args.size() - 1) {
            stream << ", ";
        }
    }
    stream << ")";
}

void IRPrinter::visit(const Allocate *op) {
    do_indent();
    stream << "allocate " << op->name << "[" << op->type;
    for (size_t i = 0; i < op->extents.size(); i++) {
        stream  << " * ";
        print(op->extents[i]);
    }
    stream << "]";
    if (!is_one(op->condition)) {
        stream << " if ";
        print(op->condition);
    }
    if (op->new_expr.defined()) {
        stream << "\n custom_new { " << op->new_expr << " }";
    }
    if (!op->free_function.empty()) {
        stream << "\n custom_delete { " << op->free_function << "(<args>); }";
    }
    stream << "\n";
    print(op->body);
}

void IRPrinter::visit(const Free *op) {
    do_indent();
    stream << "free " << op->name;
    stream << '\n';
}

void IRPrinter::visit(const IfThenElse *op) {
    do_indent();
    while (1) {
        stream << "if (" << op->condition << ") {\n";
        indent += 2;
        print(op->then_case);
        indent -= 2;

        if (!op->else_case.defined()) {
            break;
        }

        if (const IfThenElse *nested_if = op->else_case.as<IfThenElse>()) {
            do_indent();
            stream << "} else ";
            op = nested_if;
        } else {
            do_indent();
            stream << "} else {\n";
            indent += 2;
            print(op->else_case);
            indent -= 2;
            break;
        }
    }

    do_indent();
    stream << "}\n";

}

}
}


