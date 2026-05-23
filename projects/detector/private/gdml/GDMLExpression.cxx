#include "GDMLParserPrivate.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

namespace siren {
namespace detector {
namespace gdml {
// Safe attribute value accessor: returns empty string if attribute is missing
const char* SafeAttrVal(rapidxml::xml_node<>* node, const char* attr_name) {
    if(!node) return "";
    rapidxml::xml_attribute<>* attr = node->first_attribute(attr_name);
    if(!attr) return "";
    return attr->value();
}

// Trim leading and trailing whitespace from a string
std::string TrimWhitespace(std::string const & s) {
    size_t start = s.find_first_not_of(" \t\n\r");
    if(start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\n\r");
    return s.substr(start, end - start + 1);
}


// Iterative expression evaluator using the shunting-yard algorithm.
// Tokenizes the input, converts infix to postfix (RPN), and evaluates
// the RPN with an explicit value stack. No recursion, O(1) call-stack depth.

enum class TokenKind { NUMBER, OP, LPAREN, RPAREN, COMMA, FUNC, UNARY_NEG };

struct Token {
    TokenKind kind;
    double value;
    char op;
    std::string name;
    int argc;
};

static int OpPrecedence(char op) {
    if(op == '+' || op == '-') return 1;
    if(op == '*' || op == '/') return 2;
    if(op == '^') return 3;
    return 0;
}

static bool IsRightAssociative(char op) {
    return op == '^';
}

// IsIdentStart and IsIdentChar are in GDMLParserPrivate.h

static double ApplyFunc(std::string const & name, double const * args, int argc, std::string const & expr) {
    if(argc == 1) {
        double a = args[0];
        if(name == "sin") return std::sin(a);
        if(name == "cos") return std::cos(a);
        if(name == "tan") return std::tan(a);
        if(name == "asin") return std::asin(a);
        if(name == "acos") return std::acos(a);
        if(name == "atan") return std::atan(a);
        if(name == "exp") return std::exp(a);
        if(name == "log") return std::log(a);
        if(name == "sqrt") return std::sqrt(a);
        if(name == "abs") return std::fabs(a);
    } else if(argc == 2) {
        double a = args[0], b = args[1];
        if(name == "pow") return std::pow(a, b);
        if(name == "atan2") return std::atan2(a, b);
        if(name == "min") return std::min(a, b);
        if(name == "max") return std::max(a, b);
    }
    throw std::runtime_error("GDML expression error: unknown function '" + name + "' in '" + expr + "'");
}

double EvalExpression(std::string const & expr, std::unordered_map<std::string, double> const & constants,
                      std::unordered_map<std::string, std::pair<int, std::vector<double>>> const & matrices,
                      int depth);

static std::vector<Token> Tokenize(std::string const & s, std::unordered_map<std::string, double> const & constants, std::string const & expr,
                                   std::unordered_map<std::string, std::pair<int, std::vector<double>>> const & matrices = {},
                                   int depth = 0) {
    std::vector<Token> tokens;
    size_t i = 0;
    while(i < s.size()) {
        char c = s[i];
        if(c == ' ' || c == '\t' || c == '\n' || c == '\r') { ++i; continue; }

        if(c == '(') { tokens.push_back({TokenKind::LPAREN}); ++i; continue; }
        if(c == ')') { tokens.push_back({TokenKind::RPAREN}); ++i; continue; }
        if(c == ',') { tokens.push_back({TokenKind::COMMA}); ++i; continue; }

        if(c == '+' || c == '-' || c == '*' || c == '/' || c == '^') {
            bool is_unary = (c == '-' || c == '+') &&
                (tokens.empty() || tokens.back().kind == TokenKind::LPAREN ||
                 tokens.back().kind == TokenKind::OP || tokens.back().kind == TokenKind::COMMA ||
                 tokens.back().kind == TokenKind::UNARY_NEG);
            if(is_unary && c == '-') {
                tokens.push_back({TokenKind::UNARY_NEG});
            } else if(is_unary && c == '+') {
                // unary plus: skip
            } else {
                Token t; t.kind = TokenKind::OP; t.op = c;
                tokens.push_back(t);
            }
            ++i; continue;
        }

        if(std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
            char* end = nullptr;
            double val = std::strtod(s.c_str() + i, &end);
            Token t; t.kind = TokenKind::NUMBER; t.value = val;
            tokens.push_back(t);
            i += (end - (s.c_str() + i)); continue;
        }

        if(IsIdentStart(c)) {
            size_t start = i;
            while(i < s.size() && IsIdentChar(s[i])) ++i;
            std::string name = s.substr(start, i - start);
            size_t j = i;
            while(j < s.size() && (s[j] == ' ' || s[j] == '\t')) ++j;
            if(j < s.size() && s[j] == '(') {
                Token t; t.kind = TokenKind::FUNC; t.name = name; t.argc = 1;
                tokens.push_back(t);
            } else if(j < s.size() && s[j] == '[') {
                // Matrix/vector element access: name[idx] or name[row,col]
                auto mit = matrices.find(name);
                if(mit == matrices.end()) {
                    throw std::runtime_error("GDML expression error: unknown matrix '" + name + "' in '" + expr + "'");
                }
                i = j + 1; // skip '['
                // Parse first index expression (up to ',' or ']')
                int bracket_depth = 1;
                size_t idx_start = i;
                size_t comma_pos = std::string::npos;
                while(i < s.size() && bracket_depth > 0) {
                    if(s[i] == '[') ++bracket_depth;
                    else if(s[i] == ']') --bracket_depth;
                    else if(s[i] == ',' && bracket_depth == 1) comma_pos = i;
                    if(bracket_depth > 0) ++i;
                }
                if(bracket_depth != 0) {
                    throw std::runtime_error("GDML expression error: unmatched '[' in '" + expr + "'");
                }
                size_t bracket_end = i;
                ++i; // skip ']'

                int coldim = mit->second.first;
                auto const & vals = mit->second.second;
                double result;

                if(comma_pos != std::string::npos) {
                    // Two indices: name[row, col] (1-based)
                    std::string row_expr = s.substr(idx_start, comma_pos - idx_start);
                    std::string col_expr = s.substr(comma_pos + 1, bracket_end - comma_pos - 1);
                    // Recursively evaluate the index expressions
                    double row_val = EvalExpression(row_expr, constants, {}, depth + 1);
                    double col_val = EvalExpression(col_expr, constants, {}, depth + 1);
                    int row = static_cast<int>(row_val + 0.5);
                    int col = static_cast<int>(col_val + 0.5);
                    int flat_idx = (row - 1) * coldim + (col - 1);
                    if(flat_idx < 0 || flat_idx >= static_cast<int>(vals.size())) {
                        throw std::runtime_error("GDML expression error: matrix index out of range for '" + name + "[" + std::to_string(row) + "," + std::to_string(col) + "]' in '" + expr + "'");
                    }
                    result = vals[flat_idx];
                } else {
                    // Single index: name[idx] (1-based)
                    std::string idx_expr = s.substr(idx_start, bracket_end - idx_start);
                    double idx_val = EvalExpression(idx_expr, constants, {}, depth + 1);
                    int idx = static_cast<int>(idx_val + 0.5);
                    int flat_idx = idx - 1;
                    if(flat_idx < 0 || flat_idx >= static_cast<int>(vals.size())) {
                        throw std::runtime_error("GDML expression error: matrix index out of range for '" + name + "[" + std::to_string(idx) + "]' in '" + expr + "'");
                    }
                    result = vals[flat_idx];
                }
                Token t; t.kind = TokenKind::NUMBER; t.value = result;
                tokens.push_back(t);
            } else {
                auto it = constants.find(name);
                if(it == constants.end()) {
                    throw std::runtime_error("GDML expression error: unknown constant '" + name + "' in '" + expr + "'");
                }
                Token t; t.kind = TokenKind::NUMBER; t.value = it->second;
                tokens.push_back(t);
            }
            continue;
        }

        throw std::runtime_error("GDML expression error: unexpected character '" + std::string(1, c) + "' in '" + expr + "'");
    }
    return tokens;
}

double EvalExpression(std::string const & expr, std::unordered_map<std::string, double> const & constants,
                      std::unordered_map<std::string, std::pair<int, std::vector<double>>> const & matrices,
                      int depth) {
    if(depth > 64) {
        throw std::runtime_error("GDML expression error: recursion depth exceeded in '" + expr + "'");
    }
    std::string s = TrimWhitespace(expr);
    if(s.empty()) return 0.0;

    // Fast path: bare number
    {
        char* end = nullptr;
        errno = 0;
        double val = std::strtod(s.c_str(), &end);
        if(end != s.c_str() && end == s.c_str() + s.size()) {
            if(errno == ERANGE)
                throw std::runtime_error("GDML expression error: numeric value out of range in '" + expr + "'");
            return val;
        }
    }

    // Fast path: bare constant name (only if no bracket follows)
    {
        bool bare = true;
        for(size_t i = 0; i < s.size(); ++i) {
            if(!IsIdentChar(s[i])) { bare = false; break; }
        }
        if(bare && IsIdentStart(s[0])) {
            auto it = constants.find(s);
            if(it != constants.end()) return it->second;
            if(matrices.find(s) == matrices.end()) {
                throw std::runtime_error("GDML expression error: unknown constant '" + s + "'");
            }
        }
    }

    std::vector<Token> tokens = Tokenize(s, constants, expr, matrices, depth);

    // Shunting-yard: convert to RPN
    std::vector<Token> output;
    std::vector<Token> op_stack;
    output.reserve(tokens.size());
    op_stack.reserve(tokens.size());

    for(size_t i = 0; i < tokens.size(); ++i) {
        Token const & tok = tokens[i];
        switch(tok.kind) {
        case TokenKind::NUMBER:
            output.push_back(tok);
            break;
        case TokenKind::FUNC:
            op_stack.push_back(tok);
            break;
        case TokenKind::COMMA:
            while(!op_stack.empty() && op_stack.back().kind != TokenKind::LPAREN) {
                output.push_back(op_stack.back());
                op_stack.pop_back();
            }
            if(!op_stack.empty() && op_stack.size() >= 2) {
                Token & func_below = op_stack[op_stack.size() - 2];
                if(func_below.kind == TokenKind::FUNC) func_below.argc++;
            }
            break;
        case TokenKind::UNARY_NEG:
        case TokenKind::OP: {
            // Unary minus is right-associative at the same precedence as ^
            // so -2+3 = (-2)+3 = 1, -2*3 = (-2)*3 = -6, but -2^2 = -(2^2) = -4
            // (matching standard math and Geant4/CLHEP convention).
            int cur_prec = (tok.kind == TokenKind::UNARY_NEG) ? 3 : OpPrecedence(tok.op);
            bool cur_right = (tok.kind == TokenKind::UNARY_NEG) || IsRightAssociative(tok.op);
            while(!op_stack.empty()
                  && (op_stack.back().kind == TokenKind::OP || op_stack.back().kind == TokenKind::UNARY_NEG)) {
                int stack_prec = (op_stack.back().kind == TokenKind::UNARY_NEG) ? 3 : OpPrecedence(op_stack.back().op);
                if(cur_right ? (stack_prec > cur_prec) : (stack_prec >= cur_prec)) {
                    output.push_back(op_stack.back());
                    op_stack.pop_back();
                } else {
                    break;
                }
            }
            op_stack.push_back(tok);
            break;
        }
        case TokenKind::LPAREN:
            op_stack.push_back(tok);
            break;
        case TokenKind::RPAREN:
            while(!op_stack.empty() && op_stack.back().kind != TokenKind::LPAREN) {
                output.push_back(op_stack.back());
                op_stack.pop_back();
            }
            if(op_stack.empty()) {
                throw std::runtime_error("GDML expression error: mismatched parentheses in '" + expr + "'");
            }
            op_stack.pop_back(); // pop LPAREN
            if(!op_stack.empty() && op_stack.back().kind == TokenKind::FUNC) {
                output.push_back(op_stack.back());
                op_stack.pop_back();
            }
            break;
        }
    }
    while(!op_stack.empty()) {
        if(op_stack.back().kind == TokenKind::LPAREN) {
            throw std::runtime_error("GDML expression error: mismatched parentheses in '" + expr + "'");
        }
        output.push_back(op_stack.back());
        op_stack.pop_back();
    }

    // Evaluate RPN
    std::vector<double> val_stack;
    val_stack.reserve(output.size());
    for(auto const & tok : output) {
        switch(tok.kind) {
        case TokenKind::NUMBER:
            val_stack.push_back(tok.value);
            break;
        case TokenKind::OP: {
            if(val_stack.size() < 2) {
                throw std::runtime_error("GDML expression error: malformed expression '" + expr + "'");
            }
            double b = val_stack.back(); val_stack.pop_back();
            double a = val_stack.back(); val_stack.pop_back();
            switch(tok.op) {
            case '+': val_stack.push_back(a + b); break;
            case '-': val_stack.push_back(a - b); break;
            case '*': val_stack.push_back(a * b); break;
            case '/':
                if(b == 0.0) throw std::runtime_error("GDML expression error: division by zero in '" + expr + "'");
                val_stack.push_back(a / b);
                break;
            case '^':
                val_stack.push_back(std::pow(a, b));
                break;
            }
            break;
        }
        case TokenKind::UNARY_NEG: {
            if(val_stack.empty()) {
                throw std::runtime_error("GDML expression error: malformed expression '" + expr + "'");
            }
            val_stack.back() = -val_stack.back();
            break;
        }
        case TokenKind::FUNC: {
            int argc = tok.argc;
            if(argc > 2) {
                throw std::runtime_error("GDML expression error: function '" + tok.name + "' called with " + std::to_string(argc) + " arguments (max 2) in '" + expr + "'");
            }
            if(static_cast<int>(val_stack.size()) < argc) {
                throw std::runtime_error("GDML expression error: not enough arguments for '" + tok.name + "' in '" + expr + "'");
            }
            double args[2];
            for(int j = argc - 1; j >= 0; --j) {
                args[j] = val_stack.back();
                val_stack.pop_back();
            }
            val_stack.push_back(ApplyFunc(tok.name, args, argc, expr));
            break;
        }
        default:
            throw std::runtime_error("GDML expression error: unexpected token in RPN for '" + expr + "'");
        }
    }
    if(val_stack.size() != 1) {
        throw std::runtime_error("GDML expression error: malformed expression '" + expr + "'");
    }
    return val_stack[0];
}

// Parse a double from a string, returning 0 if empty.
// Optionally evaluates expressions using the provided constants map.
double SafeParseDouble(const char* val, std::unordered_map<std::string, double> const & constants,
                       std::unordered_map<std::string, std::pair<int, std::vector<double>>> const & matrices) {
    if(!val || val[0] == '\0') return 0.0;
    char* end = nullptr;
    errno = 0;
    double result = std::strtod(val, &end);
    if(end != val && *end == '\0') {
        if(errno == ERANGE)
            throw std::runtime_error("GDML expression error: numeric value out of range in '" + std::string(val) + "'");
        return result;
    }
    return EvalExpression(std::string(val), constants, matrices);
}


} // namespace gdml
} // namespace detector
} // namespace siren
