#include "flagalg/spec.hpp"

#include <cctype>
#include <sstream>
#include <stdexcept>

namespace flagalg {

mpq_class parseRational(const std::string& s) {
    mpq_class q;
    if (q.set_str(s, 10) != 0) throw std::invalid_argument("bad rational: " + s);
    q.canonicalize();
    return q;
}

namespace {

class Parser {
public:
    explicit Parser(const std::string& t) : t_(t) {}

    std::shared_ptr<Expr> expr() {
        skip();
        expect('(');
        std::string op = atom();
        auto e = std::make_shared<Expr>();
        e->op = op;
        if (op == "flag") {
            int k = std::stoi(atom()), s = std::stoi(atom());
            std::vector<int> cols(numPairs(k));
            for (int& c : cols) c = std::stoi(atom());
            e->flag = makeFlag(k, s, cols);
        } else if (op == "scale") {
            e->scalar = parseRational(atom());
            e->args.push_back(expr());
        } else if (op == "add" || op == "mul" || op == "avg") {
            while (true) {
                skip();
                if (peek() == ')') break;
                e->args.push_back(expr());
            }
            if (e->args.empty()) throw std::invalid_argument(op + " needs arguments");
            if (op == "avg" && e->args.size() != 1) throw std::invalid_argument("avg takes one argument");
        } else {
            throw std::invalid_argument("unknown operator: " + op);
        }
        skip();
        expect(')');
        return e;
    }

    std::string atom() {
        skip();
        std::size_t start = pos_;
        while (pos_ < t_.size() && !std::isspace(static_cast<unsigned char>(t_[pos_])) && t_[pos_] != '(' &&
               t_[pos_] != ')')
            ++pos_;
        if (start == pos_) throw std::invalid_argument("expected an atom");
        return t_.substr(start, pos_ - start);
    }

    void skip() {
        while (pos_ < t_.size() && std::isspace(static_cast<unsigned char>(t_[pos_]))) ++pos_;
    }
    bool done() {
        skip();
        return pos_ == t_.size();
    }

private:
    char peek() const { return pos_ < t_.size() ? t_[pos_] : '\0'; }
    void expect(char c) {
        if (peek() != c) throw std::invalid_argument(std::string("expected '") + c + "'");
        ++pos_;
    }
    const std::string& t_;
    std::size_t pos_ = 0;
};

std::shared_ptr<Expr> parseWhole(const std::string& text) {
    Parser p(text);
    auto e = p.expr();
    if (!p.done()) throw std::invalid_argument("trailing text after expression: " + text);
    return e;
}

}  // namespace

Spec parseSpec(const std::string& text) {
    Spec spec;
    std::istringstream in(text);
    std::string line;
    int lineno = 0;
    while (std::getline(in, line)) {
        ++lineno;
        auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        std::istringstream ls(line);
        std::string key;
        if (!(ls >> key)) continue;
        std::string rest;
        std::getline(ls, rest);
        try {
            if (key == "colors") {
                spec.colors = std::stoi(rest);
            } else if (key == "n") {
                spec.n = std::stoi(rest);
            } else if (key == "blowup") {
                std::istringstream ls(rest);
                std::string mode;
                ls >> spec.blowup >> mode;
                if (spec.blowup < 1) throw std::invalid_argument("blowup N must be positive");
                if (mode == "levels") spec.blowupLevels = true;
                else if (!mode.empty()) throw std::invalid_argument("unknown blowup mode " + mode);
            } else if (key == "empty_type") {
                spec.emptyType = std::stoi(rest) != 0;
            } else if (key == "zero") {
                auto e = parseWhole(rest);
                if (e->op != "flag") throw std::invalid_argument("zero must be a flag literal");
                spec.zeros.push_back(e->flag.forgetRoots());
            } else if (key == "objective") {
                spec.objective = parseWhole(rest);
            } else if (key == "eq" || key == "le") {
                std::istringstream rs(rest);
                std::string q;
                rs >> q;
                std::string exprText;
                std::getline(rs, exprText);
                spec.constraints.push_back(
                    {key == "eq" ? Constraint::Eq : Constraint::Le, parseRational(q), parseWhole(exprText)});
            } else if (key == "psd") {
                std::istringstream rs(rest);
                int m;
                rs >> m;
                std::string exprText;
                std::getline(rs, exprText);
                spec.localizing.push_back({m, parseWhole(exprText)});
            } else {
                throw std::invalid_argument("unknown directive " + key);
            }
        } catch (const std::exception& ex) {
            throw std::invalid_argument("spec line " + std::to_string(lineno) + ": " + ex.what());
        }
    }
    if (spec.n < 1) throw std::invalid_argument("spec: n must be given");
    if (spec.n > kMaxVertices) throw std::invalid_argument("spec: n exceeds kMaxVertices");
    for (const auto& L : spec.localizing)
        if (L.m > kMaxVertices) throw std::invalid_argument("spec: psd flag size exceeds kMaxVertices");
    return spec;
}

Combo evaluate(const Theory& T, const Expr& e) {
    if (e.op == "flag") {
        Combo c;
        addTerm(c, T.canon(e.flag), 1);
        return c;
    }
    if (e.op == "scale") {
        Combo c = evaluate(T, *e.args[0]);
        Combo out;
        for (auto& [f, v] : c) addTerm(out, f, v * e.scalar);
        return out;
    }
    if (e.op == "add") {
        Combo out;
        for (auto& a : e.args)
            for (auto& [f, v] : evaluate(T, *a)) addTerm(out, f, v);
        return out;
    }
    if (e.op == "mul") {
        Combo out = evaluate(T, *e.args[0]);
        for (std::size_t i = 1; i < e.args.size(); ++i) out = T.multiply(out, evaluate(T, *e.args[i]));
        return out;
    }
    if (e.op == "avg") return T.average(evaluate(T, *e.args[0]));
    throw std::invalid_argument("evaluate: unknown op " + e.op);
}

}  // namespace flagalg
