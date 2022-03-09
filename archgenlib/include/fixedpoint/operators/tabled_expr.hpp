
#include <ostream>

namespace archgenlib {

template<typename T>
struct table_emitter {
  const T& vec;
  friend std::ostream& operator<<(std::ostream& os, table_emitter<T> v) {
    for (size_t idx = 0; idx < v.vec.size(); idx++) {
      if (idx != 0)
        os << ", ";
      os << v.vec[idx];
    }
    return os;
  }
};

template<typename T>
table_emitter<T> emit_table(const T& vec) {
  return {vec};
}

template<typename T>
struct tabled_expr_emitter {
  const T& vec;
  const char* name;
  friend std::ostream &operator<<(std::ostream &os, tabled_expr_emitter<T> v) {
    os << "auto " << v.name << " = [] (auto e) {\n";
    os << "constexpr int table[] = {" << emit_table(v.vec) << "};\n";
    os << "return table[e];\n";
    os << "};\n";
    return os;
  }
};

template<typename T>
tabled_expr_emitter<T> emit_tabled_expr(const char* name, const T& vec) {
  return {vec, name};
}

}
