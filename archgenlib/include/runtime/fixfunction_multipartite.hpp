#ifndef RUNTIME_FIXFUNCTION_MULTIPARTITE_HPP
#define RUNTIME_FIXFUNCTION_MULTIPARTITE_HPP

#include <cstddef>
#include <gmpxx.h>
#include <vector>

#include "fixedpoint/fixedpoint.hpp"
#include "sollya_fix_function.hpp"

namespace archgenlib {
class MultipartiteFunction {
public:
  struct OffsetConfigGenerator {
  private:
    std::size_t total_size;
    std::vector<vecwidth_t> config;

  public:
    OffsetConfigGenerator(std::size_t nb_tables, std::size_t total_width);
    bool has_next() const;
    void next();
    std::vector<vecwidth_t> const &getConfig() { return config; }
  };

  struct ErrorTracker {
    std::uint64_t relError;
    std::uint32_t shiftComparedToLSB;
  };

  struct MultipartiteConfiguration {
    struct OffsetTableConfig {
      vecwidth_t gamma; ///< Number of MSB to address this table
      vecwidth_t beta;  ///< Number of LSB to address this table
      vecwidth_t p;     ///< Offset from LSB to address
    };
    vecwidth_t alpha;       ///< Size of the TIV address
    std::size_t guard_bits; ///< Number of required guard bits

    /// Stored in increasing order
    std::vector<OffsetTableConfig> offset_configs;
    std::vector<mpz_class> initial_values_table;
    std::vector<std::vector<mpz_class>> offset_tables;

    MultipartiteConfiguration(SollyaFunction &func, vecwidth_t alpha,
                              std::size_t guard_bits,
                              std::vector<OffsetTableConfig> &&offset_config,
                              std::size_t extra_precision,
                              std::vector<mpz_class> &precise_values);
  };

private:
  using error_t = ErrorTracker;
  vecwidth_t max_sub_tables;

  MultipartiteConfiguration find_best_config();

public:
  SollyaFunction &function;
  bitweight_t lsb_out;

  MultipartiteFunction(SollyaFunction &func, bitweight_t LSBOut);
  bool check_best_config(bitweight_t prec) const;

private:
  MultipartiteConfiguration best_config;
};
} // namespace archgenlib

#endif
