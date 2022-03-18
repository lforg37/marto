#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include "fixedpoint/fixedpoint.hpp"
#include "runtime/expression_tree.hpp"
#include "runtime/fixfunction_multipartite.hpp"
#include "runtime/operator_builder/table_builder.hpp"
#include "runtime/sollya_fix_function.hpp"

__attribute((used)) void print_mpz_class(mpz_class const &val) {
  std::cerr << val << "\n";
}

namespace {
using namespace archgenlib;
struct BinStatsHolder {
  mpz_class *min;
  mpz_class *max;
  mpz_class scaled_init_value;
};

struct GammaBinStatsHolder {
  struct ElementWiseStats {
    mpz_class scaled_min;
    mpz_class scaled_max;
  };

  std::vector<ElementWiseStats> elementwise_stats;
};

void init_bin_stats(auto &values, std::vector<BinStatsHolder> &stats) {
  for (size_t i = 0; i < (values.size() >> 1); ++i) {
    auto &low = values[2 * i];
    auto &high = values[2 * i + 1];
    bool low_is_biggest = low > high;
    auto &min = low_is_biggest ? high : low;
    auto &max = low_is_biggest ? low : high;
    auto scaled_init_value = max + min;
    stats.emplace_back(
        BinStatsHolder{&min, &max, mpz_class{scaled_init_value}});
  }
}

void merge_stats_2by2(auto &stats) {
  for (size_t i = 0; i < (stats.size() >> 1); ++i) {
    auto low = stats[2 * i];
    auto high = stats[2 * i + 1];
    auto min = *(low.min) < *(high.min) ? low.min : high.min;
    auto max = *(low.max) > *(high.max) ? low.max : high.max;
    auto &dest = stats[i];
    dest.min = min;
    dest.max = max;
    auto init_value = *max + *min;
    dest.scaled_init_value = init_value;
  }
  stats.resize(stats.size() >> 1);
}

bool init_gamma_elements(auto &values, uint64_t beta,
                         mpz_class const &error_bound, auto &stats,
                         std::vector<GammaBinStatsHolder> &bin_value) {
  // 1 for as everything is scaled, +1 as we are interested in the error for
  // half of the difference
  mpz_class scaled_bound{error_bound << 1};

  uint64_t group_by{1};
  group_by <<= beta;

  std::vector<GammaBinStatsHolder> ret;
  // We always start with gamma = alpha - 1 so
  // we have two tiv value per gamma bin
  std::uint64_t n_gamma_groups = values.size() >> (beta + 1);
  ret.reserve(n_gamma_groups);
  std::uint64_t loop_limit = n_gamma_groups;

  for (std::uint64_t i = 0; i < loop_limit; ++i) {
    std::vector<GammaBinStatsHolder::ElementWiseStats> bin_val{};
    bin_val.reserve(1 << (beta - 1));
    auto leftbin_idx = i << (beta + 1);
    auto rightbin_idx = leftbin_idx + (1 << beta);
    std::uint64_t inside_loop_limit{1};
    inside_loop_limit <<= beta;
    for (std::size_t j = 0; j < inside_loop_limit; ++j) {
      mpz_class left{(values[leftbin_idx + j] << 1) -
                     stats[2 * i].scaled_init_value};
      mpz_class right{(values[rightbin_idx + j] << 1) -
                      stats[2 * i + 1].scaled_init_value};
      auto max = left > right ? left : right;
      auto min = left < right ? left : right;

      if (max - min >= scaled_bound) {
        return false;
      }

      bin_val.emplace_back(GammaBinStatsHolder::ElementWiseStats{min, max});
    }
    ret.emplace_back(GammaBinStatsHolder{bin_val});
  }
  std::swap(ret, bin_value);
  return true;
}

auto merge_gamma_bin_2x2(auto &bin_values, mpz_class const &error_bound) {
  mpz_class scaled_bound{error_bound << 1};
  auto bin_nb_elem = bin_values[0].elementwise_stats.size();
  auto nb_new_bin = bin_values.size() >> 1;
  for (size_t i = 0; i < nb_new_bin; ++i) {
    auto &left = bin_values[2 * i];
    auto &right = bin_values[2 * i + 1];
    std::vector<GammaBinStatsHolder::ElementWiseStats> res;
    for (size_t j = 0; j < bin_nb_elem; ++j) {
      auto &lelem = left.elementwise_stats[j];
      auto &relem = right.elementwise_stats[j];
      auto &max = (lelem.scaled_max > relem.scaled_max) ? lelem.scaled_max
                                                        : relem.scaled_max;
      auto &min = (lelem.scaled_min < relem.scaled_min) ? lelem.scaled_min
                                                        : relem.scaled_min;
      assert(max >= min);
      if (max - min >= scaled_bound) {
        return false;
      }
      res.emplace_back(GammaBinStatsHolder::ElementWiseStats{min, max});
    }
    std::swap(res, bin_values[i].elementwise_stats);
  }
  return true;
}

} // namespace

namespace archgenlib {

MultipartiteFunction::OffsetConfigGenerator::OffsetConfigGenerator(
    std::size_t nb_tables, std::size_t total_width)
    : total_size{total_width}, config(nb_tables, 1) {
  assert(nb_tables > 0);
  assert(nb_tables <= total_size);
  config[0] = total_size - nb_tables + 1;
}

bool MultipartiteFunction::OffsetConfigGenerator::has_next() const {
  return config.back() != total_size - config.size() + 1;
}

void MultipartiteFunction::OffsetConfigGenerator::next() {
  //[3, 4, 1, 2] corresponds to this segmentation:
  // [xxx|xxxx|x|xx]
  //          ^ first idx that can move

  int reset_from = -1;
  for (std::size_t i = config.size() - 1; i > 0; --i) {
    // We find the first index for which we can move the boundary left
    if (config[i + 1] > 1) {
      config[i + 1] -= 1;
      reset_from = i;
      break;
    }
  }
  //
  // After exiting the previous loop we have
  // [3, 3, 1, 2] and reset_idx = 2

  if (reset_from != -1) {
    // reset from point to the space right to the one that was
    // decremented. We had one to compensate the left space
    // reduction
    vecwidth_t sum_reset{1};

    // Then shrink all the spaces on the right to reset the
    // sub space exploration
    for (size_t i = 1 + reset_from; i < config.size(); ++i) {
      // And collect the removed space
      sum_reset += config[i] - 1;
      config[i] = 1;
    }
    config[reset_from] += sum_reset;
  }

  // No we get [3, 3, 3, 1]
  // Corresponding to [xxx|xxx|xxx|x]
}

MultipartiteFunction::MultipartiteConfiguration::MultipartiteConfiguration(
    SollyaFunction &func, vecwidth_t alpha, std::size_t guard_bits,
    std::vector<OffsetTableConfig> &&offset_config, std::size_t extra_precision,
    std::vector<mpz_class> &precise_values)
    : alpha{alpha}, guard_bits{guard_bits}, offset_configs{offset_config},
      initial_values_table{} {

  vecwidth_t input_width = func.input_format.width;

  std::vector<BinStatsHolder> stats{};
  init_bin_stats(precise_values, stats);
  for (std::size_t i = 0; i < input_width - alpha - 1; ++i) {
    merge_stats_2by2(stats);
  }

  mpz_class round_val{1};
  // The init tiv value is scaled by 2, so
  // extra_precision + 1 will reach the bit that corresponds to LSBout
  round_val <<=
      extra_precision - guard_bits; // reaches below guard bit (round bit)
  initial_values_table.reserve(stats.size());

  mpz_class truncation_error_bound{1};
  // LSBout -2 but compared to scaled value
  truncation_error_bound <<= extra_precision - 1;

  for (auto &stat : stats) {
    mpz_class scaled_res{(stat.scaled_init_value + round_val - 1) >>
                         (extra_precision - guard_bits + 1)};

    initial_values_table.emplace_back(scaled_res); // TODO ADD one ?
  }

  // TODO: For now we have only one table of offset
  auto beta = offset_config[0].beta;
  auto gamma = offset_config[0].gamma;
  std::vector<mpz_class> to{};
  round_val <<= 1;
  // Will be compared with a value expressed in quarters
  truncation_error_bound <<= 1;

  for (std::size_t c = 0; c < (1 << gamma); ++c) {
    for (std::size_t b = 0; b < (1 << beta); ++b) {
      auto tiv_idx = c << (alpha - gamma);
      auto value_idx = (c << (input_width - gamma)) | b;
      mpz_class scaled_min =
          (precise_values[value_idx] << 1) - stats[tiv_idx].scaled_init_value;
      mpz_class scaled_max = scaled_min;
      for (std::size_t i = 1; i < (1 << (alpha - gamma)); ++i) {
        ++tiv_idx;
        value_idx += 1 << beta;
        auto scaled_offset = (precise_values[value_idx] << 1) -
                             (stats[tiv_idx].scaled_init_value);
        if (scaled_offset > scaled_max)
          scaled_max = scaled_offset;
        if (scaled_offset < scaled_min)
          scaled_min = scaled_offset;
      }
      mpz_class scaled_mean = scaled_max + scaled_min;
      mpz_class mean = (scaled_mean + round_val - 1) >> (extra_precision - guard_bits + 2);
      to.emplace_back(mean);
      mpz_class rescaled_to = to.back() << (extra_precision - guard_bits + 2);

      // Error computation :
      //   err = f - (Tiv + To)
      //       = f - fapprox + fapprox - (Tiv* + To*) + (Tiv* + To*) - (Tiv + To)
      // 
      // |err| <= |f - fapprox| + |fapprox - (Tiv* + To*)| + |(Tiv* + To*) - (Tiv + To)|
      //
      // We want to ensure that |err| < 2 ^ lsbout
      //
      //  |f - fapprox| < 2^(lsbout - extra_precision) by precision requirement on Sollya
      //  
      // Let's ensure that: 
      //  -) |fapprox - (Tiv* + To*)| < 2^(lsbout - 1) - 2^(lsbout - extra_precision)
      //  -) |Tiv* + To* - Tiv - To|  < 2^(lsbout - 1)
      {
        mpz_class allowed_error{1};
        // LSBOut - 1, but all comparison expressed as i * 2^(lsbout - extra_precision - 2)
        allowed_error <<= (extra_precision - 1 + 2);  

        mpz_class approx_error = allowed_error - (2 << 2); // Take care of approximation function 

        auto tiv_idx = c << (alpha - gamma);
        auto value_idx = (c << (input_width - gamma)) | b;
        for (std::size_t i = 0; i < (1 << (alpha - gamma)); ++i) {
          mpz_class scaled_reconstruction = // Tiv* + To* 
              scaled_mean + (stats[tiv_idx].scaled_init_value << 1);
          mpz_class scaled_ref = precise_values[value_idx] << 2;
          mpz_class diff = abs(scaled_ref - scaled_reconstruction);
          // Check that |TIV* + TO* - fapprox| < 2^(lsbout - 1) - approx_err
          assert(diff < approx_error);

          // Check that |TIV* + TO* - TI - TO| < 2^(lsbout - 1)
          mpz_class mp_reconstruction = ((initial_values_table[tiv_idx]) << (extra_precision - guard_bits + 2)) + rescaled_to;
          diff = abs(scaled_reconstruction - mp_reconstruction);
          if ( diff >= allowed_error) {
            std::cout << "tiv*: " << (stats[tiv_idx].scaled_init_value << 1)
                    << "\ntiv: " << (initial_values_table[tiv_idx] << (extra_precision - guard_bits + 2))
                    << "\nto*: " << scaled_mean 
                    << "\nto: " << rescaled_to
                    << "\nscaled_reconstruction: " << scaled_reconstruction
                    << "\nmp_reconstruction: " << mp_reconstruction
                    << "\ndiff: " << diff 
                    << "\nextra_prec: " << extra_precision
                    << "\nallowed_error: " << allowed_error << std::endl;
          }
          assert( diff < allowed_error);


          // Check all except truncation
          diff = abs(scaled_ref - mp_reconstruction);
          assert(diff < (2*allowed_error - 8));

          ++tiv_idx;
          value_idx += 1 << beta;
        }
      }

    }
  }
  // TODO Embed the rounding bit inside the TIV TODO


  offset_tables.emplace_back(to);
}

void MultipartiteFunction::find_best_config() {
  vecwidth_t required_output_width = function.msb_output - lsb_out + 1;
  auto input_width = function.input_format.width;
  uint64_t best_cost{required_output_width};
  best_cost <<= input_width;
  vecwidth_t best_alpha{input_width};
  vecwidth_t best_gamma{0};
  std::size_t extra_Prec = 3 * input_width;
  auto precise_value = function.faithful_at_weight(lsb_out - extra_Prec);
  mpz_class allowed_error{1};
  allowed_error <<= extra_Prec - 2; // TODO
  allowed_error -= 3;
  std::vector<BinStatsHolder> stats;
  for (vecwidth_t alpha = input_width - 1; alpha > 0; --alpha) {
    vecwidth_t beta = input_width - alpha;
    if (alpha == input_width - 1) {
      init_bin_stats(precise_value, stats);
    } else {
      merge_stats_2by2(stats);
    }
    std::vector<GammaBinStatsHolder> gamma_stats{};
    for (vecwidth_t gamma = alpha - 1; gamma > 0; --gamma) {
      if (gamma == alpha - 1) {
        if (!init_gamma_elements(precise_value, beta, allowed_error, stats,
                                 gamma_stats))
          break;
      } else {
        if (!merge_gamma_bin_2x2(gamma_stats, allowed_error))
          break;
      }
      MultipartiteConfiguration{function,   alpha,
                                1,          {{gamma, input_width - alpha, 0}},
                                extra_Prec, {precise_value}};
      std::uint64_t TIVCost{required_output_width + 2};
      TIVCost <<= alpha;
      std::uint64_t TOCost{required_output_width + 2};
      TOCost <<= gamma + beta;
      auto total_cost = TIVCost + TOCost;
      if (total_cost < best_cost) {
        best_cost = total_cost;
        best_alpha = alpha;
        best_gamma = gamma;
      }
    }
  }

  std::optional<MultipartiteConfiguration> ret;

  if (best_alpha < input_width) {
    best_config = MultipartiteConfiguration(
        function, best_alpha, 2, {{best_gamma, input_width - best_alpha, 0}},
        extra_Prec, {precise_value});
  }
};

MultipartiteFunction::MultipartiteFunction(SollyaFunction &func,
                                           bitweight_t LSBOut)
    : max_sub_tables{7}, function{func}, lsb_out{LSBOut} {
  find_best_config();
}

bool MultipartiteFunction::check_best_config(bitweight_t prec) const {
  if (best_config.has_value()) {
    assert(prec < lsb_out);
    auto reference = function.faithful_at_weight(prec);
    mpz_class error_bound{1};
    error_bound <<= (lsb_out - prec);
    error_bound += 1;

    uint64_t nb_inputs = uint64_t{1} << function.input_format.width;

    auto alpha = best_config->alpha;
    auto beta = best_config->offset_configs[0].beta;
    auto gamma = best_config->offset_configs[0].gamma;
    auto const &tiv = best_config->initial_values_table;
    auto const &to = best_config->offset_tables[0];
    auto input_width = function.input_format.width;
    auto beta_mask = (1 << beta) - 1;

    for (std::uint64_t i = 0; i < nb_inputs; ++i) {
      mpz_class iv = tiv[i >> beta];
      auto offset_idx =
          ((i >> (input_width - gamma)) << beta) | (i & beta_mask);
      mpz_class offset = to[offset_idx];
      mpz_class res = (iv + offset + 1) >> 2; // 1 for the guard bit
      mpz_class scaled_res = res << (lsb_out - prec);
      mpz_class diff = abs(scaled_res - reference[i]);
      assert(diff < error_bound);
    }
    return true;
  } else {
    return true;
  }
}

} // namespace archgenlib
