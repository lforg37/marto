
#include <compare>

#include <sycl/sycl.hpp>
#include "fixedpoint/operators.hpp"
#include "fixedpoint/fixedpoint.hpp"
#include "fixedpoint/literal.hpp"
#include "fixedpoint/evaluate.hpp"

int main() {
  double val = 0.5;
  using format_t = archgenlib::FixedFormat<1, -8, signed>;
  using fixed_t = archgenlib::FixedNumber<format_t>;
  std::vector<fixed_t> d;
  d.push_back(fixed_t::get_from_value(val));
  sycl::buffer<fixed_t> a{d.data(), {d.size()}};
  sycl::buffer<fixed_t> b{d.size()};

  sycl::queue q;
  q.submit([&](sycl::handler &cgh) {
    sycl::accessor acc_a{a, cgh, sycl::read_only};
    sycl::accessor acc_b{b, cgh, sycl::write_only};
    cgh.single_task([=] {
      auto var = archgenlib::FreeVariable(acc_a[0]);
      acc_b[0] = archgenlib::evaluate<format_t>(archgenlib::sin(var * archgenlib::pi / 0x2p0_cst));
    });
  });
  {
    sycl::host_accessor acc_b{b, sycl::read_only};
    std::cout << acc_b[0].get_as<double>() << std::endl;
  }
}
