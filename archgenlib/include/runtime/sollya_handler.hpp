#ifndef SOLLYA_HANDLER_HPP
#define SOLLYA_HANDLER_HPP

#include <atomic>
#include <memory>
#include <type_traits>

#include <sollya.h>

namespace archgenlib {

namespace detail {
static struct SollyaLifetimeManager {
  SollyaLifetimeManager();
  ~SollyaLifetimeManager();
} sollya_lib_lifetime_manager;
} // namespace detail

class SollyaHandler {
private:
  using manager_t = std::shared_ptr<std::remove_pointer_t<sollya_obj_t>>;
  manager_t manager;

public:
  SollyaHandler(sollya_obj_t);
  operator sollya_obj_t() { return manager.get(); };
};
} // namespace archgenlib

#endif
