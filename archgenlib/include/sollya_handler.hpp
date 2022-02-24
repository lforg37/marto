#ifndef SOLLYA_HANDLER_HPP
#define SOLLYA_HANDLER_HPP

#include <atomic>
#include <memory>
#include <type_traits>

#include <sollya.h>

namespace archgenlib {
class SollyaHandler {
private:
  sollya_obj_t managed;
  std::atomic_uint32_t *ref_count = nullptr;
  static struct LibLifetimeManager {
    LibLifetimeManager() { sollya_lib_init(); }
    ~LibLifetimeManager() { sollya_lib_close(); }
  } lifetime_manager;

  void dec_ref_count() {
    if (ref_count != nullptr) {
      auto ret = --(*ref_count);
      if (ret == 0) {
        delete ref_count;
        sollya_lib_clear_obj(managed);
      }
    }
  }
  void inc_ref_count(){++(*ref_count);}

public:
  SollyaHandler(sollya_obj_t so) : managed{so} {
    ref_count = new std::atomic_uint32_t{1};
  }

  SollyaHandler(SollyaHandler const &orig) {
    managed = orig.managed;
    ref_count = orig.ref_count;
    inc_ref_count();
  }

  SollyaHandler(SollyaHandler const &&orig) {
    managed = orig.managed;
    ref_count = orig.ref_count;
    inc_ref_count();
  }

  ~SollyaHandler() { dec_ref_count(); }

  SollyaHandler &operator=(SollyaHandler const &orig) {
      dec_ref_count();
      managed = orig.managed;
      ref_count = orig.ref_count;
      inc_ref_count();
      return *this;
  }
};
} // namespace archgenlib

#endif
