#ifndef PCSIM_CORE_FACTORY_LOADED
#define PCSIM_CORE_FACTORY_LOADED

#include <map>
#include <memory>
#include <mutex>
#include <pcsim/core/assert.hh>
#include <pcsim/core/logger.hh>
#include <string>

namespace pcsim {

// errors
class factory_error : public pcsim::exception {
public:
  factory_error(const std::string& msg,
                const std::string& type = "factory_error")
      : pcsim::exception{msg, type} {}
};
class factory_key_error : public factory_error {
public:
  factory_key_error(const std::string& msg)
      : factory_error{msg, "factory_key_error"} {}
};

// register a class with a (global) factory
#define FACTORY_REGISTER(gen, class, name)                                     \
  namespace {                                                                  \
  struct proxy_##class {                                                       \
    proxy_##class() { gen::factory.add<class>(name); }                         \
  };                                                                           \
  const proxy_##class p_##class;                                               \
  }

#define FACTORY_CREATE(gen, conf, path, rng)                                   \
  gen::factory.create(conf.get<std::string>(path "/type"), conf, path, rng)

// a generic factory base class
template <class T, class... Args> class factory {
private:
  struct worker_base;
  using worker_map_type = std::map<std::string, std::shared_ptr<worker_base>>;

public:
  std::shared_ptr<T> create(const std::string& name, Args... a) const {
    std::lock_guard<std::mutex> lock{mutex_};
    LOG_DEBUG("factory", "Constructing " + name);
    if (!workers_ || !workers_->count(name)) {
      LOG_ERROR("factory", "Cannot construct object " + name +
                               ": no such name registered.");
      throw factory_key_error{"Failed to construct object."};
    }
    return workers_->find(name)->second->create(a...);
  }

  template <class U> void add(const std::string& name) {
    // there is NO guarantee that the worker map will be initialized before
    // we attempt to proxy-register a class to a global factory.
    // To avoid this problem we manually construct the map when it's
    // first needed
    // A bit of a hack, but necessecary
    std::lock_guard<std::mutex> lock{mutex_};
    if (!workers_) {
      workers_ = std::make_unique<worker_map_type>();
    }
    workers_->insert(std::make_pair(name, std::make_shared<worker<U>>()));
  }

  void remove(const std::string& name) {
    std::lock_guard<std::mutex> lock{mutex_};
    if (!workers_ || !workers_->count(name)) {
      LOG_ERROR("factory", "Cannot remove worker for class " + name +
                               ": no such name registered.");
      throw factory_key_error{"Failed to remove worker."};
    }
    workers_.erase(name);
  }

private:
  // factory worker
  struct worker_base {
    virtual std::shared_ptr<T> create(Args...) const = 0;
  };
  template <class U> struct worker : worker_base {
    virtual std::shared_ptr<T> create(Args... a) const {
      return std::make_shared<U>(a...);
    }
  };

  std::unique_ptr<worker_map_type> workers_;
  mutable std::mutex mutex_;

};

}

#endif
