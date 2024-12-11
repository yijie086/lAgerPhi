// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>
// 
// This file is part of lAger.
// 
// lAger is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Shoftware Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// lAger is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with lAger.  If not, see <https://www.gnu.org/licenses/>.
// 

#ifndef LAGER_CORE_FACTORY_LOADED
#define LAGER_CORE_FACTORY_LOADED

#include <map>
#include <memory>
#include <mutex>
#include <lager/core/assert.hh>
#include <lager/core/logger.hh>
#include <string>

namespace lager {

// errors
class factory_error : public lager::exception {
public:
  factory_error(const std::string& msg,
                const std::string& type = "factory_error")
      : lager::exception{msg, type} {}
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
    proxy_##class() { gen::factory_instance.add<class>(name); }                \
  };                                                                           \
  const proxy_##class p_##class;                                               \
  }
#define FACTORY_REGISTER2(gen, class, name)                                    \
  gen::factory_instance.add<class>(name)

// construct a generator instance
#define FACTORY_CREATE(gen, conf, path, rng)                                   \
  gen::factory_instance.create(                                                \
      conf.get<std::string>(string_path(path) / "type"), conf, path, rng)

// a generic factory base class
template <class T, class... Args> class factory {
private:
  struct worker_base;
  using worker_map_type = std::map<std::string, std::shared_ptr<worker_base>>;

public:
  std::shared_ptr<T> create(const std::string& name, Args... a) const {
    // Removed lock guard because it screws up when doing nested constructions,
    // and we are not running in a multi-threaded environment anyway
    //    std::lock_guard<std::mutex> lock{mutex_};
    LOG_DEBUG("factory", "HK Constructing " + name);
    if (!workers_) {
      LOG_DEBUG("factory",
                "Worker pointer is a null pointer. This should never happen.");
    }
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
    LOG_DEBUG("factory", "Registered class: " + name);
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

} // namespace lager

#endif
