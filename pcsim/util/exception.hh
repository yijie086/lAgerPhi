#ifndef PYTHIA6M_UTIL_EXCEPTION_LOADED
#define PYTHIA6M_UTIL_EXCEPTION_LOADED

#include <exception>
#include <string>

namespace pythia6m {
class exception : public std::exception {
public:
  exception(const std::string& msg, const std::string& type = "exception")
      : msg_{msg}, type_{type} {}

  virtual const char* what() const throw() { return msg_.c_str(); }
  virtual const char* type() const throw() { return type_.c_str(); }
  virtual ~exception() throw() {}

private:
  std::string msg_;
  std::string type_;
};
} // ns pythia6m

#endif
