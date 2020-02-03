#include "composite.hh"
#include <liege/core/factory.hh>

namespace liege {
namespace detector {

composite::composite(const configuration& cf, const string_path& path,
                     std::shared_ptr<TRandom> r)
    : composite::base_type{r} {
  auto tmp_conf = cf;
  auto& conf = tmp_conf.raw_node(path / "components");
  for (const auto& child : conf) {
    string_path child_path = path / "components" / child.first.c_str();
    LOG_INFO("composite", "Constructing child detector: " + child_path.str());
    auto det = FACTORY_CREATE(detector, tmp_conf, child_path, r);
    detectors_.push_back(det);
  }
}
} // namespace detector

} // namespace liege
