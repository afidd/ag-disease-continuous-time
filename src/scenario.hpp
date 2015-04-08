#ifndef _SCENARIO_HPP_
#define _SCENARIO_HPP_ 1

#include <cstddef>
#include "spatial_process.hpp"

template<typename RNG>
class Scenario {
 public:
  Scenario(size_t herd_cnt, RNG& rng) : herd_cnt_(herd_cnt) {
    InitializeSpatial(rng);
  }
  ~Scenario() {};
  size_t herd_cnt() { return herd_cnt_; }
  double airborne_hazard(int64_t i, int64_t j) {
    return 1./distance2(i, j);
  }
  double distance2(int64_t i, int64_t j) {
    auto& ixy=locations_[i];
    auto& jxy=locations_[j];
    return std::pow(std::get<0>(ixy)-std::get<0>(jxy), 2) +
        std::pow(std::get<1>(ixy)-std::get<1>(jxy), 2);
  }
 private:
  void InitializeSpatial(RNG& rng) {
    locations_=hard_sphere_process({0., 1., 0., 1.}, herd_cnt_, 0.2, rng);
  }
  std::vector<std::array<double,2>> locations_;
  size_t herd_cnt_;
};


#endif // _SCENARIO_HPP_
