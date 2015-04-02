#ifndef _SCENARIO_HPP_
#define _SCENARIO_HPP_ 1

#include <cstddef>

class Scenario {
 public:
  Scenario(size_t herd_cnt);
  ~Scenario();
  size_t size();
 private:
  size_t herd_cnt_;
};


#endif // _SCENARIO_HPP_
