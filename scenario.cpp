#include "scenario.hpp"


Scenario::Scenario(size_t herd_cnt) : herd_cnt_(herd_cnt) {}
Scenario::~Scenario() {}

size_t Scenario::size() { return herd_cnt_; }
