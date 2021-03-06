#ifndef _SIR_EXP_HPP_
#define _SIR_EXP_HPP_ 1
#include <map>
#include <random>
#include "boost/random/mersenne_twister.hpp"
#include "naadsm_xml.hpp"
//#include "mt19937.hpp"

//using RandGen=afidd::rng::mt19937;
//using RandGen=boost::mt19937;
using RandGen=std::mt19937;

enum class ADParam { FirstFarm };
enum class DiseaseState : int { None, Susceptible, Latent,
  Subclinical, Clinical, Recovered, Immune, Vaccinated, Dead };

struct Parameter {
  ADParam kind;
  std::string name;
  double value;
  std::string description;
  Parameter(ADParam k, std::string n, double v, std::string desc)
  : kind(k), name(n), value(v), description(desc) {}
  Parameter()=default;
};

struct TrajectoryEntry {
  int64_t s;
  int64_t i;
  int64_t r;
  double t;
  TrajectoryEntry(int64_t s, int64_t i, int64_t r, double t)
  : s(s), i(i), r(r), t(t) {}
  TrajectoryEntry()=default;
};

class TrajectoryObserver
{
public:
  virtual void Step(TrajectoryEntry sirt)=0;
  virtual const std::vector<TrajectoryEntry>& Trajectory() const =0;
};


int64_t SIR_run(double time_limit,const std::vector<Parameter>& parameters,
  NAADSMScenario& scenario, std::shared_ptr<TrajectoryObserver> observer,
  RandGen& rng);

#endif
