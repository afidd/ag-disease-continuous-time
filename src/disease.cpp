
#include <tuple>
#include <map>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <memory>
#include <set>
#include <functional>
#include "stochnet.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/log/core.hpp"
#include "boost/math/constants/constants.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/mpl/vector.hpp"
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "disease.hpp"


namespace smv=afidd::smv;
using namespace smv;


struct SimpleToken
{
  SimpleToken()=default;

  inline friend
  std::ostream& operator<<(std::ostream& os, const SimpleToken& it){
    return os << "T";
  }
};


struct HerdToken
{
  DiseaseState disease_state;
  int64_t herd_size;

  HerdToken()=default;
  HerdToken(DiseaseState disease, int64_t cnt)
    : disease_state(disease), herd_size(cnt) {}

  inline friend
  std::ostream& operator<<(std::ostream& os, const HerdToken& it){
    return os << "T";
  }
};


// Animal disease place
struct ADPlace
{
  int64_t location;
  int64_t kind;
  ADPlace()=default;
  ADPlace(int64_t l, int64_t k) : location(l), kind(k) {}
  friend inline
  bool operator<(const ADPlace& a, const ADPlace& b) {
    return LazyLess(a.location, b.location, a.kind, b.kind);
  }


  friend inline
  bool operator==(const ADPlace& a, const ADPlace& b) {
    return (a.location==b.location) && (a.kind==b.kind);
  }


  friend inline
  std::ostream& operator<<(std::ostream& os, const ADPlace& cp) {
    return os << '(' << cp.location << ',' << cp.kind << ')';
  }
};


struct ADKey
{
  int kind;
  ADKey()=default;
  ADKey(int k) : kind(k) {}

  friend inline
  bool operator<(const ADKey& a, const ADKey& b) {
    return LazyLess(a.kind, b.kind);
  }

  friend inline
  bool operator==(const ADKey& a, const ADKey& b) {
    return a.kind==b.kind;
  }

  friend inline
  std::ostream& operator<<(std::ostream& os, const ADKey& cp) {
    return os << '(' << cp.kind << ')';
  }
};


// This is as much of the marking as the transition will see.
using Local=LocalMarking<Uncolored<SimpleToken>,Uncolored<HerdToken>>;
// Extra state to add to the system state. Will be passed to transitions.
struct WithParams {
  // Put our parameters here.
  std::map<ADParam,double> dparams;
  std::map<ADParam,int> iparams;
};


// The transition needs to know the local marking and any extra state.
using SIRTransition=ExplicitTransition<Local,RandGen,WithParams>;

using Dist=TransitionDistribution<RandGen>;
using ExpDist=ExponentialDistribution<RandGen>;

bool is_infectious(DiseaseState s) {
  return (s==DiseaseState::Subclinical) || (s==DiseaseState::Clinical) ||
    (s==DiseaseState::Latent);
}

bool is_infected(DiseaseState s) {
  return (s==DiseaseState::Latent) || is_infectious(s);
}

// Herd infects neighboring herd directly.
class Infect : public SIRTransition {
 private:
  double rate_;
 public:
  Infect(double rate) : rate_(rate) {}
  virtual ~Infect() {}

  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0, RandGen& rng) override {
    SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"Infect::Enabled begin");
    auto SC=lm.template GetToken<1>(0, [](const HerdToken& t)->bool {
      return t.disease_state==DiseaseState::Susceptible;
    });
    auto IC=lm.template GetToken<1>(1, [](const HerdToken& t)->bool {
      return is_infectious(t.disease_state);
    });
    if (IC.second && SC.second && SC.first && IC.first) {
      SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"Infect::Enabled enable "<<rate_);
      return {true, std::unique_ptr<ExpDist>(new ExpDist(rate_, te))};
    } else {
      SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"Infect::Enabled disable");
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) override {
    SMVLOG(BOOST_LOG_TRIVIAL(trace) << "Infect::Fire " << lm);
    const auto infector=[](HerdToken& t) {
      assert(t.disease_state==DiseaseState::Susceptible);
      t.disease_state=DiseaseState::Latent;
    };
    lm.template Move<1, 1, decltype(infector)>(0, 0, 1, infector);
  }
};




// The GSPN itself.
using SIRGSPN=
    ExplicitTransitions<ADPlace, ADKey, Local, RandGen, WithParams>;


SIRGSPN
BuildSystem(NAADSMScenario& scenario) {
  int64_t herd_cnt=scenario.herd_cnt();
  BuildGraph<SIRGSPN> bg;
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;

  int herd_kind=0;
  for (int location_idx=0; location_idx<herd_cnt; ++location_idx ) {
    bg.AddPlace({location_idx, herd_kind}, 1);
  }

  enum { infect, recover, birth, deaths, deathi, deathr };

  for (int source_idx=0; source_idx<herd_cnt; ++source_idx) {
    for (int target_idx=0; target_idx<herd_cnt; ++target_idx) {
      if (target_idx!=source_idx) {
        double rate=scenario.airborne_hazard(source_idx, target_idx);
        bg.AddTransition({infect},
          {Edge{{target_idx, herd_kind}, -1},
           Edge{{source_idx, herd_kind}, -1}},
          std::unique_ptr<SIRTransition>(new Infect(rate))
          );
      }
    }
  }

  // std::move the transitions because they contain unique_ptr.
  return std::move(bg.Build());
}



template<typename GSPN, typename SIRState>
struct SIROutput {
  GSPN& gspn_;
  std::shared_ptr<TrajectoryObserver> observer_;

  SIROutput(GSPN& gspn, std::shared_ptr<TrajectoryObserver> observer)
  : gspn_(gspn), observer_(observer) { 
  };

  int64_t step_cnt{0};

  void operator()(const SIRState& state) {
    auto transition=gspn_.VertexTransition(state.last_transition);
    SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"last transition "<<state.last_transition);

    // infection
    if (transition.kind==0) {
      auto transition_places=
          NeighborsOfTransition(gspn_, state.last_transition);
      int64_t target_internal_place=std::get<0>(transition_places[0]);
      int64_t source_internal_place=std::get<0>(transition_places[1]);
      int64_t target=gspn_.VertexPlace(
          std::get<0>(transition_places[0])).location;
      int64_t source=gspn_.VertexPlace(
          std::get<0>(transition_places[1])).location;

      size_t token_cnt=Length<1>(state.marking, target_internal_place);
      assert(token_cnt==1);

      auto res=Get<1>(state.marking, target_internal_place,
        [](const HerdToken& t)->bool {
          return is_infectious(t.disease_state);
        });
      assert(res.second==true);
      assert(res.first==true);
      auto sres=Get<1>(state.marking, source_internal_place,
        [](const HerdToken& t)->bool {
          return is_infectious(t.disease_state);
        });
      assert(sres.second==true);
      assert(sres.first==true);

      observer_->Step({transition.kind, target, source, state.CurrentTime()});
    }
    ++step_cnt;
  }

  void final(const SIRState& state) {
    BOOST_LOG_TRIVIAL(info) << "Took "<< step_cnt << " transitions.";
  }
};



int64_t SIR_run(double end_time, const std::vector<Parameter>& parameters,
    NAADSMScenario& scenario, std::shared_ptr<TrajectoryObserver> observer,
    RandGen& rng)
{
  int64_t herd_cnt=scenario.herd_cnt();
  int64_t herd_size=100;
  auto gspn=BuildSystem(scenario);

  // Marking of the net.
  static_assert(std::is_same<int64_t,SIRGSPN::PlaceKey>::value,
    "The GSPN's internal place type is int64_t.");
  using Mark=Marking<SIRGSPN::PlaceKey, Uncolored<SimpleToken>,
    Uncolored<HerdToken>>;
  using ADState=GSPNState<Mark,SIRGSPN::TransitionKey,WithParams>;

  ADState state;
  for (auto& cp : parameters) {
    state.user.dparams[cp.kind]=cp.value;
  }

  int64_t infected_herd=std::round<int64_t>(
      state.user.dparams[ADParam::FirstFarm]);
  for (int64_t init_idx=0; init_idx<herd_cnt; ++init_idx) {
    auto vertex_id=gspn.PlaceVertex({init_idx, 0});
    if (init_idx!=infected_herd) {
      Add<1>(state.marking, vertex_id, {DiseaseState::Susceptible, herd_size});
    } else {
      Add<1>(state.marking, vertex_id, {DiseaseState::Clinical, herd_size});
    }
  }

  //using Propagator=PropagateCompetingProcesses<int64_t,RandGen>;
  using Propagator=NonHomogeneousPoissonProcesses<int64_t,RandGen>;
  Propagator competing;
  using Dynamics=StochasticDynamics<SIRGSPN,ADState,RandGen>;
  Dynamics dynamics(gspn, {&competing});

  BOOST_LOG_TRIVIAL(debug) << state.marking;

  SIROutput<SIRGSPN,ADState> output_function(gspn, observer);

  dynamics.Initialize(&state, &rng);

  bool running=true;
  auto nothing=[](ADState&)->void {};
  double last_time=state.CurrentTime();
  while (running && state.CurrentTime()<end_time) {
    running=dynamics(state);
    if (running) {
      double new_time=state.CurrentTime();
      if (new_time-last_time<-1e-12) {
        BOOST_LOG_TRIVIAL(warning) << "last time "<<last_time <<" "
          << " new_time "<<new_time;
      }
      last_time=new_time;
      output_function(state);
    }
  }
  if (running) {
    BOOST_LOG_TRIVIAL(info)<<"Reached end time "<<state.CurrentTime();
  } else {
    BOOST_LOG_TRIVIAL(info)<<"No transitions left to fire at time "<<last_time;
  }
  output_function.final(state);
  return 0;
}

