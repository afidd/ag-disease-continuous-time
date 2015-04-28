#include <set>
#include <exception>
#include <iostream>
#include <fstream>
#include <locale>
#include <fcntl.h>
#include "boost/property_tree/xml_parser.hpp"
#include "boost/math/special_functions/fpclassify.hpp"
#include "naadsm_xml.hpp"
#include "smv.hpp"
namespace pt = boost::property_tree;

std::string DiseaseModel::production_type() const { return name_; }

void DiseaseModel::load_disease_model(pt::ptree& tree) {
  auto attr=tree.get_child("<xmlattr>");
  name_=attr.get<std::string>("production-type");
  id_=attr.get<int>("production-type-id");
  std::vector<std::string> periods={"latent-period",
    "infectious-subclinical-period", "infectious-clinical-period",
    "immunity-period"};
  for (const auto& period : periods) {
    transitions_.push_back(load_disease_pdf(tree.get_child(period)));
  }
}


/*! Turn what was read into a set of things to build.
 *  Keep a list of states as integers. 1 is susceptible,
 *  so maybe { 1, 2, 4, 5 } for susceptible, latent, clinical, immune.
 *  Keep a list of transitions, erasing those excluded in the description.
 *  0th is latent to clinical.
 *  1st is clinical to immune.
 *  2nd is immune back to susceptible.
 */
void DiseaseModel::build_states() {
  int current_state=1;
  states_.push_back(current_state); // susceptible=1
  int transition_idx=0;
  while (transition_idx<transitions_.size()) {
    ++current_state;
    if (has_state(transitions_[transition_idx])) {
      states_.push_back(current_state);
      ++transition_idx;
    } else {
      transitions_.erase(transitions_.begin()+transition_idx);
    }
  }
}


bool DiseaseModel::has_state(DistributionDescription& dist) {
  if (std::get<0>(dist)==DistributionEnum::point &&
    boost::math::fpclassify(std::get<1>(dist)["a"])==FP_ZERO) {
    return false;
  }
  return true;
}

DiseaseModel::DistributionDescription
DiseaseModel::load_disease_pdf(pt::ptree& tree) {
  auto pdf_tree=tree.get_child("probability-density-function");
  DistributionEnum dist_type{DistributionEnum::unknown};
  std::map<std::string,double> params;

  for (auto& ppair : pdf_tree) {
    if (ppair.first.compare("gamma")==0) {
      dist_type=DistributionEnum::gamma;
    } else if (ppair.first.compare("point")==0) {
      dist_type=DistributionEnum::point;
      // The parameter for this distribution is unnamed, stored as
      // a value.
      params["a"]=ppair.second.get_value<double>();
    } else if (ppair.first.compare("units")==0) {
      if (ppair.second.get<std::string>("xdf:unit").compare("day")!=0) {
        BOOST_LOG_TRIVIAL(warning)<<"One distribution doesn't use days. "
          "It uses "<<ppair.second.get<std::string>("xdf:unit") << " instead.";
      }
    } else if (ppair.first.compare("<xmlattr>")==0) {
      ;
    } else {
      BOOST_LOG_TRIVIAL(warning) << "Unknown distribution "<<ppair.first;
      dist_type=DistributionEnum::unknown;
    }
    if (dist_type!=DistributionEnum::unknown) {
      for (auto& param : ppair.second) {
        try {
          double value=param.second.get_value<double>();
          params[param.first]=value;
        } catch (pt::ptree_bad_data& e) {
          BOOST_LOG_TRIVIAL(error) << "Could not read parameter " << param.first
            << " for distribution " << ppair.first;
        }
      }
      return std::make_tuple(dist_type, params);
    } else {
      ; // continue looping until you find a distribution
    }
  }
  return std::make_tuple(dist_type, params);
}

std::ostream& operator<<(std::ostream& os, const DiseaseModel& dm) {
  os << "DiseaseModel("<<dm.name_<<") : "<<dm.id_ << std::endl;
  for (int s : dm.states_) {
    os << s << " ";
  }
  for (const auto& dd : dm.transitions_) {
    DistributionEnum de=std::get<0>(dd);
    const auto& params=std::get<1>(dd);
    if (de==DistributionEnum::gamma) {
      os << "gamma " << params.at("alpha") << " " << params.at("beta")
        << std::endl;
    } else if (de==DistributionEnum::point) {
      os << "point " << params.at("a") << std::endl;
    }
  }
  return os;
}

AirborneSpread::AirborneSpread(std::string source) : source_(source) {}

double AirborneSpread::hazard_per_day(
    const std::string& target, double dx) const {
  double probability=std::exp(-probability1km_.at(target)*dx);
  return -std::log(1-probability);
}

void AirborneSpread::load_target(std::string target, pt::ptree& tree) {
  double p=tree.get<double>("prob-spread-1km");
  probability1km_[target]=-std::log(p);
  assert(tree.get<int>("wind-direction-start.value")==0);
  assert(tree.get<int>("wind-direction-end.value")==360);
  assert(tree.get<int>("delay.probability-density-function.point")==0);
}

std::ostream& operator<<(std::ostream& os, const AirborneSpread& dm) {
  os << "AirborneSpread("<<dm.source_<<")"<<std::endl;
  for (const auto& entry : dm.probability1km_) {
    os << entry.first << " " << entry.second << std::endl;
  }
  return os;
}


/*! Distance computed on a spherical earth.
 *  Taken from http://williams.best.vwh.net/avform.htm.
 */
double distancekm(double lat1, double lon1, double lat2, double lon2) {
  constexpr double degrees_to_radians=boost::math::constants::pi<double>()/180;
  constexpr double radians_km=180*60*1.852/boost::math::constants::pi<double>();
  lat1=lat1*degrees_to_radians;
  lon1=lon1*degrees_to_radians;
  lat2=lat2*degrees_to_radians;
  lon2=lon2*degrees_to_radians;
  double d=2*std::asin(std::sqrt(std::pow(std::sin((lat1-lat2)/2), 2) + 
    std::cos(lat1)*std::cos(lat2)*std::pow(std::sin((lon1-lon2)/2), 2)));
  return d*radians_km;
}


double Herd::distancekm(const Herd& b) const {
  return ::distancekm(latlong.first, latlong.second,
    b.latlong.first, b.latlong.second);
}


void Herds::load(const std::string& filename) {
  std::ifstream input_file_stream;
  input_file_stream.open(filename);
  if (!input_file_stream.is_open()) {
    BOOST_LOG_TRIVIAL(error)<<"Cannot open file "<<filename;
    return;
  }
  input_file_stream.imbue(std::locale("en_US.UTF-8"));

  pt::ptree tree;
  pt::read_xml(input_file_stream, tree);
  auto& herds=tree.get_child("herds");
  for (pt::ptree::value_type& v: herds) {
    if (v.first.compare("herd")==0) {
      state_.emplace_back(
        v.second.get<int>("id"),
        v.second.get<std::string>("production-type"),
        v.second.get<int>("size"),
        std::make_pair(v.second.get<double>("location.latitude"),
                       v.second.get<double>("location.longitude")),
        v.second.get<std::string>("status")
        );
    } else if (v.first.compare("<xmlattr>")==0) {
      ; // also good
    } else {
      BOOST_LOG_TRIVIAL(warning) << "Unknown tag for herds " << v.first;
    }
  }
  // set index of ids
  for (int idx=0; idx<state_.size(); ++idx) {
    id_to_idx_[state_[idx].id]=idx;
  }
}

std::ostream& operator<<(std::ostream& os, const Herd& h) {
  os << "Herd("<<h.id<<", "<<h.production_type<<", "<<
    h.size<<", " << h.status << ", " << h.latlong.first << ", "
    << h.latlong.second << ")";
  return os;
}


int64_t Herds::size() const { return state_.size(); }

std::vector<int> Herds::herd_ids() const {
  std::vector<int> ids;
  ids.reserve(state_.size());
  for (auto& h : state_) {
    ids.push_back(h.id);
  }
  return ids;
}

std::ostream& operator<<(std::ostream& os, const Herds& h) {
  for (const auto& herd : h.state_) {
    os << herd << std::endl;
  }
  return os;
}

void NAADSMScenario::load(const std::string& scenario,
    const std::string& herd) {
  load_scenario(scenario);
  herds_.load(herd);
}

int64_t NAADSMScenario::herd_cnt() const { return herds_.size(); }

std::vector<int> NAADSMScenario::herd_ids() const { return herds_.herd_ids(); }

const std::vector<int>& NAADSMScenario::disease_states(int herd_id) const {
  int herd_idx=herds_.id_to_idx_.at(herd_id);
  const auto& prod_type=herds_.state_[herd_idx].production_type;
  return disease_model_.at(prod_type).states_;
}


const std::map<std::string,double>& NAADSMScenario::disease_transition(
  int herd_id, int transition_idx, DistributionEnum& transition_kind,
    int& start, int& finish) const {
  int herd_idx=herds_.id_to_idx_.at(herd_id);
  const auto& prod_type=herds_.state_[herd_idx].production_type;
  const DiseaseModel& dm=disease_model_.at(prod_type);
  transition_kind=std::get<0>(dm.transitions_[transition_idx]);
  start=dm.states_[transition_idx+1];
  int finish_idx=transition_idx+2;
  if (finish_idx>=dm.states_.size()) {
    finish_idx=0; // from immune to susceptible. It's a chain.
  }
  finish=dm.states_[finish_idx];
  return std::get<1>(dm.transitions_[transition_idx]);
}


int NAADSMScenario::disease_cnt(int herd_id) {
  int herd_idx=herds_.id_to_idx_.at(herd_id);
  const auto& prod_type=herds_.state_[herd_idx].production_type;
  return disease_model_.at(prod_type).transitions_.size();
}

void NAADSMScenario::load_scenario(const std::string& filename) {
  // Load the file respecting UTF-8 locale.
  std::ifstream input_file_stream;
  input_file_stream.open(filename);
  if (!input_file_stream.is_open()) {
    BOOST_LOG_TRIVIAL(error)<<"Cannot open file "<<filename;
    return;
  }
  input_file_stream.imbue(std::locale("en_US.UTF-8"));

  pt::ptree tree;
  pt::read_xml(input_file_stream, tree);
  auto models=tree.get_child("naadsm:disease-simulation.models");
  for (pt::ptree::value_type& v: models) {
    if (v.first.compare("disease-model")==0) {
      DiseaseModel flock_type;
      flock_type.load_disease_model(v.second);
      disease_model_[flock_type.production_type()]=flock_type;
    } else if (v.first.compare("airborne-spread-exponential-model")==0) {
      auto attr=v.second.get_child("<xmlattr>");
      auto from=attr.get<std::string>("from-production-type");
      auto whom=attr.get<std::string>("to-production-type");
      if (airborne_spread_.find(from)==airborne_spread_.end()) {
        airborne_spread_.emplace(from, from);
      }
      airborne_spread_[from].load_target(whom, v.second);
    }
  }
  for (auto& dm : disease_model_) {
    dm.second.build_states();
  }
}


std::vector<std::array<double,2>> NAADSMScenario::GetLocations() const {
  std::vector<std::array<double,2>> locations;
  for (const Herd& h : herds_.state_) {
    locations.push_back({h.latlong.first, h.latlong.second});
  }
  return locations;
}


double NAADSMScenario::airborne_hazard(int64_t source, int64_t target) const {
  source=herds_.id_to_idx_.at(source);
  target=herds_.id_to_idx_.at(target);
  const auto& source_prod=herds_.state_[source].production_type;
  const std::pair<double,double>& source_loc=herds_.state_[source].latlong;
  const auto& target_prod=herds_.state_[source].production_type;
  const std::pair<double,double>& target_loc=herds_.state_[source].latlong;

  double distance=::distancekm(source_loc.first, source_loc.second,
    target_loc.first, target_loc.second);
  return airborne_spread_.at(source_prod).hazard_per_day(
    target_prod, distance);
}

std::ostream& operator<<(std::ostream& os, const NAADSMScenario& s) {
  for (const auto& dm : s.disease_model_) {
    os << dm.first << " " << dm.second;
  }
  for (const auto& as : s.airborne_spread_) {
    os << as.first << " " << as.second;
  }
  os << s.herds_;
  return os;
}


