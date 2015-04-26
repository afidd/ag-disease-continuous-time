#include <set>
#include <exception>
#include <iostream>
#include <fstream>
#include <locale>
#include <fcntl.h>
#include "boost/property_tree/xml_parser.hpp"
#include "naadsm_xml.hpp"
#include "smv.hpp"
namespace pt = boost::property_tree;

std::string DiseaseModel::production_type() const { return name_; }

void DiseaseModel::load_disease_model(pt::ptree& tree) {
  auto attr=tree.get_child("<xmlattr>");
  name_=attr.get<std::string>("production-type");
  id_=attr.get<int>("production-type-id");
  std::cout << name_ << " " << id_ << std::endl;
  latent_=load_disease_pdf(tree.get_child("latent-period"));
  subclinical_=load_disease_pdf(tree.get_child("infectious-subclinical-period"));
  clinical_=load_disease_pdf(tree.get_child("infectious-clinical-period"));
  immune_=load_disease_pdf(tree.get_child("immunity-period"));
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


AirborneSpread::AirborneSpread(std::string source) : source_(source) {}

double AirborneSpread::hazard_per_day(
    const std::string& target, double dx) const {
  return std::exp(-probability1km_.at(target)*dx);
}

void AirborneSpread::load_target(std::string target, pt::ptree& tree) {
  double p=tree.get<double>("prob-spread-1km");
  probability1km_[target]=-std::log(p);
  assert(tree.get<int>("wind-direction-start.value")==0);
  assert(tree.get<int>("wind-direction-end.value")==360);
  assert(tree.get<int>("delay.probability-density-function.point")==0);
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
}

int64_t Herds::size() const { return state_.size(); }


void NAADSMScenario::load(const std::string& scenario,
    const std::string& herd) {
  load_scenario(scenario);
  herds_.load(herd);
}

int64_t NAADSMScenario::herd_cnt() const { return herds_.size(); }


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
}


std::vector<std::array<double,2>> NAADSMScenario::GetLocations() const {
  std::vector<std::array<double,2>> locations;
  for (const Herd& h : herds_.state_) {
    locations.push_back({h.latlong.first, h.latlong.second});
  }
  return locations;
}


double NAADSMScenario::airborne_hazard(int64_t source, int64_t target) const {
  const auto& source_prod=herds_.state_[source].production_type;
  const std::pair<double,double>& source_loc=herds_.state_[source].latlong;
  const auto& target_prod=herds_.state_[source].production_type;
  const std::pair<double,double>& target_loc=herds_.state_[source].latlong;

  double distance=::distancekm(source_loc.first, source_loc.second,
    target_loc.first, target_loc.second);
  return airborne_spread_.at(source_prod).hazard_per_day(
    target_prod, distance);
}


