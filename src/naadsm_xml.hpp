#ifndef _NAADSM_XML_H_
#define _NAADSM_XML_H_ 1

#include <string>
#include <set>
#include <exception>
#include <map>
#include <tuple>
#include <string>
#include "boost/property_tree/ptree.hpp"

enum class DistributionEnum : int { unknown, gamma, point };


class DiseaseModel {
 public:
  using DistributionDescription=
    std::tuple<DistributionEnum,std::map<std::string,double>>;
 private:
  DistributionDescription latent_, subclinical_, clinical_, immune_;
  std::string name_;
  int id_;
 public:
  DiseaseModel()=default;
  std::string production_type() const;

  void load_disease_model(boost::property_tree::ptree& tree);
 private:
  DistributionDescription load_disease_pdf(boost::property_tree::ptree& tree);
};


class AirborneSpread {
  std::string source_;
  std::map<std::string,double> probability1km_;
 public:
  AirborneSpread()=default;
  AirborneSpread(std::string source);
  double hazard_per_day(const std::string& target, double distance) const;

  void load_target(std::string target, boost::property_tree::ptree& tree);
};

struct Herd {
  int id;
  std::string production_type;
  int size;
  std::pair<double,double> latlong;
  std::string status;
  Herd()=default;
  Herd(int id, std::string prod, int s, std::pair<double,double> l,
    std::string status) : id(id), production_type(prod), size(s),
    latlong(l), status(status) {}
  double distancekm(const Herd& b) const;
};

class Herds {
 public:
  std::vector<Herd> state_;
  Herds()=default;
  void load(const std::string& filename);
  int64_t size() const;
};


class NAADSMScenario {
  std::map<std::string,DiseaseModel> disease_model_;
  std::map<std::string,AirborneSpread> airborne_spread_;
  Herds herds_;
 public:
  void load(const std::string& scenario, const std::string& herd);
  int64_t herd_cnt() const;
  std::vector<std::array<double,2>> GetLocations() const;
  double airborne_hazard(int64_t source, int64_t target) const;
 private:
  void load_scenario(const std::string& filename);
};


// _NAADSM_XML_H_
#endif
