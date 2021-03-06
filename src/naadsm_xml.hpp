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
  std::string name_;
  int id_;
 public:
  std::vector<int> states_;
  std::vector<DistributionDescription> transitions_;

  DiseaseModel()=default;
  std::string production_type() const;

  void load_disease_model(boost::property_tree::ptree& tree);
  void build_states();
 private:
  friend std::ostream& operator<<(std::ostream& os, const DiseaseModel& dm);
  bool has_state(DistributionDescription& dist);
  DistributionDescription load_disease_pdf(boost::property_tree::ptree& tree);
};


class AirborneSpread {
  std::string source_;
  std::map<std::string,double> probability1km_;
 public:
  AirborneSpread()=default;
  AirborneSpread(std::string source);
  double spread_factor(const std::string& target) const;

  void load_target(std::string target, boost::property_tree::ptree& tree);
 private:
  friend std::ostream& operator<<(std::ostream& os, const AirborneSpread& dm);
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
  friend std::ostream& operator<<(std::ostream& os, const Herd& h);
};


class Herds {
 public:
  std::vector<Herd> state_;
  std::map<int,int> id_to_idx_;
  std::vector<double> infection_factor_;
  Herds()=default;
  void load(const std::string& filename);
  int64_t size() const;
  std::vector<int> herd_ids() const;
  void CalculateFactor();
 private:
  friend std::ostream& operator<<(std::ostream& os, const Herds& h);
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
  std::vector<int> herd_ids() const;
  const std::vector<int>& disease_states(int herd_id) const;
  const std::map<std::string,double>& disease_transition(
    int herd_id, int transition_idx, DistributionEnum& transition_kind,
    int& start, int& finish) const;
  int disease_cnt(int herd_id);
  std::vector<double> Distances();
 private:
  friend std::ostream& operator<<(std::ostream& os, const NAADSMScenario& s);
  void load_scenario(const std::string& filename);
};


// _NAADSM_XML_H_
#endif
