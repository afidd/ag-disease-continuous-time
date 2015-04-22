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
  double hazard_per_day(const std::string& target, double distance);

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
};

class Herds {
  std::vector<Herd> state_;
 public:
  Herds()=default;
  void load(const std::string& filename);
};


void load_naadsm(const std::string& filename);


// _NAADSM_XML_H_
#endif
