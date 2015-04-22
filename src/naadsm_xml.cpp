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

void AirborneSpread::load_target(std::string target, pt::ptree& tree) {
  double p=tree.get<double>("prob-spread-1km");
  probability1km_[target]=p;
  assert(tree.get<int>("wind-direction-start.value")==0);
  assert(tree.get<int>("wind-direction-end.value")==360);
  assert(tree.get<int>("delay.probability-density-function.point")==0);
}


void load_naadsm(const std::string& filename) {
  // Load the file respecting UTF-8 locale.
  std::ifstream input_file_stream;
  input_file_stream.open(filename);
  if (!input_file_stream.is_open()) {
    BOOST_LOG_TRIVIAL(error)<<"Cannot open file "<<filename;
    return;
  }
  input_file_stream.imbue(std::locale("en_US.UTF-8"));

  std::vector<DiseaseModel> disease_model;
  std::map<std::string,AirborneSpread> airborne_spread;

  pt::ptree tree;
  pt::read_xml(input_file_stream, tree);
  auto models=tree.get_child("naadsm:disease-simulation.models");
  for (pt::ptree::value_type& v: models) {
    if (v.first.compare("disease-model")==0) {
      disease_model.push_back(DiseaseModel());
      disease_model[disease_model.size()-1].load_disease_model(v.second);
    } else if (v.first.compare("airborne-spread-exponential-model")==0) {
      auto attr=v.second.get_child("<xmlattr>");
      auto from=attr.get<std::string>("from-production-type");
      auto whom=attr.get<std::string>("to-production-type");
      if (airborne_spread.find(from)==airborne_spread.end()) {
        airborne_spread.emplace(from, from);
      }
      airborne_spread[from].load_target(whom, v.second);
    }
  }
}

