#ifndef _NAADSM_XML_H_
#define _NAADSM_XML_H_ 1

#include <string>
#include <set>
#include <exception>
#include <iostream>
#include <fstream>
#include <locale>
#include <fcntl.h>
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/xml_parser.hpp"
#include "boost/foreach.hpp"
namespace pt = boost::property_tree;

enum class DistributionEnum : int { unknown, gamma, point };


std::tuple<DistributionEnum,std::map<std::string,double>> 
load_disease_pdf(pt::ptree& tree) {
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


void load_disease_model(pt::ptree& tree) {
  auto latent_dist=load_disease_pdf(tree.get_child("latent-period"));
  auto sc_dist=load_disease_pdf(tree.get_child("infectious-subclinical-period"));
  auto cl_dist=load_disease_pdf(tree.get_child("infectious-clinical-period"));
  auto im_dist=load_disease_pdf(tree.get_child("immunity-period"));
}

void load_airborne_spread(pt::ptree& tree) {

}


void load_xml(const std::string& filename) {
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
  for (pt::ptree::value_type& v: tree.get_child("naadsm:disease-simulation.models")) {
    if (v.first.compare("disease-model")==0) {
      load_disease_model(v.second);
    } else if (v.first.compare("airborne-spread-exponential-model")==0) {
      load_airborne_spread(v.second);
    }
  }
}


// _NAADSM_XML_H_
#endif
