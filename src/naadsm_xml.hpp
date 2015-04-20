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
    std::string label = v.first;
    std::cout << label << std::endl;
  }
}


// _NAADSM_XML_H_
#endif
