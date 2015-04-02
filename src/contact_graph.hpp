#ifndef _CONTACT_GRAPH_HPP_
#define _CONTACT_GRAPH_HPP_ 1

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graphviz.hpp"
#include "boost/graph/erdos_renyi_generator.hpp"
#include "boost/graph/copy.hpp"
#include "boost/property_map/function_property_map.hpp"
#include "boost/dynamic_properties.hpp"
#include "smv.hpp"



//! Every contact graph will be of this type, a MutableGraph.
using ContactGraph=boost::adjacency_list<
                    boost::listS, // VertexList container
                    boost::setS, // OutEdgeList container
                    boost::undirectedS, // Directed/bidirectional
                    boost::no_property, // vert prop
                    boost::no_property, // Edge property
                    boost::no_property, // Graph property
                    boost::vecS // EdgeList container
                    >;


template<typename Random>
std::tuple<std::shared_ptr<ConstructGraph>,
  std::shared_ptr<boost::dynamic_properties>>
ErdosRenyiContactGraph(uint64_t node_cnt, double edge_fraction, Random& rn_gen) {
  // Use an edge list with a set in order to generate the graph without
  // having two edges connecting the same two nodes.
  using ERGen=boost::sorted_erdos_renyi_iterator<RNGen,ConstructType>;

  bool self_loops=false;
  auto graph=std::make_shared<ContactGraph>(ERGen(rn_gen, node_cnt,
    edge_fraction, self_loops), ERGen(), node_cnt);
  auto properties=std::make_shared<boost::dynamics_properties>();
  return std::make_tuple(graph, properties);
}


std::tuple<std::shared_ptr<ConstructGraph>,
  std::shared_ptr<boost::dynamic_properties>>
GraphMLContactGraph(std::string name) {
  auto g=std::make_shared<ContactGraph>();
  auto dynamic_properties=std::make_shared<boost::dynamic_properties>();
  std::ifstream input_file(name);
  boost::read_graphml(input_file, *g, *dynamic_properties);
  return std::make_tuple(graph, properties);
}



#endif
