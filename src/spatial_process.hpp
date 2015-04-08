#ifndef _SPATIAL_PROCESS_HPP_
#define _SPATIAL_PROCESS_HPP_ 1

#include <vector>
#include <array>
#include <random>
#include <cmath>
#include "boost/math/constants/constants.hpp"


/*! Generates a list of points with no covariance in a rectangle.
 * 
 *  \param bounds rectangular bounds, (left, right, down, up)
 *  \param cnt The number of points to generate
 *  \param rng A random number generator
 */
template<typename RNG>
std::vector<std::array<double,2>> complete_spatial_randomness(
        std::vector<double> bounds, int64_t cnt, RNG& rng) {
    assert(cnt>0);
    assert(bounds[1]>bounds[0]);
    assert(bounds[3]>bounds[2]);

    auto point=std::vector<std::array<double,2>>(cnt, {0, 0});
    std::uniform_real_distribution<double> left_right(bounds[0], bounds[1]);
    std::uniform_real_distribution<double> up_down(bounds[2], bounds[3]);
    auto back_ptr=point.begin();
    for (int64_t i=0; i<cnt; ++i) {
      *back_ptr++={left_right(rng), up_down(rng)};
    }
    return point;
}



/*! A hard sphere process, so no two points are too close to each other.
 *  This places points and keeps guessing until it finds a spot that
 *  isn't close to another point.
 *  If each point covered a disc on the surface, the areal fraction
 *  is the area of each point's disc multiplied by the number of points,
 *  divided by the total area. This is meant to give a rough feeling for
 *  how much points are pushed apart, even when you change the size of
 *  the domain and number of points.
 *
 *  \param bounds rectangular bounds of domain (left, right, up, down)
 *  \param cnt The number of points to create.
 *  \param areal_fraction Safe to keep below 0.5. They get packed in there.
 *  \param rng random number generator
 */
template<typename RNG>
std::vector<std::array<double,2>> hard_sphere_process(
    std::vector<double> bounds, int64_t cnt, double areal_fraction, RNG& rng) {

    assert(cnt>0);
    assert(bounds[1]>bounds[0]);
    assert(bounds[3]>bounds[2]);

    double min_distance2=4*areal_fraction*(bounds[1]-bounds[0])*
        (bounds[3]-bounds[2])/(boost::math::constants::pi<double>()*cnt);
    std::uniform_real_distribution<double> left_right(bounds[0], bounds[1]);
    std::uniform_real_distribution<double> up_down(bounds[2], bounds[3]);

    auto point=std::vector<std::array<double,2>>(cnt, {0, 0});
    int64_t draw_idx=0;
    int64_t found_idx=0;
    while (found_idx<cnt && draw_idx<cnt*10000) {
        auto xy=std::array<double,2>{left_right(rng), up_down(rng)};
        bool found{true};
        for (int64_t search_idx=0; search_idx<found_idx; ++search_idx) {
            auto x=point[search_idx][0];
            auto y=point[search_idx][1];
            double distance2=pow(x-xy[0], 2) + pow(y-xy[1], 2);
            if (distance2<min_distance2) {
                found=false;
                break;
            }
        }
        if (found) {
            point[found_idx]=xy;
            ++found_idx;
        }
        ++draw_idx;
    }
  return point;
}



#endif //_SPATIAL_PROCESS_HPP_
