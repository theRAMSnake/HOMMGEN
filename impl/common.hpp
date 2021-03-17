#pragma once

#include <boost/polygon/voronoi.hpp>
#include <boost/polygon/polygon.hpp>
#include <boost/geometry.hpp>

using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;

typedef boost::geometry::model::d2::point_xy<int> Coord;
typedef boost::geometry::model::polygon<Coord> Polygon;

namespace boost 
{
namespace polygon 
{

template <>
struct geometry_concept<Coord> { typedef point_concept type; };

template <>
struct point_traits<Coord> 
{
    typedef int coordinate_type;

    static inline coordinate_type get(const Coord& point, orientation_2d orient) {
        return (orient == HORIZONTAL) ? point.x() : point.y();
    }
};

}
}