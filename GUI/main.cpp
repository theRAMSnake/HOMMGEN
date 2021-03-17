#include <nana/gui.hpp>
#include <nana/gui/place.hpp>
#include <nana/gui/screen.hpp>
#include <nana/gui/widgets/form.hpp>

#include "impl/generator.hpp"
#include "impl/rng.hpp"

#include <thread>
#include <iostream>
#include <mutex>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

typedef boost::geometry::model::d2::point_xy<int> point;
typedef boost::geometry::model::polygon<point> polygon;

const nana::color COLORS[] = {
   nana::colors::aqua	,
   nana::colors::beige	,
   nana::colors::blue	,
   nana::colors::brown	,
   nana::colors::coral ,
   nana::colors::gold ,
   nana::colors::gray ,
   nana::colors::green ,
   nana::colors::indigo	,
   nana::colors::ivory ,
   nana::colors::magenta ,
   nana::colors::maroon	,
   nana::colors::olive	,
   nana::colors::orange	,
   nana::colors::pink	,
   nana::colors::purple	,
   nana::colors::red		,
   nana::colors::salmon ,
   nana::colors::silver	,
   nana::colors::teal	,
   nana::colors::turquoise ,
   nana::colors::yellow	,
};

int main()
{
   nana::size sz = nana::screen().primary_monitor_size();
#ifdef WINDOWS
   nana::form fm(nana::rectangle{sz});
#else
   sz.width -= 50;
   nana::form fm(nana::rectangle{sz});
#endif

   fm.caption("Snake Map Generator");

   auto map = snakegen::genMap(214, Rng::gen32());

   nana::drawing draw(fm);
   draw.draw([&](nana::paint::graphics& graph){

      /*for(std::size_t i = 0; i < map->zones.size(); ++i)
      {
         auto& z = map->zones[i].polygon.outer();
         if(z.size() < 2)
         {
            continue;
         }

         if(map->zones[i].isStart)
         {
            for(int x = 0; x < 214; x++)
            {
               for(int y = 0; y < 214; y++)
               {
                  if(boost::geometry::intersects(map->zones[i].polygon, Coord{x, y}))
                  {
                     graph.set_pixel( 100 + x, 100 + y, nana::colors::blue);
                  }
               }
            }
         }
      }*/

      /*for(std::size_t i = 0; i < map->zones.size(); ++i)
      {
         auto& z = map->zones[i].polygon.outer();
         if(z.size() < 2)
         {
            continue;
         }
         for(std::size_t j = 0; j < z.size() - 1; ++j)
         {
            graph.line({100 + z[j].x(), 100 + z[j].y()}, {100 + z[j + 1].x(), 100 + z[j + 1].y()}, nana::colors::black);
         }
      }*/

      for(int x = 0; x < 214; x++)
      {
         for(int y = 0; y < 214; y++)
         {
            if(map->obstacleLayer(x, y) == snakegen::Presence::Occupied)
            {
               graph.set_pixel( 500 + x, 100 + y, nana::colors::black);
            }
         }
      }

      /*for(auto& c : map->connections)
      {
         Coord c1;
         Coord c2;
         boost::geometry::centroid(map->zones[c.first()].polygon, c1);
         boost::geometry::centroid(map->zones[c.second()].polygon, c2);
         graph.line({500 + c1.x(), 100 + c1.y()}, {500 + c2.x(), 100 + c2.y()}, nana::colors::red);
      }*/

      for(auto& o : map->objects)
      {
         nana::rectangle rect{500 + o.pos.x(), 100 + o.pos.y(), o.size.x(), o.size.y()};
         graph.rectangle(rect, true, nana::colors::blue);
      }

      for(auto& o : map->monsters)
      {
         graph.set_pixel(500 + o.pos.x(), 100 + o.pos.y(), nana::colors::red);
         graph.set_pixel(500 + o.pos.x() + 1, 100 + o.pos.y(), nana::colors::red);
         graph.set_pixel(500 + o.pos.x() + 1, 100 + o.pos.y() + 1, nana::colors::red);
         graph.set_pixel(500 + o.pos.x(), 100 + o.pos.y() + 1, nana::colors::red);
      }

      /*for(int x = 0; x < 214; x++)
      {
         for(int y = 0; y < 214; y++)
         {
            if(map->roadLayer(x, y) == snakegen::Presence::Occupied)
            {
               graph.set_pixel(500 + x, 100 + y, nana::colors::red);
            }
         }
      }*/

      for(int x = 0; x < 214; x++)
      {
         for(int y = 0; y < 214; y++)
         {
            auto d = map->distanceLayer(x, y);
            auto xo = 1000 + x;
            auto yo = 100 + y;
            auto color = nana::colors::black;

            if(d > 0 && d < 75)
            {
               color = nana::colors::green;
            }
            else if(d > 74 && d < 150)
            {
               color = nana::colors::yellow;
            }
            if(d > 149 && d < 225)
            {
               color = nana::colors::red;
            }
            if(d > 224 && d < 2000)
            {
               color = nana::colors::purple;
            }

            graph.set_pixel( xo, yo, color);
         }
      }
   });
   
   fm.show();

   nana::exec();
}