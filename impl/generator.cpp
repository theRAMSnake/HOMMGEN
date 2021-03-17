#include "generator.hpp"
#include "rng.hpp"
#include <limits>
#include <cmath>
#include <iostream>
#include "monsterlib.hpp"

namespace std
{

bool operator<(const Coord& a, const Coord& b)
{
    return std::make_pair(a.x(), a.y()) < std::make_pair(b.x(), b.y());
}

}

namespace snakegen
{

static ObjectLib gObjectLibrary;
static MonsterLib gMonsterLibrary;
const int BORDER_LENGTH = 7;

std::vector<Polygon> voronoiToPolygonList(
    const voronoi_diagram<double>& vd, const std::vector<Coord>& points, const int size)
{
    std::vector<Polygon> result;

    for (auto it = vd.cells().begin(); it != vd.cells().end(); ++it) 
    {
        const voronoi_diagram<double>::cell_type &cell = *it;
        const voronoi_diagram<double>::edge_type *edge = cell.incident_edge();
        
        Polygon p;

        // This is convenient way to iterate edges around Voronoi cell.
        do 
        {
            //do
            if (edge->is_finite())
            {
                p.outer().push_back({static_cast<int>(edge->vertex0()->x()), static_cast<int>(edge->vertex0()->y())});
                p.outer().push_back({static_cast<int>(edge->vertex1()->x()), static_cast<int>(edge->vertex1()->y())});
            }
            else
            {
                auto p1 = points[it->source_index()];
                auto p2 = points[edge->twin()->cell()->source_index()];
                Coord origin;
                Coord direction;

                origin.x((p1.x() + p2.x()) * 0.5);
                origin.y((p1.y() + p2.y()) * 0.5);

                direction.x(p1.y() - p2.y());
                direction.y(p2.x() - p1.x());

                auto koef = size / std::max(fabs(direction.x()), fabs(direction.y()));

                if (edge->vertex0() == NULL)
                {
                    p.outer().push_back({origin.x() - static_cast<int>(direction.x() * koef), 
                        origin.y() - static_cast<int>(direction.y() * koef)});
                }
                else
                {
                    p.outer().push_back({static_cast<int>(edge->vertex0()->x()), static_cast<int>(edge->vertex0()->y())});
                }

                if (edge->vertex1() == NULL)
                {
                    p.outer().push_back({origin.x() + static_cast<int>(direction.x() * koef), 
                        origin.y() + static_cast<int>(direction.y() * koef)});
                }
                else
                {
                    p.outer().push_back({static_cast<int>(edge->vertex1()->x()), static_cast<int>(edge->vertex1()->y())});
                }
            }

            edge = edge->next();
        } while (edge != cell.incident_edge());

        if(boost::geometry::equals(p.outer().back(), p.outer().front()))
        {
            p.outer().push_back(p.outer().front());
        }

        boost::geometry::correct(p);
        result.push_back(p);
    }

    Polygon mapClipArea{{{0, 0}, {0, size}, {size, size}, {size, 0}, {0, 0}}};

    //Clip zones
    for(auto& p : result)
    {
        std::vector<Polygon> out;
        boost::geometry::intersection(p, mapClipArea, out);

        if(out.size() != 0)
        {
            p = out[0];
        }
        else
        {
            continue;
        }
    }

    return result;
}

Coord getShiftedPoint(const Coord& pt1, const Coord& pt2, const int distance)
{
    Coord middlePoint{(pt1.x() + pt2.x()) / 2, (pt1.y() + pt2.y()) / 2};

    //1. Find out linear equation
    auto A = pt1.y() - pt2.y();
    auto B = pt2.x() - pt1.x();
    //auto C = pt1.x * pt2.y - pt2.x * pt1.y;

    auto newM = static_cast<double>(B) / (static_cast<double>(A) + 0.01);

    //2. Get parameter function
    auto mag = sqrt(1 + newM * newM);
    std::pair<double, double> N {1.0 / mag, newM / mag};

    return {middlePoint.x() + static_cast<int>(distance * N.first), middlePoint.y() + static_cast<int>(distance * N.second)};
}

struct PairedEdge
{
    std::size_t pIndex;
    std::size_t iIndex;
};

bool fuzzyCompare(const Coord& pt1, const Coord& pt2)
{
    return std::abs(pt1.x() - pt2.x()) < 2 && std::abs(pt1.y() - pt2.y()) < 2;
}

std::optional<PairedEdge> findPairedEdge(const std::vector<Polygon>& zones, const Coord& pt1, const Coord& pt2)
{
    for(std::size_t pIndex = 0; pIndex < zones.size(); ++pIndex)
    {
        if(zones[pIndex].outer().size() != 0)
        {
            for(std::size_t i = 0; i < zones[pIndex].outer().size() - 1; ++i)
            {
                if(fuzzyCompare(pt1, zones[pIndex].outer()[i + 1]) && fuzzyCompare(pt2, zones[pIndex].outer()[i]))
                {
                    PairedEdge result{pIndex, i};
                    return result;
                }
            }
        }
    }

    return {};
}

Coord getRandomShiftedPoint(const Coord& pt1, const Coord& pt2)
{
    auto shiftPercentage = static_cast<int>(Rng::genChoise(50)) - 25;
    auto len = sqrt(pow(pt2.x() - pt1.x(), 2) + pow(pt2.y() - pt1.y(), 2));
    auto d = shiftPercentage / 100.0 * len;

    return getShiftedPoint(pt1, pt2, d);
}

std::vector<Polygon> morphEdges(const std::vector<Polygon>& zones)
{
    std::vector<Polygon> result = zones;

    struct MorphOperation
    {
        Coord pt1;
        Coord pt2;

        Coord insertPoint;
    };

    std::vector<MorphOperation> operations; 

    for(std::size_t pIndex = 0; pIndex < zones.size(); ++pIndex)
    {
        if(zones[pIndex].outer().size() != 0)
        {
            for(std::size_t i = 0; i < zones[pIndex].outer().size() - 1; ++i)
            {
                auto pt1 = zones[pIndex].outer()[i];
                auto pt2 = zones[pIndex].outer()[i + 1];
                auto pairedEdge = findPairedEdge(zones, pt1, pt2);
                if(pairedEdge)
                {
                    Coord shifted = getRandomShiftedPoint(pt1, pt2);

                    if(std::find_if(operations.begin(), operations.end(), [&](auto x)
                    {
                        return fuzzyCompare(x.pt1, pt2) && fuzzyCompare(x.pt2, pt1);
                    }) == operations.end())
                    {
                        operations.push_back({pt1, pt2, shifted});
                    }
                }
            }
        }
    }

    for(auto& op : operations)
    {
        for(std::size_t pIndex = 0; pIndex < result.size(); ++pIndex)
        {
            if(result[pIndex].outer().size() != 0)
            {   
                for(std::size_t i = 0; i < result[pIndex].outer().size() - 1; ++i)
                {
                    if((fuzzyCompare(op.pt1, result[pIndex].outer()[i]) && fuzzyCompare(op.pt2, result[pIndex].outer()[i + 1])) ||
                        (fuzzyCompare(op.pt2, result[pIndex].outer()[i]) && fuzzyCompare(op.pt1, result[pIndex].outer()[i + 1])))
                    {
                        result[pIndex].outer().insert(result[pIndex].outer().begin() + i + 1, op.insertPoint);
                        break;
                    }
                }
            }
        }
    }

    return result;
}

void allocateZones(Map& map, const unsigned int size)
{
    auto numAdditionalZones = 2 + 3 + Rng::genChoise(20);

    std::vector<Coord> points;
    
    for(unsigned int i = 0; i < numAdditionalZones; ++i)
    {
        points.push_back({static_cast<int>(Rng::genChoise(size)), static_cast<int>(Rng::genChoise(size))});//Maybe worth to distribute points
    }

    voronoi_diagram<double> vd;
    construct_voronoi(points.begin(), points.end(), &vd);

    auto zones = voronoiToPolygonList(vd, points, size);

    //Morph edges (twice)
    zones = morphEdges(zones);
    zones = morphEdges(zones);
    zones = morphEdges(zones);

    for(auto& z : zones)
    {
        map.zones.push_back({z});
    }

    for(std::size_t zIndex = 0; zIndex < map.zones.size(); ++zIndex)
    {
        auto& z = map.zones[zIndex];
        for(std::size_t i = 0; i < z.polygon.outer().size() - 1; ++i)
        {
            if((fuzzyCompare(z.polygon.outer()[i], z.polygon.outer()[i + 1])))
            {
                z.polygon.outer().erase(z.polygon.outer().begin() + i + 1);
                --i;
            }
        }
    }
}

double calculateArea(const Polygon& p)
{
    return boost::geometry::area(p);
}

void selectSpawns(Map& map, const unsigned int size)
{
    auto spawnScore = [size](double area)
    {
        const double respIdealSize = static_cast<double>(size * size) * 0.1;
        if(area < respIdealSize)
        {
            return 99999.0;
        }
        else
        {
            return std::abs(area - respIdealSize);
        }
    };

    std::vector<std::pair<std::size_t, double>> indexToAreaMap;
    for(std::size_t i = 0; i < map.zones.size(); ++i)
    {
        indexToAreaMap.push_back({i, calculateArea(map.zones[i].polygon)});
    }

    std::sort(indexToAreaMap.begin(), indexToAreaMap.end(), [=](auto x, auto y)
    {
        return spawnScore(x.second) < spawnScore(y.second);
    });

    map.zones[indexToAreaMap[0].first].isStart = true;
    map.zones[indexToAreaMap[1].first].isStart = true;
}

void renderDot(const int X, const int Y, Matrix<Presence>& target, const int size, const Presence brush)
{
    auto n = size / 2 - 1;
    
    int minX = std::max(0, X - n);
    int maxX = std::min(target.width(), X + n);
    int minY = std::max(0, Y - n);
    int maxY = std::min(target.heigth(), Y + n);

    for(int i = minX; i <= maxX; ++i)
    {
        for(int j = minY; j <= maxY; ++j)
        {
            target(i, j) = brush;
        }
    }
}

void renderLine(const Coord& p1, const Coord& p2, Matrix<Presence>& target, const int size, const Presence brush)
{
    auto dx = p2.x() - p1.x();
    auto dy = p2.y() - p1.y();
    auto len = std::max(std::abs(dx), std::abs(dy));

    double deltaX = static_cast<double>(dx) / len;
    double deltaY = static_cast<double>(dy) / len;
    double X = p1.x();
    double Y = p1.y();

    renderDot(X, Y, target, size, brush);

    auto i = 1;
    while(i <= len)
    {
        X += deltaX;
        Y += deltaY;

        renderDot(X, Y, target, size, brush);

        i++;
    }
}

void renderRect(const Coord& topLeft, const Coord& sz, Matrix<Presence>& target, const Presence brush)
{
    for(int x = topLeft.x(); x < topLeft.x() + sz.x(); ++x)
    {
        for(int y = topLeft.y(); y < topLeft.y() + sz.y(); ++y)
        {
            target(x, y) = brush;
        }
    }
}

std::optional<std::size_t> findZoneThroughBorder(const std::vector<Zone>& zones, const Coord& pt1, const Coord& pt2)
{
    for(std::size_t zIndex = 0; zIndex < zones.size(); ++zIndex)
    {
        auto& z = zones[zIndex];
        for(std::size_t i = 0; i < z.polygon.outer().size() - 1; ++i)
        {
            if((fuzzyCompare(pt2, z.polygon.outer()[i]) && fuzzyCompare(pt1, z.polygon.outer()[i + 1])))
            {
                return zIndex;
            }
        }
    }
    
    return {};
}

void paintBorders(Map& map)
{
    for(auto& z : map.zones)
    {
        for(std::size_t i = 0; i < z.polygon.outer().size() - 1; ++i)
        {
            int lineWidth = BORDER_LENGTH;

            if(z.isStart)
            {
                auto otherZ = findZoneThroughBorder(map.zones, z.polygon.outer()[i], z.polygon.outer()[i + 1]);
                if(otherZ && map.zones[(*otherZ)].isStart)
                {
                    lineWidth = 15;
                }
            }
            renderLine(z.polygon.outer()[i], z.polygon.outer()[i + 1], map.obstacleLayer, lineWidth, Presence::Occupied);
        }
    }
}

std::vector<Coord> getSharedBorder(const std::vector<Zone>& zones, const std::size_t index1, const std::size_t index2)
{
    std::vector<Coord> result;

    for(std::size_t i = 0; i < zones[index1].polygon.outer().size() - 1; ++i)
    {
        auto otherZ = findZoneThroughBorder(zones, zones[index1].polygon.outer()[i], zones[index1].polygon.outer()[i + 1]);
        if(otherZ && (*otherZ) == index2)
        {
            if(result.empty())
            {
                result.push_back(zones[index1].polygon.outer()[i]);
                result.push_back(zones[index1].polygon.outer()[i + 1]);
            }
            else
            {
                result.push_back(zones[index1].polygon.outer()[i + 1]);
            }
        }
    }

    return result;
}

void connectZones(Map& map)
{
    std::map<Connection, int> totalLengthPerBorder;

    //Connect zones if shared border is at least X long
    for(std::size_t zIndex = 0; zIndex < map.zones.size(); ++zIndex)
    {
        auto& z = map.zones[zIndex];
        for(std::size_t i = 0; i < z.polygon.outer().size() - 1; ++i)
        {
            auto otherZ = findZoneThroughBorder(map.zones, z.polygon.outer()[i], z.polygon.outer()[i + 1]);
            if(otherZ)
            {
                if(zIndex == *otherZ)
                {
                    //Can happen due to fuzzy compare
                    continue;
                }

                if(z.isStart && map.zones[*otherZ].isStart)
                {
                    continue;
                }

                totalLengthPerBorder[{zIndex, *otherZ}] += boost::geometry::distance(z.polygon.outer()[i], z.polygon.outer()[i + 1]);
            }
        }
    }

    for(auto& i : totalLengthPerBorder)
    {
        if(i.second > 39)
        {
            map.connections.insert(i.first);
        }
    }

    for(auto& c : map.connections)
    {
        auto sharedBorder = getSharedBorder(map.zones, c.first(), c.second());
        for(int i = 0; i < 20; ++i) //20 attempts
        {
            auto randomSegment = Rng::genChoise(sharedBorder.size() - 1);
            auto seg = std::make_pair(sharedBorder[randomSegment], sharedBorder[randomSegment + 1]); 

            Coord shifted1 = getShiftedPoint(seg.first, seg.second, BORDER_LENGTH);
            Coord shifted2 = getShiftedPoint(seg.first, seg.second, -BORDER_LENGTH);

            if(map.obstacleLayer(shifted1.x(), shifted1.y()) == Presence::Free &&
                map.obstacleLayer(shifted2.x(), shifted2.y()) == Presence::Free)
            {
                renderLine(shifted1, shifted2, map.obstacleLayer, 3, Presence::Free);

                Coord middlePoint{(shifted1.x() + shifted2.x()) / 2, (shifted1.y() + shifted2.y()) / 2};
                map.zones[c.first()].exitPoints.push_back(middlePoint);
                map.zones[c.second()].exitPoints.push_back(middlePoint);

                break;
            }
        }
    }

    //Verify that each zone has at least one connection
    for(std::size_t i = 0; i < map.zones.size(); ++i)
    {
        bool cFound = false;

        for(auto& c : map.connections)
        {
            if(c.first() == i || c.second() == i)
            {
                cFound = true;
            }
        }

        if(!cFound)
        {
            throw std::runtime_error("Invalid map generated: disconnected zone");
        }
    }
}

Coord genRandomPoint(const Polygon& p)
{
    boost::geometry::model::box<Coord> box;
    boost::geometry::envelope(p, box);

    auto minPt = box.min_corner();
    auto maxPt = box.max_corner();

    Coord ranPt{-1, -1};
    while(!boost::geometry::intersects(p, ranPt))
    {
        ranPt.x(Rng::genChoise(maxPt.x() - minPt.x()) + minPt.x());
        ranPt.y(Rng::genChoise(maxPt.y() - minPt.y()) + minPt.y());
    }

    return ranPt;
}

bool isFreeArea(const Matrix<Presence>& m, const Coord& topLeft, const int area)
{
    for(int x = topLeft.x(); x < topLeft.x() + area; ++x)
    {
        for(int y = topLeft.y(); y < topLeft.y() + area; ++y)
        {
            if(!m.contains(x, y) || m(x, y) == Presence::Occupied)
            {
                return false;
            }
        }
    }

    return true;
}

void generateStrategicObstacles(Map& map)
{
    auto minZoneSizeToHaveObstacles = 1000;
    for(auto& z : map.zones)
    {
        if(boost::geometry::area(z.polygon) > minZoneSizeToHaveObstacles)
        {
            const int numAttempts = boost::geometry::area(z.polygon) / 100;//For now
            for(int i = 0; i < numAttempts; ++i)
            {
                auto pt1 = genRandomPoint(z.polygon);
                auto pt2 = genRandomPoint(z.polygon);

                if(isFreeArea(map.obstacleLayer, pt1, 13) && isFreeArea(map.obstacleLayer, pt2, 13))
                {
                    std::vector<Coord> pts {pt1, pt2};
                    
                    for(int k = 0; k < 7; ++k)
                    {
                        std::vector<Coord> newPts;

                        for(std::size_t ptPos = 0; ptPos < pts.size() - 1; ptPos++)
                        {
                            newPts.push_back(pts[ptPos]);

                            Coord shifted = getRandomShiftedPoint(pts[ptPos], pts[ptPos + 1]);

                            if(isFreeArea(map.obstacleLayer, shifted, 10) && boost::geometry::intersects(z.polygon, shifted))
                            {
                                newPts.push_back(shifted);
                            }
                        }

                        newPts.push_back(pts.back());

                        pts = newPts;
                    }

                    for(std::size_t ptPos = 0; ptPos < pts.size() - 1; ptPos++)
                    {
                        renderLine(pts[ptPos], pts[ptPos + 1], map.obstacleLayer, 10, Presence::Occupied);
                    }
                    
                }
            }
        }
    }
}

Coord offset(const Coord pt, const int off)
{
    return {pt.x() + off, pt.y() + off};
}

void placeObject(Map& map, const Coord& pt, Object&& o)
{
    const int AURA_SIZE = 3;
    o.pos = pt;
    
    map.objects.push_back(o);

    auto topLeft = offset(pt, -AURA_SIZE);
    auto sz = gObjectLibrary.getDef(o.id).size;
    renderRect(topLeft, {sz.x() + AURA_SIZE * 2, sz.y() + AURA_SIZE * 2}, map.availabilityLayer, Presence::Occupied);
}

bool isFarFromOthers(const std::vector<Object>& objects, const Coord& pt)
{
    for(auto& o : objects)
    {
        if(boost::geometry::distance(o.pos, pt) < 50)
        {
            return false;
        }
    }

    return true;
}

void placeTowns(Map& map, const std::size_t size)
{
    //Place 2 starting towns and other 7 in biggest zones
    for(auto& z : map.zones)
    {
        if(z.isStart)
        {
            bool townPlaced = false;
            const int numAttempts = 25;
            for(int i = 0; i < numAttempts; ++i)
            {
                auto pt = genRandomPoint(z.polygon);
                if(isFreeArea(map.availabilityLayer, pt, 13))
                {
                    placeObject(map, pt, gObjectLibrary.create("town"));
                    townPlaced = true;
                    break;
                }
            }

            if(!townPlaced)
            {
                throw std::runtime_error("cannot place starting town");
            }
        }
    }

    //If zones < numTowns - place multiple towns in biggest zones
    auto numTowns = size / 30;

    std::vector<std::pair<std::size_t, double>> indexToAreaMap;
    for(std::size_t i = 0; i < map.zones.size(); ++i)
    {
        indexToAreaMap.push_back({i, calculateArea(map.zones[i].polygon)});
    }

    std::sort(indexToAreaMap.begin(), indexToAreaMap.end(), [=](auto x, auto y)
    {
        return x.second > y.second;
    });

    int townsLeft = numTowns;
    int attempts = 1000;
    std::size_t i = 0;
    while(townsLeft > 0 && attempts != 0)
    {
        auto& z = map.zones[indexToAreaMap[i].first];
        if(!z.isStart)
        {
            auto pt = genRandomPoint(z.polygon);
            if(isFreeArea(map.availabilityLayer, pt, 13) && isFarFromOthers(map.objects, pt))
            {
                placeObject(map, pt, gObjectLibrary.create("town"));
                townsLeft--;
            }
        }

        attempts--;

        i++;
        if(i == indexToAreaMap.size())
        {
            i = 0;
        }
    }
}

Coord getExitPoint(const Object& o)
{
    auto sz = gObjectLibrary.getDef(o.id).size;

    switch(o.rotation)
    {
        case Rotation::CW180:
            return {o.pos.x() + sz.x() / 2, o.pos.y()};
        ;
        case Rotation::No:
            return {o.pos.x() + sz.x() / 2, o.pos.y() + sz.y()};
        ;
        case Rotation::CW270:
            return {o.pos.x() + sz.x(), o.pos.y() + sz.y() / 2};
        ;
        case Rotation::CW90:
            return {o.pos.x(), o.pos.y() + sz.y() / 2};
        ;
    }

    throw std::runtime_error("Invalid rotation");
}

std::vector<Coord> getNeighbors(const Coord pt)
{
    return {
        {pt.x() - 1, pt.y() - 1},
        {pt.x(), pt.y() - 1},
        {pt.x() + 1, pt.y() - 1},
        {pt.x() - 1, pt.y()},
        {pt.x() + 1, pt.y()},
        {pt.x() - 1, pt.y() + 1},
        {pt.x(), pt.y() + 1},
        {pt.x() + 1, pt.y() + 1}
    };
}

double diff(const Coord pt1, const Coord pt2)
{
    if(pt1.x() != pt2.x() && pt1.y() != pt2.y())
    {
        return 1.41;
    }
    else
    {
        return 1.0;
    }
}

void buildDistanceMap(
    Matrix<double>& distances, 
    const std::set<Coord>& initialSet, 
    std::function<bool(const Coord)> filter,
    std::function<void(const Coord, const Coord)> onSmallestNeighbor
    )
{
    std::set<Coord> ptsToConsider = initialSet;

    for(auto& c : initialSet)
    {
        distances(c.x(), c.y()) = 0;
    }
    
    while(!ptsToConsider.empty())
    {
        std::set<Coord> newPtsToConsider;

        for(auto& pt : ptsToConsider)
        {
            auto nbrs = getNeighbors(pt);
            for(auto& nbr : nbrs)
            {
                auto newDistance = distances(pt.x(), pt.y()) + diff(pt, nbr);
                if(newDistance < distances(nbr.x(), nbr.y()) && filter(nbr))
                {
                    distances(nbr.x(), nbr.y()) = newDistance;
                    newPtsToConsider.insert(nbr);
                    onSmallestNeighbor(nbr, pt);
                }
            }
        }

        ptsToConsider = newPtsToConsider;
    }
}

void buildDistanceMap(Map& map)
{
    map.distanceLayer.setAll(9999);

    auto startPt1 = getExitPoint(map.objects[0]/*tmp*/);
    auto startPt2 = getExitPoint(map.objects[1]/*tmp*/);

    std::set<Coord> ptsToConsider;
    ptsToConsider.insert(startPt1);
    ptsToConsider.insert(startPt2);

    buildDistanceMap(map.distanceLayer, ptsToConsider, [&map](auto x){
        return map.obstacleLayer(x.x(), x.y()) == Presence::Free;
    },
    [](auto x, auto y){}
    );
}

bool getFreeSpot(Map& map, const Coord size, Coord& spot)
{
    auto sz = std::max(size.x(), size.y()) + 2;

    for(int x = 0; x < map.availabilityLayer.width() - sz; ++x)
    {
        for(int y = 0; y < map.availabilityLayer.heigth() - sz; ++y)
        {
            if(isFreeArea(map.availabilityLayer, {x, y}, sz))
            {
                spot = {x + 1, y + 1};
                return true;
            }
        }
    }
    return false;
}

void populateObjects(Map& map, const std::size_t size)
{
    Coord spot;
    auto maxSize = gObjectLibrary.getMaxObjectSize();
    auto seekSize = std::max(maxSize.x(), maxSize.y());

    while(seekSize != 0)
    {
        while(getFreeSpot(map, {seekSize, seekSize}, spot))
        {
            placeObject(map, spot, gObjectLibrary.generate({seekSize, seekSize}, map.distanceLayer(spot.x(), spot.y())));
        }

        seekSize--;
    }
}

void initAvailability(Matrix<Presence>& availabilityLayer, const Matrix<Presence>& obstacleLayer)
{
    for(int x = 0; x < obstacleLayer.width(); ++x)
    {
        for(int y = 0; y < obstacleLayer.heigth(); ++y)
        {
            if(obstacleLayer(x, y) == Presence::Occupied)
            {
                availabilityLayer(x, y) = Presence::Occupied;
            }
        }
    }
}

void updateBlockedAreas(Map& map)
{
    //All areas with distance {init} - should be obstacles and unavailable
    for(int x = 0; x < map.distanceLayer.width(); ++x)
    {
        for(int y = 0; y < map.distanceLayer.heigth(); ++y)
        {
            if(map.distanceLayer(x, y) > 9998)
            {
                map.availabilityLayer(x, y) = Presence::Occupied;
                map.obstacleLayer(x, y) = Presence::Occupied;
            }
        }
    }
}

bool isSameZone(const std::vector<Zone>& zones, const Coord& pt1, const Coord& pt2)
{
    for(auto& z : zones)
    {
        if(boost::geometry::intersects(z.polygon, pt1) &&
            boost::geometry::intersects(z.polygon, pt2))
            {
                return true;
            }
    }

    return false;
}

Coord clip(const Coord& a, const int upperBound)
{
    return {
        std::max(0, std::min(a.x(), upperBound)),
        std::max(0, std::min(a.y(), upperBound))
    };
}

void placeEssentials(Map& map, const Coord& pt)
{
    //2 Basic mines and a gold mine
    std::vector<Object> objectsToPlace{
        gObjectLibrary.create("sawmill"),
        gObjectLibrary.create("ore_pit"),
        gObjectLibrary.create("gold_mine"),
    };

    int attempts = 1000;
    for(auto& o : objectsToPlace)
    {
        while(attempts != 0)
        {
            attempts++;
            Coord ranpt{pt.x() + static_cast<int>(Rng::genChoise(50)) - 25, pt.y() + static_cast<int>(Rng::genChoise(50)) - 25};
            ranpt = clip(ranpt, map.size);

            if(map.distanceLayer(ranpt.x(), ranpt.y()) < 50 && 
                isFreeArea(map.availabilityLayer, ranpt, gObjectLibrary.getDef(o.id).size.y() + 2) &&
                isSameZone(map.zones, pt, ranpt))
            {
                placeObject(map, ranpt, std::move(o));
                break;
            }
        }
    }
}

std::vector<Object> findCastles(const Zone& z, const Map& m)
{
    std::vector<Object> result;

    for(auto& o : m.objects)
    {
        if(o.id == "town" && boost::geometry::intersects(z.polygon, o.pos))
        {
            result.push_back(o);
        }
    }

    return result;
}

std::vector<Coord> createIntersectionPoints(const Zone& z, const int amount, const Map& m)
{
    std::vector<Coord> result;

    for(int i = 0; i < amount; ++i)
    {
        int attempts = 1000;
        while(attempts != 0)
        {
            auto pt = genRandomPoint(z.polygon);
            if(isFreeArea(m.obstacleLayer, pt, 1))
            {
                result.push_back(pt);
                break;
            }

            attempts--;
        }
    }

    return result;
}

Coord findApproximatePoint(const Matrix<double>& distances, const Coord& destination)
{
    int radius = 1;
    while(radius < 10)
    {
        std::map<Coord, double> distancesOfRadius;

        for(int x = destination.x() - radius; x < destination.x() + radius; ++x)
        {
            for(int y = destination.y() - radius; y < destination.y() + radius; ++y)
            {
                if(distances.contains(x, y))
                {
                    distancesOfRadius[{x, y}] = distances(x, y);
                }
            }
        }

        Coord result = destination;
        double minimal = 9000;
        for(auto& c : distancesOfRadius)
        {
            if(c.second < minimal)
            {
                minimal = c.second;
                result = c.first;
            }
        }

        if(minimal != 9000)
        {
            return result;
        }

        radius++;
    }

    throw std::runtime_error("cannot find even approximate path");
}

std::vector<Coord> shortestPath(const Map& map, const Coord& startPt, const Coord& endPt)
{
    std::map<Coord, Coord> prevs;
    Matrix<double> distances(std::size_t(map.size), std::size_t(map.size));
    distances.setAll(9999);

    std::set<Coord> ptsToConsider;
    ptsToConsider.insert(startPt);

    buildDistanceMap(distances, ptsToConsider, [&map](auto x){
        return map.obstacleLayer(x.x(), x.y()) == Presence::Free;
    },
    [&prevs](auto x, auto y){
        prevs[x] = y;
    }
    );

    auto destination = endPt;
    if(distances(endPt.x(), endPt.y()) > 9000)
    {
        destination = findApproximatePoint(distances, endPt);
    }

    std::vector<Coord> result;

    auto unwindPt = destination;
    while(!boost::geometry::equals(unwindPt, startPt))
    {
        result.push_back(unwindPt);
        unwindPt = prevs[unwindPt];
    }

    result.push_back(startPt);

    return result;
}

void renderRoad(Matrix<Presence>& roads, const Map& map, const Coord& startPt, const Coord& endPt)
{
    if(!boost::geometry::equals(startPt, endPt))
    {
        auto path = shortestPath(map, startPt, endPt);
        for(auto& pt: path)
        {
            roads(pt.x(), pt.y()) = Presence::Occupied;
        }
    }
}

Coord findClosestRoadPoint(const Matrix<Presence>& roads, const Map& map, const Coord& pt)
{
    Matrix<double> distances(std::size_t(map.size), std::size_t(map.size));
    distances.setAll(9999);

    std::set<Coord> ptsToConsider;
    ptsToConsider.insert(pt);

    buildDistanceMap(distances, ptsToConsider, [&map](auto x){
        return map.obstacleLayer(x.x(), x.y()) == Presence::Free;
    },
    [](auto x, auto y){
    }
    );

    Coord result = pt;
    double curDist = 9999;
    for(int x = 0; x < distances.width(); ++x)
    {
        for(int y = 0; y < distances.heigth(); ++y)
        {
            if(distances(x, y) < curDist && roads(x, y) == Presence::Occupied && pt.x() != x && pt.y() != y)
            {
                curDist = distances(x, y);
                result.x(x);
                result.y(y);
            }
        }
    }

    return result;
}

void drawRoads(Map& map)
{
    //For each zone
    for(auto& z : map.zones)
    {
        Matrix<Presence> localRoads(map.size, map.size);

        //  Find exit points
        auto exitPts = z.exitPoints;

        //  Find castles
        auto castles = findCastles(z, map);
        for(auto& c : castles)
        {
            exitPts.push_back(getExitPoint(c));
        }

        //  Create intersection points
        auto intersectionPoints = createIntersectionPoints(z, std::max(std::size_t(2), exitPts.size() - 1), map);
        if(intersectionPoints.size() < 2)
        {
            std::cout << "intersection points less than 2";
        }
        renderRoad(localRoads, map, intersectionPoints[0], intersectionPoints[1]);

        for(std::size_t i = 2; i < intersectionPoints.size(); ++i)
        {
            auto pt = findClosestRoadPoint(localRoads, map, intersectionPoints[i]);
            renderRoad(localRoads, map, pt, intersectionPoints[i]);
        }
        for(std::size_t i = 0; i < exitPts.size(); ++i)
        {
            auto pt = findClosestRoadPoint(localRoads, map, exitPts[i]);
            renderRoad(localRoads, map, pt, exitPts[i]);
        }

        map.roadLayer.merge(localRoads);
    }
}

void placeGrail(Map& map)
{
    Coord spot;
    
    getFreeSpot(map, {1, 1}, spot);
    placeObject(map, spot, gObjectLibrary.create("grail"));
}

void generateSmallObstacles(Map& map, const Matrix<Presence>& spots)
{
    const double CHANCE = 1.0 / 20;
    for(int x = 0; x < map.size; ++x)
    {
        for(int y = 0; y < map.size; ++y)
        {
            if(spots(x, y) == Presence::Free && Rng::genProbability(CHANCE))
            {
                map.obstacleLayer(x, y) = Presence::Occupied;
            }
        }
    }
}

void generateFreeResources(Map& map, const Matrix<Presence>& spots)
{
    const double CHANCE = 1.0 / 25;
    for(int x = 0; x < map.size; ++x)
    {
        for(int y = 0; y < map.size; ++y)
        {
            if(spots(x, y) == Presence::Free && Rng::genProbability(CHANCE))
            {
                placeObject(map, {x, y}, gObjectLibrary.create("resource"));
            }
        }
    }
}

Matrix<Presence> compileFreeSpotsMap(const Map& map)
{
    //obstacleLayer + objects + guardians
    Matrix<Presence> result = map.obstacleLayer;
    
    for(auto &o : map.objects)
    {
        renderRect(o.pos, o.size, result, Presence::Occupied);
    }

    //guards here

    return result;
}

Side rotate(const Side& original, const Rotation& rotation)
{
    if(original == Side::Middle || rotation == Rotation::No)
    {
        return original;
    }

    if(rotation == Rotation::CW90)
    {
        switch (original)
        {
        case Side::Top:
            return Side::Right;

        case Side::Left:
            return Side::Top;

        case Side::Right:
            return Side::Bottom;

        case Side::Bottom:
            return Side::Left;

        case Side::TopLeft:
            return Side::TopRight;

        case Side::TopRight:
            return Side::BottomRight;

        case Side::BottomLeft:
            return Side::TopLeft;

        case Side::BottomRight:
            return Side::BottomLeft;
        }
    }
    else if(rotation == Rotation::CW180)
    {
        switch (original)
        {
        case Side::Top:
            return Side::Bottom;

        case Side::Left:
            return Side::Right;

        case Side::Right:
            return Side::Left;

        case Side::Bottom:
            return Side::Top;

        case Side::TopLeft:
            return Side::BottomRight;

        case Side::TopRight:
            return Side::BottomLeft;

        case Side::BottomLeft:
            return Side::TopRight;

        case Side::BottomRight:
            return Side::TopLeft;
        }
    }
    else /*rotation == Rotation::CW270*/
    {
        switch (original)
        {
        case Side::Top:
            return Side::Left;

        case Side::Left:
            return Side::Right;

        case Side::Right:
            return Side::Top;

        case Side::Bottom:
            return Side::Right;

        case Side::TopLeft:
            return Side::BottomLeft;

        case Side::TopRight:
            return Side::TopLeft;

        case Side::BottomLeft:
            return Side::BottomRight;

        case Side::BottomRight:
            return Side::TopRight;
        }
    }
}

Coord getGuardPos(const Object& o)
{
    auto middleX = o.pos.x() + o.size.x() / 2;
    auto middleY = o.pos.y() + o.size.y() / 2;

    auto& def = gObjectLibrary.getDef(o.id);

    switch (rotate(def.guardPos, o.rotation))
    {
    case Side::Top:
        return {middleX, o.pos.y() - 1};

    case Side::Left:
        return {o.pos.x() - 1, middleY};

    case Side::Right:
        return {o.pos.x() + o.size.x(), middleY};

    case Side::Bottom:
        return {middleX, o.pos.y() + o.size.y()};

    case Side::TopLeft:
        return {o.pos.x() - 1, o.pos.y() - 1};

    case Side::TopRight:
        return {o.pos.x() + o.size.x(), o.pos.y() - 1};

    case Side::BottomLeft:
        return {o.pos.x() - 1, o.pos.y() + o.size.y()};

    case Side::BottomRight:
        return {o.pos.x() + o.size.x(), o.pos.y() + o.size.y()};

    case Side::Middle:
        return {middleX, middleY};
    }

    throw std::runtime_error("Invalid guard pos");
}

void setGuardians(Map& map)
{
    //Find all zone guards pos
    std::map<Coord, Connection> exitPointToConnectionMap;

    for(std::size_t i = 0; i < map.zones.size(); ++i)
    {
        auto& z = map.zones[i];
        auto exitPts = z.exitPoints;
        for(auto& e : exitPts)
        {
            if(exitPointToConnectionMap.contains(e))
            {   
                continue;
            }

            //Find paired zone
            for(std::size_t j = 0; j < map.zones.size(); ++j)
            {
                if(i == j)
                {
                    continue;
                }

                for(auto& e2 : map.zones[j].exitPoints)
                {
                    if(e2.x() == e.x() && e2.y() == e.y())
                    {
                        exitPointToConnectionMap.insert(std::make_pair(e, Connection(i, j)));
                        break;
                    }
                }
            }
        }
    }

    assert(exitPointToConnectionMap.size() == map.connections.size());

    //Determine total zone value
    std::vector<int> totalValuePerZone;
    for(auto& z : map.zones)
    {
        int totalValue = 0;

        for(int x = 0; x < map.size; ++x)
        {
            for(int y = 0; y < map.size; ++y)
            {
                if(boost::geometry::intersects(z.polygon, Coord{x, y}) && map.obstacleLayer(x, y) == Presence::Free)
                {
                    totalValue += map.distanceLayer(x, y);
                }
            }
        }

        totalValuePerZone.push_back(totalValue);
    }

    // for(std::size_t i = 0; i < map.zones.size(); ++i)
    // {
    //     std::cout << "Zone: " << i << " value: " << totalValuePerZone[i] << std::endl;
    // }

    //Place zones guards
    for(auto& x : exitPointToConnectionMap)
    {
        auto& z1 = map.zones[x.second.first()];
        auto& z2 = map.zones[x.second.second()];

        auto avgZoneValue = totalValuePerZone[x.second.first()] + totalValuePerZone[x.second.second()];
        auto modifier = (z1.isStart || z2.isStart) ? 0.1 : 0.25;

        map.monsters.push_back(gMonsterLibrary.genMonster(x.first, static_cast<int>(avgZoneValue * modifier)));
    }
    
    //Place objects guards
    for(std::size_t i = 2; i < map.objects.size(); ++i)
    {
        auto& o = map.objects[i];

        if(o.value != 0)
        {
            auto pos = getGuardPos(o);
            map.monsters.push_back(gMonsterLibrary.genMonster(pos, o.value));
        }
    }
}

std::unique_ptr<Map> genMap(const std::size_t size, const std::size_t seed)
{
    Rng::seed(seed);

    auto result = std::make_unique<Map>(size);

    std::cout << "Start" << std::endl;
    allocateZones(*result, size); 
    selectSpawns(*result, size);
    paintBorders(*result);
    connectZones(*result);
    generateStrategicObstacles(*result); 
    initAvailability(result->availabilityLayer, result->obstacleLayer);
    placeTowns(*result, size);
    buildDistanceMap(*result);
    updateBlockedAreas(*result);

    placeEssentials(*result, getExitPoint(result->objects[0]/*tmp*/));
    placeEssentials(*result, getExitPoint(result->objects[1]/*tmp*/));

    placeGrail(*result);

    std::cout << "Populating Objects" << std::endl;
    populateObjects(*result, size); 

    std::cout << "Drawing Roads" << std::endl;
    drawRoads(*result);

    std::cout << "Setting Guardians" << std::endl;
    setGuardians(*result);

    //auto freeSpotsMap = compileFreeSpotsMap(*result);
    //generateSmallObstacles(*result, freeSpotsMap);

    //freeSpotsMap = compileFreeSpotsMap(*result);
    //generateFreeResources(*result, freeSpotsMap);

    //setBiomes(*result);

    std::cout << "Done" << std::endl;

    return result;
}

}