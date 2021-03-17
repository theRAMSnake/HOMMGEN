#include "objectlib.hpp"
#include "rng.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace snakegen
{

Side parseSide(const std::string& src)
{
    if(src == "top") return Side::Top;
    if(src == "left") return Side::Left;
    if(src == "right") return Side::Right;
    if(src == "bottom") return Side::Bottom;
    if(src == "top_left") return Side::TopLeft;
    if(src == "top_right") return Side::TopRight;
    if(src == "bottom_left") return Side::BottomLeft;
    if(src == "bottom_right") return Side::BottomRight;
    if(src == "middle") return Side::Middle;

    throw std::runtime_error("Cannot parse side");
}

ObjectLib::ObjectLib()
{
    namespace pt = boost::property_tree;
    pt::ptree root;
    pt::read_json("objects.json", root);

    auto cats = root.get_child("categories");
    for(auto c : cats)
    {
        std::vector<ObjectDef> items;

        for(auto i : c.second.get_child("items"))
        {
            items.push_back({
                i.second.get<std::string>("id"),
                i.second.get<std::string>("func"),
                {i.second.get<int>("sizex"), i.second.get<int>("sizey")}
                });

            if(i.second.get_optional<double>("modifier"))
            {
                items.back().modifier = i.second.get<double>("modifier");
            }
            if(i.second.get_optional<int>("value"))
            {
                items.back().value = i.second.get<int>("value");
            }
            if(i.second.get_optional<std::string>("guard"))
            {
                items.back().guardPos = parseSide(i.second.get<std::string>("guard"));
            }
        }

        mDefs[c.second.get<std::string>("category")] = items;
    }

    /*for(auto c : mDefs)
    {
        std::cout << c.first << ":" << std::endl;
        for(auto i : c.second)
        {
            std::cout << "\t" << i.id << std::endl;
        }
    }*/

    mDefs["ungeneratable"].push_back({"town", "", {11, 11}, 0.0, 75000, Side::Bottom});
    mDefs["ungeneratable"].push_back({"grail", "", {1, 1}});
}

const ObjectDef& ObjectLib::getDef(const std::string& id) const
{
    for(auto& c : mDefs)
    {
        for(auto& i : c.second)
        {
            if(i.id == id)
            {
                return i;
            }
        }
    }

    throw std::runtime_error("Object not found");
}

Object ObjectLib::create(const std::string& id)
{   
    return {id, static_cast<Rotation>(Rng::genChoise(4)), {0, 0}, getDef(id).size, getDef(id).value};
}

double probFunc(const std::string& funcName, const double distance)
{
    if(funcName == "green_strict")
    {
        return distance < 75 ? 1.0 : 0.0;
    }
    else if(funcName == "identity")
    {
        return 1.0;
    }
    else if(funcName == "yellow_smooth")
    {
        if(distance > 74 && distance < 150)
        {
            return 1.0;
        }
        else
        {
            return std::max(0.05, 1.0 - std::abs(112 - distance) * 0.01);
        }
    }
    else if(funcName == "green_smooth")
    {
        if(distance < 75)
        {
            return 1.0;
        }
        else
        {
            return std::max(0.05, 1.0 - std::abs(75 - distance) * 0.02);
        }
    }
    else if(funcName == "linear")
    {
        return distance * 0.005;
    }
    else if(funcName == "linear++")
    {
        return distance * 0.0025;
    }
    else if(funcName == "red++")
    {
        return distance > 225 ? 1.0 : 0.0;
    }
    else
    {
        throw std::runtime_error("Unknown distance func");
    }
}

double calcValue(const double distance, const ObjectDef& def)
{
    if(def.id == "deposit" || def.id == "elemental_prison")
    {
        return 500 * distance;
    }

    return def.value;
}

Object ObjectLib::generate(const Coord& sz, const double distance)
{
    while(true)
    {
        auto iter = mDefs.end();
        while(iter == mDefs.end() || iter->first == "ungeneratable")
        {
            iter = std::next(mDefs.begin(), Rng::genChoise(mDefs.size()));
        }
        
        std::vector<std::pair<std::size_t, double>> index2chance;
        for(std::size_t i = 0; i < iter->second.size(); ++i)
        {
            auto& o = iter->second[i];
            if(o.size.x() <= sz.x() && o.size.y() <= sz.y())
            {
                index2chance.push_back({i, Rng::genReal() * probFunc(o.func, distance) * o.modifier});
            }
        }

        if(index2chance.empty())
        {
            continue;
        }

        std::sort(index2chance.begin(), index2chance.end(), [](auto x, auto y){return x.second < y.second;});

        auto& o = iter->second[index2chance.back().first];
        Object result{o.id, static_cast<Rotation>(Rng::genChoise(4)), {}, o.size, calcValue(distance, o)};

        return result;
    }
}

Coord ObjectLib::getMaxObjectSize() const
{
    Coord sz{0, 0};

    for(auto& c : mDefs)
    {
        if(c.first == "ungeneratable")
        {
            continue;
        }

        for(auto i : c.second)
        {
            if(sz.x() < i.size.x())
            {
                sz = i.size;
            }
        }
    }

    return sz;
}

}