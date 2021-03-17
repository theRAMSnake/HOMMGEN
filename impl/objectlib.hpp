#pragma once 
#include <string>
#include "common.hpp"

namespace snakegen
{

enum class Side
{
    Top,
    Left,
    Bottom,
    Right,
    TopLeft,
    TopRight,
    BottomLeft,
    BottomRight,
    Middle
};

struct ObjectDef
{
    std::string id;
    std::string func;
    Coord size{0, 0};
    double modifier = 1.0;
    int value = 0;
    Side guardPos;
};

enum class Rotation
{
    No,
    CW90,
    CW180,
    CW270
};

struct Object
{
    std::string id;
    Rotation rotation;
    Coord pos; //topleft
    Coord size;
    double value;
};

class ObjectLib
{
public:
    ObjectLib(); //+Town

    const ObjectDef& getDef(const std::string& id) const;
    Object create(const std::string& id);
    Object generate(const Coord& sz, const double distance);

    Coord getMaxObjectSize() const;

private:
    std::map<std::string, std::vector<ObjectDef>> mDefs;
};

}