#pragma once 
#include <string>
#include "common.hpp"

namespace snakegen
{

struct MonsterDef
{
    std::string id;
    int power;
};

struct Monster
{
    std::string id;
    Coord pos;
    int amount;
};

class MonsterLib
{
public:
    MonsterLib();

    Monster genMonster(const Coord& pos, const int power);

private:
    std::vector<MonsterDef> mDefs;
};

}