#include "monsterlib.hpp"
#include "rng.hpp"
#include <fstream>

namespace snakegen
{

MonsterLib::MonsterLib()
{
    std::ifstream infile("monsters.def");
    std::string id;
    int power = 0;

    while(infile)
    {
        infile >> id;
        infile >> power;

        mDefs.push_back(MonsterDef{id, power});
    }
}

Monster MonsterLib::genMonster(const Coord& pos, const int power)
{
    assert(power < 4500 * 250);

    while(true)
    {
        auto& ran = mDefs[Rng::genChoise(mDefs.size())];
        auto count = power / ran.power;
        if(count < 251)
        {
            return {ran.id, pos, count};
        }
    }

    throw std::runtime_error("Cannot generate guard");
}

}