#define BOOST_TEST_MODULE a
#define BOOST_TEST_DYN_LINK
#include <numeric>
#include <filesystem>
#include <sstream>
#include <fstream>
#include <string>
#include <regex>
#include <algorithm>

#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>

#include "impl/generator.hpp"
#include "impl/rng.hpp"

class Test
{
public:
   Test()
   {

   }
};

std::string transformName(const std::string& src)
{
    std::string result;
    if(src.find('_') == std::string::npos)
    {
        result += src[0];
        for(auto iter = src.begin() + 1; iter != src.end(); ++iter)
        {
            if(isupper(*iter))
            {
                result += "_";
            }
            result += *iter;
        }   
    }
    else
    {
        result = src;
    }

    std::transform(result.begin(), result.end(), result.begin(), ::toupper);

    return result;
}

void parseCreature(const std::string& lines, std::map<std::string, int>& resultMap)
{
    std::string tag;
    int power = 0;
    {
        std::regex word_regex("/MapObjects/.*/(\\w*)\\..*");
        auto words_begin = std::sregex_iterator(lines.begin(), lines.end(), word_regex);
        auto words_end = std::sregex_iterator();
        for (std::sregex_iterator i = words_begin; i != words_end; ++i) 
        {
            std::smatch match = *i;
            std::string match_str = match[1];
            tag = transformName(match_str);
        }
    }
    {
        std::regex word_regex("\\<Power\\>(.*)\\</Power\\>");
        auto words_begin = std::sregex_iterator(lines.begin(), lines.end(), word_regex);
        auto words_end = std::sregex_iterator();
        for (std::sregex_iterator i = words_begin; i != words_end; ++i) 
        {
            std::smatch match = *i;
            std::string match_str = match[1];
            power = boost::lexical_cast<int>(match_str);
        }
    }

    resultMap[tag] = power;
}

BOOST_FIXTURE_TEST_CASE( TestSplit, Test ) 
{
    /*for(int i = 0; i < 10; ++i)
    {
        auto map = snakegen::genMap(snakegen::XL, Rng::gen32());

        BOOST_CHECK(map->zones.size() > 3);

        auto totalArea = std::accumulate(map->zones.begin(), map->zones.end(), 0, [](auto acc, auto x){
            return acc + x.area;
        });

        BOOST_CHECK_EQUAL(snakegen::XL * snakegen::XL, totalArea);

        std::cout << "-----" << std::endl;
        for(auto z : map->zones)
        {
            std::cout << z.area << " ";
        }
        std::cout << std::endl;
    }*/

    using namespace std;
    using namespace std::filesystem;

    std::map<std::string, int> resultMap; 

    for(const auto& entry : recursive_directory_iterator("Creatures")) 
    {
        // Is it a file / directory?
        bool isNormalFile = is_regular_file(entry);
        bool isDirectory = is_directory(entry);

        auto path = entry.path();

        if(isNormalFile)
        {
            std::ifstream infile(path);

            std::string lines;
            std::string line;
            while (std::getline(infile, line))
            {
                lines += line;
            }

            if(lines.find("</Creature>") != std::string::npos)
            {
                parseCreature(lines, resultMap);
            }
        }
    }

    std::ofstream outfile("monsters.def");
    for(auto& c : resultMap)
    {
        outfile << c.first << " " << c.second << std::endl;
    }
}