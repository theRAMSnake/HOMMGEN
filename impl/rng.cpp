#include "rng.hpp"
#include <boost/random.hpp>
#include <random>

std::random_device dv;

thread_local boost::random::mt19937 rng(dv());
thread_local boost::uniform_real<> doubleGen(-1.0, 1.0);
thread_local boost::uniform_real<> perturbationGen(-0.05, 0.05);
thread_local boost::uniform_real<> realGen(0, 1.0);

void Rng::seed(const unsigned int seed)
{
   rng.seed(seed);
}

bool Rng::genProbability(const double chance)
{
   return realGen(rng) < chance;
}

unsigned int Rng::genProbabilities(const double chance, const unsigned int number)
{
   unsigned int result = 0;
   for(unsigned int i = 0; i < number; ++i)
   {
      if(genProbability(chance))
      {
         result++;
      }
   }
   return result;
}

double Rng::genWeight()
{
   return doubleGen(rng);
}

double Rng::genReal()
{
   return realGen(rng);
}

double Rng::genPerturbation()
{
   return perturbationGen(rng);
}

unsigned int Rng::genChoise(const unsigned int numOptions)
{
   return rng() % numOptions;
}

unsigned int Rng::gen32()
{
    return rng();
}
