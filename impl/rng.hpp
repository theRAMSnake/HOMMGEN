#pragma once

class Rng
{
public:
   static void seed(const unsigned int seed);
   static bool genProbability(const double chance);
   static double genWeight();
   static double genPerturbation();
   static double genReal();
   static unsigned int genChoise(const unsigned int numOptions);
   static unsigned int gen32();
   static unsigned int genProbabilities(const double chance, const unsigned int number);
};