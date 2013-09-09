#include "TRandom.h"
#include <iostream>


void DoRan(const int t)
{
	TRandom ran;

	for (int i=0;i<t;++i)
	{
		std::cout << ran.Uniform() << std::endl;
	}
}


void TRand() 
{
	for (int i=0;i <10; i++)
	{
		std::cout << "\n\n loop = " << i << std::endl;
		DoRan(2);
	}

}
