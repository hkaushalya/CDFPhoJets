#include <iostream>
#include <vector>


int main(const int loopto=100000000, const int mtd=2)
{
	assert (loopto>=1 && "loopto must be >=1");
	
	
	if (mtd == 1)
	{
		std::cout << "method 1" << std::endl;
		int val[loopto];
		for (int i=0; i < loopto; ++i) val[i] = i;
	} else if (mtd == 2)
	{
		std::cout << "method 2" << std::endl;
		std::vector<int> vec;
		for (int i=0; i < loopto; ++i) vec.push_back(i);

	}
return 0;
};
