#include<iostream>
#include<set>
#include<utility>


int main()
{

	std::set<int> myset;
	std::set<int>::iterator it;
	std::set<std::pair<int,int> > myset2;
	std::set<std::pair<int,int> >::iterator it2;

	myset.insert(50);
	myset.insert(10);
	myset.insert(20);
	myset.insert(40);
	myset.insert(20);

	for (it = myset.begin(); it!=myset.end(); it++)
	{
		std::cout << *it << std::endl;
	}

	std::cout << "\t" << std::endl;

	std::pair<unsigned, unsigned> pair1(1000,12);
	std::pair<unsigned, unsigned> pair2(2000,110);
	std::pair<unsigned, unsigned> pair3(1500,121);
	std::pair<unsigned, unsigned> pair4(100,50);
	myset2.insert(pair1);
	myset2.insert(pair2);
	myset2.insert(pair3);
	myset2.insert(pair4);

	for (it2 = myset2.begin(); it2!=myset2.end(); it2++)
	{
		 std::cout << it2->first << "\t" << it2->second <<std::endl;
	}
	
	return 0;
};
