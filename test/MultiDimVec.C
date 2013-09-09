#include<iostream>
#include<vector>

using namespace std;

//study muti dimention vector array
int main()
{

	vector< vector<int> > vec;

	for (int i=0; i<10; ++i)
	{
		vector<int> row;		//create am empty row
		for (int j=0; j<20; ++j)
		{
			row.push_back(i*j);
		}
		vec.push_back(row);
	}

	std::cout << "size = " << vec.size() << "\t" << vec.at(0).size() << std::endl;
	
	for (int i=0; i<10; ++i)
	{
		for (int j=0; j<20; ++j)
		{
			//std::cout << vec[i][j] << "  ";
			std::cout << vec.at(i).at(j) << "  ";
		}
		std::cout << " " << std::endl;
	}

	vec.clear();
	std::cout << "size = " << vec.size() << std::endl;
	return 0;
}
