#include <iostream>


// 0.0 <= x < 0.2 =  0
// 0.2 <= x < 0.4 =  1
// 0.4 <= x < 0.6 =  2
// 0.6 <= x < 0.8 =  3
// 0.8 <= x < 1.0 =  4
// 1.0 <= x < 1.2 =  5 
// 1.2 <= x < 1.4 =  6
// 1.4 <= x < 1.6 =  7
// 1.6 <= x < 1.8 =  8
// 1.8 <= x < 2.0 =  9
// 2.0 <= x < 2.2 = 10
// 2.2 <= x < 2.4 = 11
// 2.4 <= x < 2.6 = 12
// 2.6 <= x < 2.8 = 13
// 2.8 >=         = 14

void EtaBin(const double fEta)
{
	int scale=10;
	double fEtaBinSize = 0.2;
	//if (fEta >1.0) scale = 10;
	
	float fbin = (fEta/fEtaBinSize);  
	std::cout << "fbin = " << fbin << " Int(fbin) = "<< (int)fbin << std::endl;
	int bin = (fEta/fEtaBinSize);  
	//int bin = -1;
	if (fEta>2.8) bin = 14;
	else 
	{
	//	bin = tbin/2;   
	}

	
	std::cout << "fEta["<<fEta << "] Etabin[" << bin << "] fbin[" << fbin << "]"  << std::endl;
}


void RunEtaBin()
{
	for (int i=1; i< 30; ++i)
	{
		EtaBin(i*0.1324);
	}
}
