#include<iostream>
#include<time.h>


//std::string time()

int main()
{
	std::cout << "hi" << std::endl; 
	std::cout << "clks/sec = " << CLOCKS_PER_SEC << std::endl;
	std::cout << "clks     = " << clock() << std::endl;
	clock_t start = clock();
	std::cout << "start     = " << start << std::endl;
	//for (int i=0;i < 100; ++i)
	//{
	//	std::cout << i << "\t" << clock() <<  std::endl;
		
	//}

	time_t t1 = time(NULL);
	for (int i=0;i < 100000; ++i) std::cout << i << "\t" << clock() <<  std::endl;
	//time_t t2 = time(&t1);
	time_t t2 = time(NULL);
	clock_t end = clock();

	std::cout << t1  << "\t" << t2 << std::endl;
	std::cout << "time = " << t2-t1  << " (sec) " << std::endl;
	std::cout << "clk = " << (end-start)/CLOCKS_PER_SEC  << " (sec) " << std::endl;


return 0;
}
