#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <istream>
#include <vector>
#include <set>

/*
I am writing this to search for duplciate strings in a given file.
All entries are assumed to be in one column.

*/


int find_duplicate_files (std::string file)
{

	if (file.length()==0) {
		std::cout << "pl. specify a file."<< std::endl;
		return 1;
	}
	
	ifstream in(file.c_str());
	if (! in) {
		std::cout << "File not found. pl. check."<<  std::endl;
		return 1;
	}

	std::string temp;
	std::vector<std::string> wordlist;

	while (in >> temp) wordlist.push_back(temp);
	std::cout << "Words found = " << wordlist.size() << std::endl;
	
	std::set<std::string> uniq, duplicates;

	for (std::vector<std::string>::iterator it = wordlist.begin(); it != wordlist.end(); it++) {
		//std::cout << "Searching for  " << *it << std::endl;
		bool found = false;
		std::vector<std::string>::iterator it2;
		it2 = it;
		it2++;
		for (; it2 != wordlist.end(); it2++) {
			if (*it2 == *it) {
				duplicates.insert(*it);
				//std::cout << "Duplicate files = " << *it  << " and " << *it2<< std::endl;
				found = true;
				break;
			}
		}
		if (! found) {
			uniq.insert(*it);
			//std::cout << "Unique file = " << *it << std::endl;
		}
	}
	
		std::cout << "================= Unique words found " << uniq.size() << " =========== " << std::endl;
		ofstream out1 ("uniq.txt");
		for (std::set<std::string>::iterator IT = uniq.begin(); IT != uniq.end(); IT++) {
			//std::cout << *IT << std::endl;
			out1 << *IT << "\n";
		}
		out1.close();

		std::cout << "================= Duplicate words found " << duplicates.size() << " =========== " << std::endl;
		ofstream out2 ("duplicates.txt");
		for (std::set<std::string>::iterator IT = duplicates.begin(); IT != duplicates.end(); IT++) {
		//	std::cout << *IT << std::endl;
			out2 << *IT << "\n";
		}
 		out2.close();

	
return 0;
};
