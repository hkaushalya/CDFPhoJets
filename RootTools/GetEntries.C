#include<iostream>
#include<fstream>
#include<TFile.h>
#include<TTree.h>
#include<vector>
#include<string>

void GetEntries(const std::string& filename) {

ifstream fp(filename.c_str());

if (!fp) {
	std::cout << "Input file not found!" <<std::endl;
	return;
}

std::vector<std::string> file;
std::string str;

while(fp >> str) {
	file.push_back(str);
}

for (unsigned i=0; i < file.size(); i++) std::cout << i << "= "<< file[i] << std::endl;

	Double_t sument =0;
	for (unsigned i=0; i < file.size(); i++) {
		TFile f(file[i].c_str(),"READ");
		if (!f.IsZombie()) {
			Double_t ent = dynamic_cast<TTree*>(f.Get("STNTUPLE"))->GetEntries();
			sument += ent;
			std::cout << file[i] << "  entries = " << ent <<std::endl;
			f.Close();
		} else std::cout << file[i] << " NOT FOUND!" << std::endl;
	}
	std::cout << "  Sum Entries = " << sument <<std::endl;

}
