#include <iostream>
#include <vector>
#include <algorithm>
#include "TLorentzVector.h"

bool comp(int i, int j) { return (i<j); }
bool compTL(TLorentzVector tl1, TLorentzVector tl2) { return (tl1.Pt()>tl2.Pt()); }


void StdSort()
{
int myints[] = {32,71,12,45,26,80,53,33};
std::vector<int> myvector (myints, myints+8);
std::vector<int>::iterator it;
std::sort (myvector.begin(), myvector.end(), comp);
std::cout << "myvector contains:";
  for (it=myvector.begin(); it!=myvector.end(); ++it)
	      std::cout << " " << *it;

  std::cout << std::endl;

  std::vector<TLorentzVector> tlVec;
  TLorentzVector tl (4,5,10,20);
  TLorentzVector tl2 (2,3,1,2);
  TLorentzVector tl3 (20,30,1,2);
  tlVec.push_back(tl);
  tlVec.push_back(tl2);
  tlVec.push_back(tl3);

  std::vector <TLorentzVector>::iterator it2;
  for (it2 = tlVec.begin(); it2 != tlVec.end(); ++it2)
	  std::cout << " " << (*it2).Pt();
  std::cout << std::endl;

  std::sort (tlVec.begin(), tlVec.end(), compTL);

  for (it2 = tlVec.begin(); it2 != tlVec.end(); ++it2)
	  std::cout << " " << (*it2).Pt();
  std::cout << std::endl;

}
