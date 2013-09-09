/*
what this should do
have all the photons id cuts
sets tight id for the superphotons in InitSuperPhotons module

todo
return id-word for a given TStnPhoton
create histos for all cut values
*/


#ifndef FREEFUNCTIONS_HH
#define FREEFUNCTIONS_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#include <iostream>
#include <string>
#include "Stntuple/data/TCalTower.hh"
#include "TFolder.h"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/loop/TStnModule.hh"
#include <fstream>
#include <vector>
#include <sstream>
#include "Stntuple/obj/TGenpBlock.hh"
#endif


bool SortTowersByEnergy(TCalTower* c1, TCalTower* c2);
void StdOut(const std::string& msg, const std::string& moduleName);
void StdOut(const std::string& msg, const double num, const std::string& moduleName);
void StdErr(const std::string& fileName, const unsigned int lineNumber,
								const int severity, const std::string& message);
void StdOut(const std::string& fileName, const unsigned int lineNumber,
								const int severity, const std::string& message);
void ToFile(const ofstream* ofFile, const std::string& fileName, const unsigned int lineNumber,
								const int severity, const std::string& message);
void ToFile(ofstream* ofFile, const std::string& msg, const std::string& moduleName);
void ToFile(ofstream* ofFile, const std::string& msg, const double num, const std::string& moduleName);
void ToFile(const ofstream* ofFile, const std::string& message);
void ToFile(ofstream* ofFile, const char* message);
TFolder* GetHistFolder(TStnModule* const mod, const std::string& folderName,
									const std::string& folderTitle);
void PrintInBinaryForm(const int& number);
void binary(const int& number);
double DelPhi(double p1, double p2);
std::string ToStr(float num);

TLorentzVector FindMatchingHEPGPar(TGenpBlock* fGenpBlock, const TLorentzVector tlObj,
									const float fDelR, const int iPDGcode, const int iStatus=0,
									const int iLoopTo=-1, const bool bDebug=0);
TLorentzVector FindHEPGPar(TGenpBlock* fGenpBlock, const int iPDGcode, 
									const int iStatus, const int iLoopTo, 
									const bool bDebug);

bool CompareTLorentzVectorPt(TLorentzVector tl1, TLorentzVector tl2);
#endif
