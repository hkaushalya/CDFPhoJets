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

#include <iostream>
#include <string>
#include "TFolder.h"
#include <fstream>
#include <vector>
#include <sstream>
#include "Stuple.hh"

void StdErr(const std::string& fileName, const unsigned int lineNumber,
								const int severity, const std::string& message);
void StdOut(const std::string& fileName, const unsigned int lineNumber,
								const int severity, const std::string& message);

std::string ToStr(float num);

void CorrectPhotonEnergy(Stuple& stuple, int JES);

#endif
