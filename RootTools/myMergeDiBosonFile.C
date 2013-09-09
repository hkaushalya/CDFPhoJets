#include <iomanip>
#include <iostream>
#include <fstream>
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TFolder.h"
#include <Stntuple/Stntuple/val/stntuple_val_functions.hh>

void myMergeDiBosonFile()
{
  gROOT->Reset();
  char fileList[1200];

  char f1[100]="results/myhisto_DiBosonAna_PythWW_it0sww_MetStudy_022309.root";
  char f2[100]="results/myhisto_DiBosonAna_PythWZ_it0swz_MetStudy_022309.root";
  char f3[100]="results/myhisto_DiBosonAna_PythZZ_it0szz_MetStudy_022309.root";

  sprintf(fileList,"%s %s %s",f1,f2,f3);

  merge_stn_hist(fileList,"results/myhisto_Pythia_DiBoson_MetStudy_standardCuts_NOdPhiClosest_031809.root");

  return;
}
