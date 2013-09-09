#include <sstream>
#include "samantha/Pho/EmTimeCorrection.hh"
#include "TFile.h"
#include "TH1.h"
//ClassImp("EmTimeCorrection")

/*-------------------------------------------------------------------*/
EmTimeCorrection::EmTimeCorrection(OperMode mode, std::string varname)
{
  _opmode 		= mode;
  _varName 		= varname;
  _curentRun 	= -1;
  _curentValue = -999.;
  _curentHist 	= NULL;
  _histNbins 	= 0;
  _histXmin 	= 0.;
  _histXmax 	= 1.;
  _printCount 	= 0;
  _hasCor 		= false;
  _rfinder 		= EXACT;
}

/*-------------------------------------------------------------------*/
void EmTimeCorrection::setHisto(int n, float x1, float x2)
{
  _histNbins = n;
  _histXmin  = x1;
  _histXmax  = x2;
}

/*-------------------------------------------------------------------*/
void EmTimeCorrection::addEvent(int runNumber, float val)
{
  if (_opmode != WRITEM) return;

  if (runNumber != _curentRun) {
    _curentRun = runNumber;

    RHMI hpoint = _runHists.find(runNumber);
    if (hpoint == _runHists.end()) {
      std::stringstream hname;
      hname<<_varName<<"_"<<runNumber;
      //std::cout<<"created new run hist\n:"<<hname.str()<<"\n";

      TH1F* nhist = new TH1F(hname.str().c_str(),hname.str().c_str(),
			     _histNbins,_histXmin,_histXmax);
      _runHists.insert(RunHistPair(runNumber,nhist));
      _curentHist = nhist;
    } else {
      _curentHist = hpoint->second;
    }

  }
  _curentHist->Fill(val);
}

/*-------------------------------------------------------------------*/
void EmTimeCorrection::endJob()
{
  if (_opmode != WRITEM) return;
  if (_runHists.size() < 1) return;

  std::vector<int> allRuns;
  std::vector<float> allVals;
  for (RHMI ii = _runHists.begin(); ii != _runHists.end(); ii++){
    allRuns.push_back(ii->first);
    allVals.push_back(ii->second->GetMean());
  }
  std::sort(allRuns.begin(),allRuns.end());
  std::sort(allVals.begin(),allVals.end());

  int minRun = allRuns[0];
  int maxRun = allRuns[allRuns.size()-1];
  float minVal = allVals[0];
  float maxVal = allVals[allVals.size()-1];

  TH1* runDepMean = new TH1F("Mean Vs Run",_varName.c_str(),maxRun-minRun+1,
			     minRun-0.5,maxRun+0.5);
  TH1* runDepRms = new TH1F("Rms Vs Run",_varName.c_str(),maxRun-minRun+1,
			    minRun-0.5,maxRun+0.5);

  TH1* mean1d = new TH1F("Means",_varName.c_str(),100,minVal-0.1,maxVal+0.1);
  TH1* rms1d = new TH1F("Rmss",_varName.c_str(),100,minVal-0.1,maxVal+0.1);

  std::ofstream theFile;
  theFile.open(_varName.c_str());
  if (! theFile.is_open() ){
	  StdOut(__FILE__,__LINE__,3,"File " + _varName + "did not open");
    //std::cout<<" : file "<<_varName <<" did not open\n";
    return;
  }

  std::string fname = _varName + ".root";
  TFile*  myFile = new TFile(fname.c_str(),"RECREATE","EmTimeCorrection");
  myFile->cd();
  
  int totEvents = 0;

  std::vector<int>::iterator vri = allRuns.begin();
  for ( ; vri != allRuns.end(); vri++) {
    int run = *vri;
    TH1* hist = (_runHists.find(run))->second;

    totEvents = totEvents + static_cast<int>(hist->GetEntries());

    float mean = hist->GetMean();
    float rms = hist->GetRMS();
    
    int bin1 = hist->FindBin(mean-2.*rms);
    int bin2 = hist->FindBin(mean+2.*rms);
    
    if (hist->Integral(bin1,bin2) < 20) continue;
    
    TF1* fg = (TF1*)gROOT->GetFunction("gaus");
    hist->Fit(fg,"LQ","",mean-3*rms,mean+3*rms);
    mean = fg->GetParameter(1);
    rms = fg->GetParameter(2);
    
    hist->Fit(fg,"LQ","",mean-3*rms,mean+3*rms);
    mean = fg->GetParameter(1);
    float ermean = fg->GetParError(1);
    rms = fg->GetParameter(2);
    float errms = fg->GetParError(2);
    
    int ibin = runDepMean->FindBin(run);
    theFile<<run<<" "<<fg->GetParameter(1)<<" "<<fg->GetParError(1)
	   <<" "<<fg->GetParameter(2)<<" "<<hist->GetEntries()<<"\n";
    
    runDepMean->SetBinContent(ibin,mean);
    runDepMean->SetBinError(ibin,ermean);
    
    runDepRms->SetBinContent(ibin,rms);
    runDepRms->SetBinError(ibin,errms);
    
    mean1d->Fill(mean,1./(ermean*ermean));
    rms1d->Fill(rms,1./(errms*errms));
    
    hist->Write();
  }
  runDepMean->Write();runDepRms->Write();
  mean1d->Write();rms1d->Write();
  
  myFile->Close();
  delete myFile;

  theFile.close();

  std::cout<<"\t"<<_varName<<" accumulated "<<totEvents<<" events\n";
}

/*-------------------------------------------------------------------*/
bool EmTimeCorrection::readInFile()
{
  if (_opmode == WRITEM) return false;
  
  if (_runCors.size() < 1) {
    std::ifstream theFile;
    theFile.open(_varName.c_str());
    if (! theFile.is_open() ) {
      std::cout<<" did not open file "<<_varName<<" for reading\n";
	  	StdOut(__FILE__,__LINE__,3,"File " + _varName + "did not open for reading.");
      return false;
    }
    while (! theFile.eof() ) {     
      int run, nevents;
      float mean, ermean,rms;
      
      theFile >>run>>mean>>ermean>>rms>>nevents;
      
      _runCors.insert(RunCorPair(run,mean));
      if (theFile.eof()) {break;}
    }
    theFile.close();
  }

  return true;
}

/*-------------------------------------------------------------------*/
bool EmTimeCorrection::hasRunCorrection(int run)
{
  if(! readInFile() ) return false;
  
  bool newRun = false;
  if (_curentRun != run) {
    newRun = true;
    _curentRun = run;
    _curentValue = -999.;
    _hasCor = false;
  }
  if (newRun) {
    RCVI rp = _runCors.find(run);
    if (rp == _runCors.end()) {
      
      if (_rfinder != NEAREST) return _hasCor;
      
     	RCVI rpl = _runCors.lower_bound(run);rpl--;
	   RCVI rph = _runCors.lower_bound(run);
      int clrun = 0;
      if (rpl != _runCors.end() ) {
			clrun = rpl->first;
      }
      if (rph != _runCors.end() ) {
			clrun = rph->first;
			if(fabs(rph->first - run) < fabs(clrun - run)) clrun = rph->first;
      }
      if (clrun != 0) {
			_curentValue = _runCors.find(clrun)->second;
			_printCount = 0;
			_hasCor = true;
      } else {
			if(_printCount < 1)
			std::cout<< __FILE__ << ":" << __LINE__ << ":" << _varName<<": no correction for run "<<run<< std::endl;
			_printCount++;
      }
    } else {
      _curentValue = rp->second;
      _printCount  = 0;
      _hasCor      = true;
    } 
  }

  return _hasCor;
}

/*-------------------------------------------------------------------*/
float EmTimeCorrection::getRunCorrection(int run)
{
  if ( _opmode == DUMMY) return 0.;
  if (! hasRunCorrection(run) ) return -999.;
  return _curentValue;
}
