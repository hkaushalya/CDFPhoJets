#include <list>

/*
 * Run on ths output of the scan_caf.C which generated the minprob_xx.txt.
 * to derive the observed and expected limits.
 * 04-30-2010 - Sam 
 */

void getlimits(){

  TH1D *hminprob = new TH1D("hminprob","hminprob",100,0,0.3);
  
  list<double> vminp;
  int n = 0;
  ifstream in;
  double x;
  for (int i = 1; i<=200; i++){
    in.open(Form("scanresult_1gev/minprob_%d.txt",i));
    while (1){
      in>>x;
      if (!in.good()) break;
      hminprob->Fill(x);
      vminp.push_back(x);
      n++;
    }
    in.close();
    in.clear();
  }
  hminprob->DrawCopy();
  gPad->SetLogy();
  vminp.sort();
  double p1s=0;
  double m1s=0;
  double s3s=0;
  //for (int i = 0; i<n; i++) cout<<vminp[i]<<endl;
  list<double>::iterator p = vminp.begin();
  int ind = 0;
  while(p!=vminp.end()){
    //cout<<*p<<endl;
    ind++; 
    if (ind==int((1-0.683)/2*n)) {
      p1s = *p;
      //cout<<*(--p)<<" "<<*(++p)<<" "<<*(++p)<<" "<<*(++p)endl;
    }
    if (ind==int(((1-0.683)/2+0.683)*n)) m1s = *p;
    if (ind==int((1-0.997)*n)) s3s = *p;
    p++;
  }
  cout<<p1s<<" "<<m1s<<" "<<s3s<<endl;
//  cout<<int((1-0.683)/2*n)<<endl;
//  cout<<int(((1-0.683)/2+0.683)*n)<<endl;
//  cout<<int((1-0.997)*n)<<endl;
}
