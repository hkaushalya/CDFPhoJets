TCanvas *THStackDraw() {
// Example of stacked histograms: class THStack
//
//  Author: Rene Brun
   
   THStack *hs = new THStack("hs","Stacked 1D histograms");
   //create three 1-d histograms
   TH1F *h1st = new TH1F("h1st","test hstack",100,-4,4);
   h1st->FillRandom("gaus",20000);
   h1st->SetFillColor(kRed);
   h1st->SetMarkerStyle(21);
   h1st->SetMarkerColor(kRed);
   hs->Add(h1st);
   TH1F *h2st = new TH1F("h2st","test hstack",100,-4,4);
   h2st->FillRandom("gaus",15000);
   h2st->SetFillColor(kBlue);
   h2st->SetMarkerStyle(21);
   h2st->SetMarkerColor(kBlue);
   hs->Add(h2st);
   TH1F *h3st = new TH1F("h3st","test hstack",100,-4,4);
   h3st->FillRandom("gaus",10000);
   h3st->SetFillColor(kGreen);
   h3st->SetMarkerStyle(21);
   h3st->SetMarkerColor(kGreen);
   hs->Add(h3st);
   
   TCanvas *cst = new TCanvas("cst","stacked hists",10,10,700,700);
   cst->SetFillColor(41);

	TPad *p1 = new TPad("log","Log Plot",0.01,0.4,0.99,0.99); 	//in NDC cdts
	p1->SetFillColor(kYellow);
	p1->SetBorderMode(0);
	p1->SetLineColor(kBlack);
	p1->Draw();
	p1->cd();

   hs->Draw();


	cst->cd();
	
	TPad *p2 = new TPad("ratio","Ratio Plot",0.01,0.01,0.99,0.4); 	//in NDC cdts
	p2->SetFillColor(kBlue);
	p2->SetBorderMode(2);
	p2->Draw();
   gPad->SetGrid();
	p2->cd();
	TH1F *rat = new TH1F("rat","rat",80,-4,4);
	rat->FillRandom("gaus",100);
	rat->Draw();

   return cst;
}
