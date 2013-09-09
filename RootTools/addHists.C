{
	const Int_t N = ;   // numebr of root files to loop over
	TFile* ff[N];
	TH1F* hist[N];
	char filename[200];
	char path[200];
	
	for (int i = 0; i< N; i++) {
		sprintf(filename,"PhoData_%i.root",i); 		// common file name differ by int(job section)
		ff[i] = new TFile(filename);
		gDirectory->cd("/Ana/TestModule/Hist/PhosWithStubs/");   // path to the hist interested in
		hist[i] = (TH1F*) gROOT->FindObject("PhosVsRun");      //hist to draw
		if (i>0) hist[0]->Add(hist[i]);
	}

	//TGaxis::MaxDigits(2);
	//hist[0]->Rebin(1000);
	
	hist[0]->Draw();

	//this par tis to get sum of a range of bin contents
	double xmin = hist[0]->GetXaxis()->GetXmin();
	double xmax = hist[0]->GetXaxis()->GetXmax();
	double binsize = hist[0]->GetBinWidth(1);
	int Nbins = (xmax - xmin)/binsize;

	int sum = 0;

	for (int i=1 ; i < Nbins; i++) {
		if (hist[0]->GetBinLowEdge(i)<190851) continue;
		sum+= hist[0]->GetBinContent(i);
	}

	std::cout << "pho= " << sum << std::endl;



}
