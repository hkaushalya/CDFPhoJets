{

	TF1* f1 = new TF1("f1","300*exp(-0.5*x)",0,100);
	f1->SetLineColor(kBlue);
	TF1* f2 = new TF1("f2","300*exp(-2*x)",0,100);
	f2->SetLineColor(kRed);

	new TCanvas();
	gPad->SetGridx();
	gPad->SetGridy();
	//gPad->SetLogy();
	f1->Draw();
	f2->Draw("same");
}
