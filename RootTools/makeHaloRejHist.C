{

TFile *f = new TFile("PhoJets_data.root");
TH1* hb4 = (TH1*) (f->Get("Hist/haloPhiWedge_1j_b4"));
TH1* ha4 = (TH1*) (f->Get("Hist/haloPhiWedge_1j_a4"));

assert (hb4 != NULL && "hb4 not found");
assert (ha4 != NULL && "ha4 not found");

//hb4->SetLineColor(kBlue);
hb4->SetLineWidth(2);
//ha4->SetLineColor(kRed);
ha4->SetLineWidth(2);
hb4->SetMinimum(0);
ha4->SetMinimum(0);

hb4->SetTitle("Beam Halo photon candidate #phi wedge distribution before cuts;photon #phi wedge;Events");
ha4->SetTitle("Beam Halo photon candidate #phi wedge distribution after cuts;photon #phi wedge;Events");

gStyle->SetOptStat("e");
new TCanvas();
hb4->Draw();
new TCanvas();
ha4->Draw();

}
