/********************************************************************
 * This will overlay the photon-jet pt balance plots
 * for the signal mc/ data sideband /data iso-added sideband
 * and data iso-added reweighed sideband
 ********************************************************************
 * $Id: MakePhoJetPtBalanceHist.C,v 1.2 2010/02/13 21:40:04 samantha Exp $
 *
 * $Log: MakePhoJetPtBalanceHist.C,v $
 * Revision 1.2  2010/02/13 21:40:04  samantha
 * MODIFIED: 1. Removed the reweighted histogram. Just making pho mc, data
 * sideband and iso-added data sideband are plotted. see elog#1158 for the
 * results.
 *
 * Revision 1.1  2010/02/13 19:34:45  samantha
 * Created this to make the photon-jet pt balance plots to see if adding
 * the iso to the sideband photon improves the bakground modeling
 * done using the sideband.
 *
 *
 *******************************************************************/

void MakePhoJetPtBalanceHist() {
std::string datapath("/Hist/SIDEBAND/1Jet/PhotonLeadJet/PhoJetPtBalance");
std::string mcpath("/Hist/CENTRAL/1Jet/PhotonLeadJet/PhoJetPtBalance");
std::string title("#gamma^{E_{T}>40GeV,|#eta|<1.1}+==1 Jet^{E_{T}>15GeV,|#eta|<3.2} : #Delta#Phi(#gamma, jet)>2.7 && veto extra jets with (lev6) P_{T}>8GeV");

gStyle->SetOptStat("eimr");
TFile *g = new TFile("PhoJets_phomc_gammaEt40_PhoJetBalance.root");
TH1* hnoiso_mc = (TH1*) g->Get(mcpath.c_str());
assert (hnoiso_mc != NULL && "mc hist is null");
hnoiso_mc->SetLineColor(kBlack);
hnoiso_mc->SetLineWidth(2);
TCanvas *c = new TCanvas();
hnoiso_mc->Draw();
gPad->Update();

TPaveStats *tps = (TPaveStats*) hnoiso_mc->FindObject("stats");
assert (tps != NULL && "tps is null");
tps->Print();
double X1 = tps->GetX1NDC();
double Y1 = tps->GetY1NDC();
double X2 = tps->GetX2NDC();
double Y2 = tps->GetY2NDC();


TFile *f = new TFile("PhoJets_data_gammaEt40_PhoJetBalance_IsoAddedSideband.root");
TH1* hiso_data = (TH1*) f->Get(datapath.c_str());
hiso_data->SetLineColor(kRed);
hiso_data->SetLineWidth(2);
hiso_data->Scale(hnoiso_mc->Integral()/(double) hiso_data->Integral());
c1 = new TCanvas();
hiso_data->Draw();
gPad->Update();
//gPad->ls();
//TVirtualPad *pad1 = (TVirtualPad*) c1->GetPad(0);
//TPaveStats *tps1 = (TPaveStats*) pad1->GetPrimitive("stats");
TPaveStats *tps1 = (TPaveStats*) hiso_data->FindObject("stats");
tps1->SetTextColor(kRed);
tps1->SetLineColor(kRed);
//tps1->Print();
tps1->SetX1NDC(X1);
tps1->SetX2NDC(X2);
tps1->SetY1NDC(Y1-(Y2-Y1));
tps1->SetY2NDC(Y1);
tps1->Print();
tps1->Draw("same");
gPad->Update();

TFile *f1 = new TFile("PhoJets_data_gammaEt40_PhoJetBalance.root");
TH1* hnoiso_data = (TH1*) f1->Get(datapath.c_str());
hnoiso_data->SetLineColor(kBlue);
hnoiso_data->SetLineWidth(2);
hnoiso_data->Scale(hnoiso_mc->Integral()/(double) hnoiso_data->Integral());
c2 = new TCanvas();
hnoiso_data->Draw();
gPad->Update();
TPaveStats *tps2 = (TPaveStats*) hnoiso_data->FindObject("stats");
tps2->SetTextColor(kBlue);
tps2->SetLineColor(kBlue);

tps2->SetX1NDC(X1);
tps2->SetX2NDC(X2);
tps2->SetY1NDC(Y1 - 2 * (Y2-Y1));
tps2->SetY2NDC(Y1 - (Y2-Y1));


c3 = new TCanvas();
hnoiso_mc->SetTitle(title.c_str());
hnoiso_mc->Draw();
hiso_data->Draw("same");
hnoiso_data->Draw("same");
tps->Draw("same");
tps1->Draw("same");
tps2->Draw("same");
//tps3->Draw("same");

TLegend *tl = new TLegend(X1-2 * (X2-X1), Y1-0.15 ,X1,Y2 - (Y2-Y1)/2.,"LEGEND","NDC");
tl->SetLineColor(kBlack);
tl->AddEntry(hnoiso_mc,"tight #gamma MC (no iso)");
tl->AddEntry(hiso_data,"sideband #gamma DATA (iso added)");
tl->AddEntry(hnoiso_data,"sideband #gamma DATA (no iso)");
tl->Draw("same");


return;

}

