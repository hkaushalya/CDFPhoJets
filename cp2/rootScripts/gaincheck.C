{
gROOT->Reset();

TChain inChain("STNTUPLE");

//Samarium
inChain.Add("../P19_a265914_265923.stn");
inChain.Add("../P19_a265924_265933.stn");
inChain.Add("../P19_a265934_265943.stn");
inChain.Add("../P19_a265934_265943.stn_1");
inChain.Add("../P19_a265944_265953.stn");
inChain.Add("../P19_a265954_265955.stn");
inChain.Add("../P19_a265954_265963.stn");
inChain.Add("../P19_a265956_265957.stn");
inChain.Add("../P19_a265958_265959.stn");
inChain.Add("../P19_a265960_265961.stn");
inChain.Add("../P19_a265962_265963.stn");


TStnAna x(&inChain);
TH1::AddDirectory(1);
x.GetInputModule()->SetPrintLevel(1);   // print file name as they are opened
gSystem->CompileMacro("cp2.cc","=");
 cp2* mycp2 = new cp2;
  mycp2->SetLogFile("P20HV_A265914_265963.txt");
x.AddModule(mycp2);
//x.SetNEventsToReport(10000);
x.Run();
gROOT->ProcessLine(".q");

}
