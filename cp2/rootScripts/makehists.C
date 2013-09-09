{
gROOT->Reset();

TChain inChain("STNTUPLE");

//inChain.Add("/data/nbay02/a/samantha/STNTUPLES/CP2/P20_a267488_267547.stn");
//inChain.Add("/data/nbay02/a/samantha/STNTUPLES/CP2/P20_a267548_267607.stn");
//inChain.Add("/data/nbay02/a/samantha/STNTUPLES/CP2/P20_a267548_267607.stn_1");
//inChain.Add("/data/nbay02/a/samantha/STNTUPLES/CP2/P20_a267548_267607.stn_2");

//inChain.Add("/data/nbay02/a/samantha/STNTUPLES/CP2/P19_a265364_265413.stn");
//inChain.Add("/data/nbay02/a/samantha/STNTUPLES/CP2/P19_a265364_265413.stn_1");
//inChain.Add("/data/nbay02/a/samantha/STNTUPLES/CP2/P19_a265364_265413.stn_2");


inChain.Add("/data/nbay02/a/samantha/STNTUPLES/CP2/a256904.stn");
inChain.Add("/data/nbay02/a/samantha/STNTUPLES/CP2/a256904.stn_1");

TStnAna x(&inChain);
TH1::AddDirectory(1);
x.GetInputModule()->SetPrintLevel(1);   // print file name as they are opened
gSystem->CompileMacro("cp2.cc","=");
 cp2* mycp2 = new cp2;
  mycp2->SetNewDataFileName("a256904codetest.txt");
  mycp2->SetLastDBDataFileName("cp2traler-256904_LastDBEntry.txt");
  mycp2->SetLastDBLer(100);				//required
x.AddModule(mycp2);
//x.SetNEventsToReport(10000);
x.Run();
//gROOT->ProcessLine(".q");

}
