{
gROOT->Reset();
                                               // database ler1 was teststand
  float lerdec = 219.87;  // December             database ler2
  float ler2 = 211.0635;  // 195624 March 23rd    database ler3
  float ler3 = 209.4476;  // 196367  April 10th
  float ler4 = 207.6802;  // 196664  April 17th
  float ler5 = 206.9211;  // 196989  April 26th   database ler4
  float ler6 = 207.4479;  // 197657  May 7th
  float ler7 = 207.1190;  // 198613  May 26th
  float ler8 = 203.1341;  // 199189  June 10th    database ler5
  float ler9 = 204.5386;  // 200027  June 26th
  float ler10 = 209.1750;  // 200710  July 8th
  float ler11 = 202.4663;  // 201132  July 16th
  float ler12 = 199.9296;  // 201658  July 24th   database ler6 never put in
  float ler14 = 204.9342;   // 203348  August 28th
  float ler15 = 204.3491;  // 204550  Sep 20th
  float ler16 = 203.2912;   // 206723   Nov 7th JKim version
  //float ler16 = 204.47;    //206723  Ai Nagano version
  float ler17 = 202.5118;   //218658 June 19 2006
  float ler18 = 200.9812;   //222957 Oct 1 2006  database ler7
  float ler19 = 198.36;    //232522 Jan 15 2007 
  float ler20 = 206.555;   //233764 Feb 7  2007
  float ler21 = 208.1363;   //246125 Aug 1 2007
  float ler22 = 192.61;   //255202 Dec 17 2007
  float ler23 = 189.827;   //256904 Jan 29, 2008  database ler8
  float ler24 = 192.254;    // p21:267719 - 271000

  float ler25 = 191.447;    // p22:271071_272214
  float ler26 = 192.359; 	//p23:271071 - 272214
  
  float lernew = ler26;   //change this for new DB entry 

  Float_t cp2lerin[48][54];
  float bgg,error;
  int err,id0,ipad,iwed,iside;

  Float_t cp2lerin1[48][54];
  float bgg1,error1;
  int err1,id01,ipad1,iwed1,iside1;
  int err1;
  FILE *fpp1=fopen("cp2traler-p23_272470-274055.txt","r");    // change this
  FILE *fppout=fopen("cp2ler11.txt","w");           // change this
  for(int i=0;i<2592;i++){
    err1 = fscanf(fpp1,"%d %f %f\n",&id0,&bgg1,&error1);
    if(err1<0) break;
    ipad = (id0 & 0x0000003f);
    iwed = ((id0 >> 8) & 0x0000001f);
    iside = ((id0 >> 6) & 0x00000001);
    cp2lerin1[iwed+iside*24][ipad] = bgg1*lernew/lerdec;
    err1 = fprintf(fppout,"%d %f %f\n",id0,bgg1*lerdec/lernew,error1);
    if(err1<0) break;
    if(i<30) cout<<bgg1<<" "<<cp2lerin1[iwed+iside*24][ipad]<<" "<<bgg1*lerdec/lernew<<endl;
  }
  fclose(fpp1);
  fclose(fppout);

}
