#include <iostream>
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
#include "TMath.h"

NewJetList::NewJetList()
{
};

//-------------------------------------------------------------------
NewJetList::NewJetList(Stuple& stuple, int iObjType)		//iObjType = 0 (photon), 1(electrons
{
	//if opt= ele , put all the electrons back in the jet list
	//if opt= pho, put all the photons back in the jet list
//	std::cout << "===============================" <<std::endl;

	Stuple stuple_temp = stuple;		//make a temp copy
	
/*	
	std::cout << "Njets=" << stuple.jet_num << std::endl;
	for (unsigned int i =0; i < stuple.jet_num; ++i) {
		std::cout <<  "\t" << i <<"=" << stuple.jet_Pt[i] << std::endl;
	}

	std::cout << "NTele=" << stuple.ele_Ntight << std::endl;
	for (unsigned int i =0; i < stuple.ele_num; i++) {
			if (stuple.ele_TightId[i] != 0) continue;		//select electrons only
			std::cout << "\t"<< i << "=" << stuple.ele_Etc[i] << std::endl;
	}
*/
	// here is what this loop does
	// we are adding the obj back to jet list
	// two scenarios
	// 1. a perfect match (that jet is completely removed, so all we need is to add this obj back
	// 2. a partial match (i.e. after removing obj, there is still enough energy in the cluster to form a jet. so we need to find it and add this obj.

	switch (iObjType) {
		case 1:
		{
			std::vector<int> EleNoMatch;
			bool bNeedSort = false;
			for (unsigned int i =0; i < stuple.ele_num; i++) {
				if (stuple.ele_TightId[i] != 0) continue;		//select electrons only
//				std::cout << "ele,  jetindex, etc=" << i << "," << stuple.ele_matchJetIndex[i] << "\t" << stuple.ele_Etc[i] << std::endl;
				if (stuple.ele_matchJetIndex[i] >=0) {
					// std::cout << "jets=" << stuple.jet_num << std::endl;
					bool bFoundMatch = false;
					for (unsigned int j=0; j < stuple.jet_num; j++) {
						//std::cout <<"\tlooking at jet index = " << stuple.jet_Index[j] << std::endl;
						if (stuple.ele_matchJetIndex[i] != stuple.jet_Index[j]) continue;
						bFoundMatch = true;
						bNeedSort = true;
//						std::cout << "Run,Evt::"<< stuple.evt_RunNumber << "," << stuple.evt_EventNumber << std::endl;
//						std::cout << "\tele,E = " << stuple.ele_Index[i] << "," << stuple.ele_E[i]
//									<<  "jet,E = " << stuple.jet_Index[j] << "," << stuple.jet_E[j]
//									<< std::endl;
						stuple_temp.jet_E[j]  = stuple.ele_E[i] + stuple.jet_E[j];
						stuple_temp.jet_Px[j] = stuple.ele_Px[i] + stuple.jet_Px[j];
						stuple_temp.jet_Py[j] = stuple.ele_Py[i] + stuple.jet_Py[j];
						stuple_temp.jet_Pz[j] = stuple.ele_Pz[i] + stuple.jet_Pz[j];
						TLorentzVector vJet;
						vJet.SetPxPyPzE(stuple_temp.jet_Px[j],stuple_temp.jet_Py[j],stuple_temp.jet_Pz[j],stuple_temp.jet_E[j]);
						stuple_temp.jet_Pt[j] = vJet.Pt();
//						std::cout << "\tafter" << std::endl;
//						std::cout << "\tele,E = " << stuple.ele_Index[i] << "," << stuple.ele_E[i]
//								<<  "jet,E = " << stuple.jet_Index[j] << "," << stuple.jet_E[j]
//								<< std::endl;
						
						break;
					} // for

					if (!bFoundMatch) EleNoMatch.push_back(i);

				} else EleNoMatch.push_back(i);
	
			}	// for

			if (EleNoMatch.size()>0) { 
				// add all ele to jets
				AddEleToJetList(stuple_temp, EleNoMatch);
				bNeedSort = true;
			}
		
			if (bNeedSort) {
				ReorderJets(stuple_temp);
				stuple = stuple_temp;
			}
		
//			for (unsigned int i =0; i < stuple.jet_num; ++i) {
//				std::cout <<  "\t" << i <<"=" << stuple.jet_Pt[i] << std::endl;
//			}
			break;
		} // case 1


		case 2:
		{
			std::vector<int> PhoNoMatch;
			bool bNeedSort = false;
			for (unsigned int i =0; i < stuple.pho_num; i++) {
				if (stuple.pho_TightId[i] != 0) continue;		//select electrons only
//				 std::cout << "pho,  jetindex, etc=" << i << "," << stuple.pho_matchJetIndex[i] << "\t" << stuple.pho_Etc[i] << std::endl;
				if (stuple.pho_matchJetIndex[i] >=0) {
					// std::cout << "jets=" << stuple.jet_num << std::endl;
					bool bFoundMatch = false;
					for (unsigned int j=0; j < stuple.jet_num; j++) {
						//std::cout <<"\tlooking at jet index = " << stuple.jet_Index[j] << std::endl;
						if (stuple.pho_matchJetIndex[i] != stuple.jet_Index[j]) continue;
						bFoundMatch = true;
						bNeedSort = true;
//						std::cout << "Run,Evt::"<< stuple.evt_RunNumber << "," << stuple.evt_EventNumber << std::endl;
//						std::cout << "\tele,E = " << stuple.pho_Index[i] << "," << stuple.pho_E[i]
//									<<  "jet,E = " << stuple.jet_Index[j] << "," << stuple.jet_E[j]
//									<< std::endl;
						stuple_temp.jet_E[j]  = stuple.pho_E[i] + stuple.jet_E[j];
						stuple_temp.jet_Px[j] = stuple.pho_Px[i] + stuple.jet_Px[j];
						stuple_temp.jet_Py[j] = stuple.pho_Py[i] + stuple.jet_Py[j];
						stuple_temp.jet_Pz[j] = stuple.pho_Pz[i] + stuple.jet_Pz[j];
						TLorentzVector vJet;
						vJet.SetPxPyPzE(stuple_temp.jet_Px[j],stuple_temp.jet_Py[j],stuple_temp.jet_Pz[j],stuple_temp.jet_E[j]);
						stuple_temp.jet_Pt[j] = vJet.Pt();
//						std::cout << "\tafter" << std::endl;
//						std::cout << "\tele,E = " << stuple.pho_Index[i] << "," << stuple.pho_E[i]
//								<<  "jet,E = " << stuple.jet_Index[j] << "," << stuple.jet_E[j]
//								<< std::endl;
						
						break;
					} // for

					if (!bFoundMatch) PhoNoMatch.push_back(i);

				} else PhoNoMatch.push_back(i);
	
			}	// for

			if (PhoNoMatch.size()>0) { 
				// add all ele to jets
				AddPhoToJetList(stuple_temp, PhoNoMatch);
				bNeedSort = true;
			}
		
			if (bNeedSort) {
				ReorderJets(stuple_temp);
				stuple = stuple_temp;
			}
		
//			for (unsigned int i =0; i < stuple.jet_num; ++i) {
//				std::cout <<  "\t" << i <<"=" << stuple.jet_Pt[i] << std::endl;
//			}
			break;
		} // case 2



		
	} // if obj type == 0
	
}


void NewJetList::AddUnused(Stuple& stuple, std::vector<int>UsedPho, std::vector<int>UsedEle)
{
	Stuple stuple_temp = stuple;
	
	
	bool bNeedSort = false;
	
	if (UsedPho.size() < stuple.pho_num) {
		//add these phos
			std::vector<int> PhoNoMatch;
			for (unsigned int i =0; i < stuple.pho_num; i++) {
				if (IsUsed(i, UsedPho)) continue;
//				 std::cout << "pho,  jetindex, etc=" << i << "," << stuple.pho_matchJetIndex[i] << "\t" << stuple.pho_Etc[i] << std::endl;
				if (stuple.pho_matchJetIndex[i] >=0) {
					// std::cout << "jets=" << stuple.jet_num << std::endl;
					bool bFoundMatch = false;
					for (unsigned int j=0; j < stuple.jet_num; j++) {
						//std::cout <<"\tlooking at jet index = " << stuple.jet_Index[j] << std::endl;
						if (stuple.pho_matchJetIndex[i] != stuple.jet_Index[j]) continue;
						bFoundMatch = true;
						bNeedSort = true;
//						std::cout << "Run,Evt::"<< stuple.evt_RunNumber << "," << stuple.evt_EventNumber << std::endl;
//						std::cout << "\tele,E = " << stuple.pho_Index[i] << "," << stuple.pho_E[i]
//									<<  "jet,E = " << stuple.jet_Index[j] << "," << stuple.jet_E[j]
//									<< std::endl;
						stuple_temp.jet_E[j]  = stuple.pho_E[i] + stuple.jet_E[j];
						stuple_temp.jet_Px[j] = stuple.pho_Px[i] + stuple.jet_Px[j];
						stuple_temp.jet_Py[j] = stuple.pho_Py[i] + stuple.jet_Py[j];
						stuple_temp.jet_Pz[j] = stuple.pho_Pz[i] + stuple.jet_Pz[j];
						stuple_temp.jet_Pt[j]  = stuple.pho_Etc[i] + stuple.jet_Pt[j];
						TLorentzVector vJet;
						vJet.SetPxPyPzE(stuple_temp.jet_Px[j],stuple_temp.jet_Py[j],stuple_temp.jet_Pz[j],stuple_temp.jet_E[j]);
						stuple_temp.jet_Pt[j] = vJet.Pt();
//						std::cout << "\tafter" << std::endl;
//						std::cout << "\tele,E = " << stuple.pho_Index[i] << "," << stuple.pho_E[i]
//								<<  "jet,E = " << stuple.jet_Index[j] << "," << stuple.jet_E[j]
//								<< std::endl;
						
						break;
					} // for

					if (!bFoundMatch) PhoNoMatch.push_back(i);

				} else PhoNoMatch.push_back(i);
	
			}	// for

			if (PhoNoMatch.size()>0) { 
				AddPhoToJetList(stuple_temp, PhoNoMatch);
				bNeedSort = true;
			}
		
		
	}
	
	if (UsedEle.size() < stuple.ele_num) {
		std::vector<int> EleNoMatch;
		//add these eles
		for (unsigned int i =0; i < stuple.ele_num; i++) {
//				std::cout << "ele,  jetindex, etc=" << i << "," << stuple.ele_matchJetIndex[i] << "\t" << stuple.ele_Etc[i] << std::endl;
			if (IsUsed(i, UsedEle)) continue;
			if (stuple.ele_matchJetIndex[i] >=0) {
				// std::cout << "jets=" << stuple.jet_num << std::endl;
				bool bFoundMatch = false;
				for (unsigned int j=0; j < stuple.jet_num; j++) {
					//std::cout <<"\tlooking at jet index = " << stuple.jet_Index[j] << std::endl;
					if (stuple.ele_matchJetIndex[i] != stuple.jet_Index[j]) continue;
					bFoundMatch = true;
					bNeedSort = true;
//						std::cout << "Run,Evt::"<< stuple.evt_RunNumber << "," << stuple.evt_EventNumber << std::endl;
//						std::cout << "\tele,E = " << stuple.ele_Index[i] << "," << stuple.ele_E[i]
//									<<  "jet,E = " << stuple.jet_Index[j] << "," << stuple.jet_E[j]
//									<< std::endl;
					stuple_temp.jet_E[j]  = stuple.ele_E[i] + stuple.jet_E[j];
					stuple_temp.jet_Px[j] = stuple.ele_Px[i] + stuple.jet_Px[j];
					stuple_temp.jet_Py[j] = stuple.ele_Py[i] + stuple.jet_Py[j];
					stuple_temp.jet_Pz[j] = stuple.ele_Pz[i] + stuple.jet_Pz[j];
					stuple_temp.jet_Pt[j] = stuple.ele_Etc[i] + stuple.jet_Pt[j];
					TLorentzVector vJet;
					vJet.SetPxPyPzE(stuple_temp.jet_Px[j],stuple_temp.jet_Py[j],stuple_temp.jet_Pz[j],stuple_temp.jet_E[j]);
					stuple_temp.jet_Pt[j] = vJet.Pt();
//						std::cout << "\tafter" << std::endl;
//						std::cout << "\tele,E = " << stuple.ele_Index[i] << "," << stuple.ele_E[i]
//								<<  "jet,E = " << stuple.jet_Index[j] << "," << stuple.jet_E[j]
//								<< std::endl;
					
					break;
				} // for

				if (!bFoundMatch) EleNoMatch.push_back(i);

			} else EleNoMatch.push_back(i);

		}	// for

		if (EleNoMatch.size()>0) { 
			// add all ele to jets
			AddEleToJetList(stuple_temp, EleNoMatch);
			bNeedSort = true;
		}
	
	}
		
		
	if (bNeedSort) {
		ReorderJets(stuple_temp);
		stuple = stuple_temp;
	}

}



//-------------------------------------------------------------------
bool NewJetList::IsUsed(int ind, std::vector<int> list)
{
	for (unsigned i = 0; i < list.size(); ++i) {
		if (ind == list[i]) return true;
	}
	return false;
}



//-------------------------------------------------------------------
void NewJetList::AddEleToJetList(Stuple& stuple, std::vector<int> EleNoMatch)
{
		for (unsigned int i=0; i < EleNoMatch.size(); i++) {
			if (stuple.jet_num < stuple.Njet) {
				stuple.jet_E[stuple.jet_num] = stuple.ele_E[EleNoMatch[i]]; 
				stuple.jet_Px[stuple.jet_num] = stuple.ele_Px[EleNoMatch[i]]; 
				stuple.jet_Py[stuple.jet_num] = stuple.ele_Py[EleNoMatch[i]]; 
				stuple.jet_Pz[stuple.jet_num] = stuple.ele_Pz[EleNoMatch[i]]; 
				TLorentzVector vJet;
				vJet.SetPxPyPzE(stuple.jet_Px[stuple.jet_num],
									 stuple.jet_Py[stuple.jet_num],
									 stuple.jet_Pz[stuple.jet_num],
									 stuple.jet_E[stuple.jet_num]);
				stuple.jet_Pt[stuple.jet_num] = vJet.Pt();
//				std::cout << "new jet Pt= " << stuple.jet_Pt[stuple.jet_num] <<std::endl;
				stuple.jet_num++;
				//follwing is for safety. even if i call this again in the same module by mistake, this will prevent any duplication! hopfully!
				//stuple.ele_num--;
				//if (stuple.ele_TightId[i] == 0) stuple.ele_Ntight--;
				//if (stuple.ele_LooseId[i] == 0) stuple.ele_Nloose--;
			} else {
				std::cout <<__FILE__ <<"::" << __LINE__ << ":: Maximum number of jets " << stuple.Njet << " reached while adding stuff to jet list!" << std::endl;
				exit (1);
			}
		}

}

//-------------------------------------------------------------------
void NewJetList::AddPhoToJetList(Stuple& stuple, std::vector<int> PhoNoMatch)
{
		for (unsigned int i=0; i < PhoNoMatch.size(); i++) {
			if (stuple.jet_num < stuple.Njet) {
				stuple.jet_E[stuple.jet_num] = stuple.pho_E[PhoNoMatch[i]]; 
				stuple.jet_Px[stuple.jet_num] = stuple.pho_Px[PhoNoMatch[i]]; 
				stuple.jet_Py[stuple.jet_num] = stuple.pho_Py[PhoNoMatch[i]]; 
				stuple.jet_Pz[stuple.jet_num] = stuple.pho_Pz[PhoNoMatch[i]]; 
				TLorentzVector vJet;
				vJet.SetPxPyPzE(stuple.jet_Px[stuple.jet_num],
									 stuple.jet_Py[stuple.jet_num],
									 stuple.jet_Pz[stuple.jet_num],
									 stuple.jet_E[stuple.jet_num]);
				stuple.jet_Pt[stuple.jet_num] = vJet.Pt();
//				std::cout << "new jet Pt= " << stuple.jet_Pt[stuple.jet_num] <<std::endl;
				stuple.jet_num++;
				//follwing is for safety. even if i call this again in the same module by mistake, this will prevent any duplication! hopfully!
				//stuple.pho_num--;
				//if (stuple.pho_TightId[i] == 0) stuple.pho_Ntight--;
				//if (stuple.pho_LooseId[i] == 0) stuple.pho_Nloose--;
			} else {
				std::cout <<__FILE__ <<"::" << __LINE__ << ":: Maximum number of jets " << stuple.Njet << " reached while adding stuff to jet list!" << std::endl;
				exit (1);
			}
		}

}


//-------------------------------------------------------------------
void NewJetList::ReorderJets(Stuple& stuple)
{
	//when there are no jets to begin with
	if (stuple.jet_num == 0) return;
	
	int loopupto = stuple.jet_num - 1;

	while (loopupto !=0) {
		for (int i=0; i<loopupto;i++) {
			if (stuple.jet_Pt[i] < stuple.jet_Pt[i+1]) {
				Stuple stuple_temp;

				stuple_temp.jet_Index[i] 	= stuple.jet_Index[i];
				stuple_temp.jet_Pt[i] 		= stuple.jet_Pt[i];
				stuple_temp.jet_E[i] 		= stuple.jet_E[i];
				stuple_temp.jet_Px[i] 		= stuple.jet_Px[i];
				stuple_temp.jet_Py[i] 		= stuple.jet_Py[i];
				stuple_temp.jet_Pz[i] 		= stuple.jet_Pz[i];
				stuple_temp.jet_DetEta[i] 	= stuple.jet_DetEta[i];
				stuple_temp.jet_DetPhi[i] 	= stuple.jet_DetPhi[i];
				stuple_temp.jet_HadEm[i] 	= stuple.jet_HadEm[i];
				stuple_temp.jet_Emfr[i] 	= stuple.jet_Emfr[i];
				stuple_temp.jet_Ntowers[i] 	= stuple.jet_Ntowers[i];
				stuple_temp.jet_Ntracks[i] 	= stuple.jet_Ntracks[i];
				stuple_temp.jet_SeedIPhi[i] 	= stuple.jet_SeedIPhi[i];
				stuple_temp.jet_SeedIEta[i] 	= stuple.jet_SeedIEta[i];
				stuple_temp.jet_EmTime[i] 	= stuple.jet_EmTime[i];
				stuple_temp.jet_SecVtxTag[i]     = stuple.jet_SecVtxTag[i];
				stuple_temp.jet_SecVtxppb[i]     = stuple.jet_SecVtxppb[i] ;
				stuple_temp.jet_SecVtxnpb[i]     = stuple.jet_SecVtxnpb[i];
				stuple_temp.jet_SecVtxTrkmass[i] = stuple.jet_SecVtxTrkmass[i];
		

				stuple.jet_Index[i] 	= stuple.jet_Index[i+1];
				stuple.jet_Pt[i] 		= stuple.jet_Pt[i+1];
				stuple.jet_E[i] 		= stuple.jet_E[i+1];
				stuple.jet_Px[i] 		= stuple.jet_Px[i+1];
				stuple.jet_Py[i] 		= stuple.jet_Py[i+1];
				stuple.jet_Pz[i] 		= stuple.jet_Pz[i+1];
				stuple.jet_DetEta[i] = stuple.jet_DetEta[i+1];
				stuple.jet_DetPhi[i] = stuple.jet_DetPhi[i+1];
				stuple.jet_HadEm[i] 	= stuple.jet_HadEm[i+1];
				stuple.jet_Emfr[i] 	= stuple.jet_Emfr[i+1];
				stuple.jet_Ntowers[i] 	= stuple.jet_Ntowers[i+1];
				stuple.jet_Ntracks[i] 	= stuple.jet_Ntracks[i+1];
				stuple.jet_SeedIPhi[i] 	= stuple.jet_SeedIPhi[i+1];
				stuple.jet_SeedIEta[i] 	= stuple.jet_SeedIEta[i+1];
				stuple.jet_EmTime[i] 	= stuple.jet_EmTime[i+1];
				stuple.jet_SecVtxTag[i]     = stuple.jet_SecVtxTag[i+1];
				stuple.jet_SecVtxppb[i]     = stuple.jet_SecVtxppb[i+1] ;
				stuple.jet_SecVtxnpb[i]     = stuple.jet_SecVtxnpb[i+1];
				stuple.jet_SecVtxTrkmass[i] = stuple.jet_SecVtxTrkmass[i+1];


				stuple.jet_Index[i+1] 	= stuple_temp.jet_Index[i];
				stuple.jet_Pt[i+1] 		= stuple_temp.jet_Pt[i];
				stuple.jet_E[i+1] 		= stuple_temp.jet_E[i];
				stuple.jet_Px[i+1] 		= stuple_temp.jet_Px[i];
				stuple.jet_Py[i+1] 		= stuple_temp.jet_Py[i];
				stuple.jet_Pz[i+1] 		= stuple_temp.jet_Pz[i];
				stuple.jet_DetEta[i+1] 	= stuple_temp.jet_DetEta[i];
				stuple.jet_DetPhi[i+1] 	= stuple_temp.jet_DetPhi[i];
				stuple.jet_HadEm[i+1] 	= stuple_temp.jet_HadEm[i];
				stuple.jet_Emfr[i+1] 	= stuple_temp.jet_Emfr[i];
				stuple.jet_Ntowers[i+1] = stuple_temp.jet_Ntowers[i];
				stuple.jet_Ntracks[i+1] = stuple_temp.jet_Ntracks[i];
				stuple.jet_SeedIPhi[i+1] = stuple_temp.jet_SeedIPhi[i];
				stuple.jet_SeedIEta[i+1] = stuple_temp.jet_SeedIEta[i];
				stuple.jet_EmTime[i+1]   = stuple_temp.jet_EmTime[i];
				stuple.jet_SecVtxTag[i+1]     = stuple_temp.jet_SecVtxTag[i];
				stuple.jet_SecVtxppb[i+1]     = stuple_temp.jet_SecVtxppb[i] ;
				stuple.jet_SecVtxnpb[i+1]     = stuple_temp.jet_SecVtxnpb[i];
				stuple.jet_SecVtxTrkmass[i+1] = stuple_temp.jet_SecVtxTrkmass[i];
					
			}
		} // for
		loopupto--;
	}
	
	for (unsigned int i=0 ; i < stuple.jet_num; i++) {
		if (i+1 < stuple.jet_num) assert(stuple.jet_Pt[i] >= stuple.jet_Pt[i+1]);
	}
}


// NEW FUNCTIONS TO BE USED FROM 05-16-2008.

// ------------------- CENTRAL JET -----------------------------------------------------------------
void NewJetList::AddUnusedCentral(Stuple& stuple, std::vector<int>UsedPho, std::vector<int>UsedEle)
{
	Stuple stuple_temp = stuple;
	
	bool bNeedSort = false;
	
	// add the un-used photons to jet list
	if (UsedPho.size() <= stuple.pho_num) {			// do not check for size>0, I am adding the un-used one in this loop.
		//add these phos
			std::vector<int> PhoNoMatch;
			for (unsigned int i =0; i < stuple.pho_num; i++) {
				if (IsUsed(i, UsedPho)) continue;						//check if this is the photon that was picked as THE photon
				if (stuple.pho_matchJetIndex[i] >=0) {
					bool bFoundMatch = false;
					for (unsigned int j=0; j < stuple.jet_num; j++) {
						if (stuple.pho_matchJetIndex[i] != stuple.jet_Index[j]) continue;
						bFoundMatch = true;
						bNeedSort = true;
						stuple_temp.jet_E[j]  = stuple.pho_E[i] + stuple.jet_E[j];
						stuple_temp.jet_Px[j] = stuple.pho_Px[i] + stuple.jet_Px[j];
						stuple_temp.jet_Py[j] = stuple.pho_Py[i] + stuple.jet_Py[j];
						stuple_temp.jet_Pz[j] = stuple.pho_Pz[i] + stuple.jet_Pz[j];
						stuple_temp.jet_Pt[j]  = stuple.pho_Etc[i] + stuple.jet_Pt[j];
						TLorentzVector vJet;
						vJet.SetPxPyPzE(stuple_temp.jet_Px[j],stuple_temp.jet_Py[j],stuple_temp.jet_Pz[j],stuple_temp.jet_E[j]);
						stuple_temp.jet_Pt[j] = vJet.Pt();

						//following will work only if the EM obj hold the large fraction of energ from the cluster.
						// need to thinks of a solution to this.
						// this also effects when adding completely removed (fully macthed) objects
						stuple.jet_Index[stuple.jet_num] 	= stuple.pho_matchJetIndex[i];			//using the WRONG index? correct all these function!!
						stuple.jet_DetEta[stuple.jet_num] 	= stuple.pho_DetEta[i];
						stuple.jet_DetPhi[stuple.jet_num] 	= stuple.pho_DetPhi[i];
						stuple.jet_HadEm[stuple.jet_num] 	= stuple.pho_HadEm[i];
						stuple.jet_Emfr[stuple.jet_num] 		= 0;
						stuple.jet_Ntowers[stuple.jet_num] 	= 0;
						stuple.jet_Ntracks[stuple.jet_num] 	= 0;
						stuple.jet_SeedIPhi[stuple.jet_num]	= 0;
						stuple.jet_SeedIEta[stuple.jet_num] = 0;
						stuple.jet_EmTime[stuple.jet_num] = stuple.pho_EmTime[i];
						stuple.jet_SecVtxTag[stuple.jet_num] = -999999;
						stuple.jet_SecVtxppb[stuple.jet_num] = -999999;
						stuple.jet_SecVtxnpb[stuple.jet_num] = -999999;
						stuple.jet_SecVtxTrkmass[stuple.jet_num] = -999999;
			
						break;
					} // for

					if (!bFoundMatch) PhoNoMatch.push_back(i);

				} else PhoNoMatch.push_back(i);
	
			}	// for

			if (PhoNoMatch.size()>0) { 
				AddPhoToCentralJetList(stuple_temp, PhoNoMatch);
				bNeedSort = true;
			}
		
	} else {
		StdOut(__FILE__,__LINE__,3,"Number of used photons exceeds the number of photons in the list!");
		exit (1);
	}
	

	// add the electrons
	if (UsedEle.size() <= stuple.ele_num) {
		//add these phos
			std::vector<int> EleNoMatch;
			for (unsigned int i =0; i < stuple.ele_num; i++) {
				if (IsUsed(i, UsedEle)) continue;						//skip if this is the electron that was picked as THE electron
				if (stuple.ele_matchJetIndex[i] >=0) {
					bool bFoundMatch = false;
					for (unsigned int j=0; j < stuple.jet_num; j++) {
						if (stuple.ele_matchJetIndex[i] != stuple.jet_Index[j]) continue;
						bFoundMatch = true;
						bNeedSort = true;
						stuple_temp.jet_E[j]  = stuple.ele_E[i] + stuple.jet_E[j];
						stuple_temp.jet_Px[j] = stuple.ele_Px[i] + stuple.jet_Px[j];
						stuple_temp.jet_Py[j] = stuple.ele_Py[i] + stuple.jet_Py[j];
						stuple_temp.jet_Pz[j] = stuple.ele_Pz[i] + stuple.jet_Pz[j];
						stuple_temp.jet_Pt[j]  = stuple.ele_Etc[i] + stuple.jet_Pt[j];
						TLorentzVector vJet;
						vJet.SetPxPyPzE(stuple_temp.jet_Px[j],stuple_temp.jet_Py[j],stuple_temp.jet_Pz[j],stuple_temp.jet_E[j]);
						stuple_temp.jet_Pt[j] = vJet.Pt();

						//following will work only if the EM obj hold the large fraction of energ from the cluster.
						// need to thinks of a solution to this.
						// this also effects when adding completely removed (fully macthed) objects
						stuple.jet_Index[stuple.jet_num] 	= stuple.ele_matchJetIndex[i];
						stuple.jet_DetEta[stuple.jet_num] 	= stuple.ele_DetEta[i];
						stuple.jet_DetPhi[stuple.jet_num] 	= stuple.ele_DetPhi[i];
						stuple.jet_HadEm[stuple.jet_num] 	= stuple.ele_HadEm[i];
						stuple.jet_Emfr[stuple.jet_num] 		= 0;
						stuple.jet_Ntowers[stuple.jet_num] 	= 0;
						stuple.jet_Ntracks[stuple.jet_num] 	= 0;
						stuple.jet_SeedIPhi[stuple.jet_num]	= 0;
						stuple.jet_SeedIEta[stuple.jet_num] = 0;
						stuple.jet_EmTime[stuple.jet_num] = stuple.ele_EmTime[i];
						stuple.jet_SecVtxTag[stuple.jet_num] = -999999;
						stuple.jet_SecVtxppb[stuple.jet_num] = -999999;
						stuple.jet_SecVtxnpb[stuple.jet_num] = -999999;
						stuple.jet_SecVtxTrkmass[stuple.jet_num] = -999999;
			
						break;
					} // for

					if (!bFoundMatch) EleNoMatch.push_back(i);

				} else EleNoMatch.push_back(i);
	
			}	// for

			if (EleNoMatch.size()>0) { 
				AddEleToCentralJetList(stuple_temp, EleNoMatch);
				bNeedSort = true;
			}
		
	} else {
		StdOut(__FILE__,__LINE__,3,"Number of used electron exceeds the number of electrons in the list!");
		exit (1);
	}
	
		
	if (bNeedSort) {
		ReorderCentralJets(stuple_temp);
		stuple = stuple_temp;
	}

}

//-------------------------------------------------------------------
void NewJetList::AddPhoToCentralJetList(Stuple& stuple, std::vector<int> PhoNoMatch)
{

		//since I am adding a em obj, it may not have all the information (or variables)
		// defined for a jet. So I'll set those to zero for the em obj being added.
	
		for (unsigned int i=0; i < PhoNoMatch.size(); i++) {
			if (stuple.jet_num < stuple.Njet) {
				TLorentzVector vJet;
				vJet.SetPxPyPzE(stuple.jet_Px[stuple.jet_num],
									 stuple.jet_Py[stuple.jet_num],
									 stuple.jet_Pz[stuple.jet_num],
									 stuple.jet_E[stuple.jet_num]);
				stuple.jet_Pt[stuple.jet_num] = vJet.Pt();
				stuple.jet_E[stuple.jet_num] = stuple.pho_E[PhoNoMatch[i]]; 
				stuple.jet_Px[stuple.jet_num] = stuple.pho_Px[PhoNoMatch[i]]; 
				stuple.jet_Py[stuple.jet_num] = stuple.pho_Py[PhoNoMatch[i]]; 
				stuple.jet_Pz[stuple.jet_num] = stuple.pho_Pz[PhoNoMatch[i]]; 
//				std::cout << "new jet Pt= " << stuple.jet_Pt[stuple.jet_num] <<std::endl;

				stuple.jet_Index[stuple.jet_num] 	= stuple.pho_matchJetIndex[i];
				stuple.jet_DetEta[stuple.jet_num] 	= stuple.pho_DetEta[i];
				stuple.jet_DetPhi[stuple.jet_num] 	= stuple.pho_DetPhi[i];
				stuple.jet_HadEm[stuple.jet_num] 	= stuple.pho_HadEm[i];
				stuple.jet_Emfr[stuple.jet_num] 		= 0;
				stuple.jet_Ntowers[stuple.jet_num] 	= 0;
				stuple.jet_Ntracks[stuple.jet_num] 	= 0;
				stuple.jet_SeedIPhi[stuple.jet_num]	= 0;
				stuple.jet_SeedIEta[stuple.jet_num] = 0;
				stuple.jet_EmTime[stuple.jet_num] = stuple.pho_EmTime[i];
				stuple.jet_SecVtxTag[stuple.jet_num] = -999999;
				stuple.jet_SecVtxppb[stuple.jet_num] = -999999;
				stuple.jet_SecVtxnpb[stuple.jet_num] = -999999;
				stuple.jet_SecVtxTrkmass[stuple.jet_num] = -999999;
				stuple.jet_num++;
			} else {
				std::cout <<__FILE__ <<"::" << __LINE__ << ":: Maximum number of jets " << stuple.Njet << " reached while adding stuff to jet list!" << std::endl;
				exit (1);
			}
		}

}

//-------------------------------------------------------------------
void NewJetList::AddEleToCentralJetList(Stuple& stuple, std::vector<int> EleNoMatch)
{

		//since I am adding a em obj, it may not have all the information (or variables)
		// defined for a jet. So I'll set those to zero for the em obj being added.
	
		for (unsigned int i=0; i < EleNoMatch.size(); i++) {
			if (stuple.jet_num < stuple.Njet) {
				TLorentzVector vJet;
				vJet.SetPxPyPzE(stuple.jet_Px[stuple.jet_num],
									 stuple.jet_Py[stuple.jet_num],
									 stuple.jet_Pz[stuple.jet_num],
									 stuple.jet_E[stuple.jet_num]);
				stuple.jet_Pt[stuple.jet_num] = vJet.Pt();
				stuple.jet_E[stuple.jet_num] = stuple.ele_E[EleNoMatch[i]]; 
				stuple.jet_Px[stuple.jet_num] = stuple.ele_Px[EleNoMatch[i]]; 
				stuple.jet_Py[stuple.jet_num] = stuple.ele_Py[EleNoMatch[i]]; 
				stuple.jet_Pz[stuple.jet_num] = stuple.ele_Pz[EleNoMatch[i]]; 
//				std::cout << "new jet Pt= " << stuple.jet_Pt[stuple.jet_num] <<std::endl;

				stuple.jet_Index[stuple.jet_num] 	= stuple.ele_matchJetIndex[i];
				stuple.jet_DetEta[stuple.jet_num] 	= stuple.ele_DetEta[i];
				stuple.jet_DetPhi[stuple.jet_num] 	= stuple.ele_DetPhi[i];
				stuple.jet_HadEm[stuple.jet_num] 	= stuple.ele_HadEm[i];
				stuple.jet_Emfr[stuple.jet_num] 		= 0;
				stuple.jet_Ntowers[stuple.jet_num] 	= 0;
				stuple.jet_Ntracks[stuple.jet_num] 	= 0;
				stuple.jet_SeedIPhi[stuple.jet_num]	= 0;
				stuple.jet_SeedIEta[stuple.jet_num] = 0;
				stuple.jet_EmTime[stuple.jet_num] = stuple.ele_EmTime[i];
				stuple.jet_SecVtxTag[stuple.jet_num] = -999999;
				stuple.jet_SecVtxppb[stuple.jet_num] = -999999;
				stuple.jet_SecVtxnpb[stuple.jet_num] = -999999;
				stuple.jet_SecVtxTrkmass[stuple.jet_num] = -999999;

				stuple.jet_num++;
			} else {
				std::cout <<__FILE__ <<"::" << __LINE__ << ":: Maximum number of jets " << stuple.Njet << " reached while adding stuff to jet list!" << std::endl;
				exit (1);
			}
		}

}



//-------------------------------------------------------------------
void NewJetList::ReorderCentralJets(Stuple& stuple)
{
	int loopupto = stuple.jet_num - 1;

	while (loopupto !=0) {
		for (int i=0; i<loopupto;i++) {
			if (stuple.jet_Pt[i] < stuple.jet_Pt[i+1]) {
				Stuple stuple_temp;

				stuple_temp.jet_Index[i] 	= stuple.jet_Index[i];
				stuple_temp.jet_Pt[i] 		= stuple.jet_Pt[i];
				stuple_temp.jet_E[i] 		= stuple.jet_E[i];
				stuple_temp.jet_Px[i] 		= stuple.jet_Px[i];
				stuple_temp.jet_Py[i] 		= stuple.jet_Py[i];
				stuple_temp.jet_Pz[i] 		= stuple.jet_Pz[i];
				stuple_temp.jet_DetEta[i] 	= stuple.jet_DetEta[i];
				stuple_temp.jet_DetPhi[i] 	= stuple.jet_DetPhi[i];
				stuple_temp.jet_HadEm[i] 	= stuple.jet_HadEm[i];
				stuple_temp.jet_Emfr[i] 	= stuple.jet_Emfr[i];
				stuple_temp.jet_Ntowers[i] 	= stuple.jet_Ntowers[i];
				stuple_temp.jet_Ntracks[i] 	= stuple.jet_Ntracks[i];
				stuple_temp.jet_SeedIPhi[i] 	= stuple.jet_SeedIPhi[i];
				stuple_temp.jet_SeedIEta[i] 	= stuple.jet_SeedIEta[i];
				stuple_temp.jet_EmTime[i] = stuple.jet_EmTime[i];
				stuple_temp.jet_SecVtxTag[i] = stuple.jet_SecVtxTag[i];
				stuple_temp.jet_SecVtxppb[i] = stuple.jet_SecVtxppb[i];
				stuple_temp.jet_SecVtxnpb[i] = stuple.jet_SecVtxnpb[i];
				stuple_temp.jet_SecVtxTrkmass[i] = stuple.jet_SecVtxTrkmass[i];
		

				stuple.jet_Index[i] 	= stuple.jet_Index[i+1];
				stuple.jet_Pt[i] 		= stuple.jet_Pt[i+1];
				stuple.jet_E[i] 		= stuple.jet_E[i+1];
				stuple.jet_Px[i] 		= stuple.jet_Px[i+1];
				stuple.jet_Py[i] 		= stuple.jet_Py[i+1];
				stuple.jet_Pz[i] 		= stuple.jet_Pz[i+1];
				stuple.jet_DetEta[i] = stuple.jet_DetEta[i+1];
				stuple.jet_DetPhi[i] = stuple.jet_DetPhi[i+1];
				stuple.jet_HadEm[i] 	= stuple.jet_HadEm[i+1];
				stuple.jet_Emfr[i] 	= stuple.jet_Emfr[i+1];
				stuple.jet_Ntowers[i] 	= stuple.jet_Ntowers[i+1];
				stuple.jet_Ntracks[i] 	= stuple.jet_Ntracks[i+1];
				stuple.jet_SeedIPhi[i] 	= stuple.jet_SeedIPhi[i+1];
				stuple.jet_SeedIEta[i] 	= stuple.jet_SeedIEta[i+1];
				stuple.jet_EmTime[i]    = stuple.jet_EmTime[i+1];
				stuple.jet_SecVtxTag[i] = stuple.jet_SecVtxTag[i+1];
				stuple.jet_SecVtxppb[i] = stuple.jet_SecVtxppb[i+1];
				stuple.jet_SecVtxnpb[i] = stuple.jet_SecVtxnpb[i+1];
				stuple.jet_SecVtxTrkmass[i] = stuple.jet_SecVtxTrkmass[i+1];


				stuple.jet_Index[i+1] 	= stuple_temp.jet_Index[i];
				stuple.jet_Pt[i+1] 		= stuple_temp.jet_Pt[i];
				stuple.jet_E[i+1] 		= stuple_temp.jet_E[i];
				stuple.jet_Px[i+1] 		= stuple_temp.jet_Px[i];
				stuple.jet_Py[i+1] 		= stuple_temp.jet_Py[i];
				stuple.jet_Pz[i+1] 		= stuple_temp.jet_Pz[i];
				stuple.jet_DetEta[i+1] 	= stuple_temp.jet_DetEta[i];
				stuple.jet_DetPhi[i+1] 	= stuple_temp.jet_DetPhi[i];
				stuple.jet_HadEm[i+1] 	= stuple_temp.jet_HadEm[i];
				stuple.jet_Emfr[i+1] 	= stuple_temp.jet_Emfr[i];
				stuple.jet_Ntowers[i+1] = stuple_temp.jet_Ntowers[i];
				stuple.jet_Ntracks[i+1] = stuple_temp.jet_Ntracks[i];
				stuple.jet_SeedIPhi[i+1] = stuple_temp.jet_SeedIPhi[i];
				stuple.jet_SeedIEta[i+1] = stuple_temp.jet_SeedIEta[i];
				stuple.jet_EmTime[i] = stuple_temp.jet_EmTime[i];
				stuple.jet_SecVtxTag[i] = stuple_temp.jet_SecVtxTag[i];
				stuple.jet_SecVtxppb[i] = stuple_temp.jet_SecVtxppb[i];
				stuple.jet_SecVtxnpb[i] = stuple_temp.jet_SecVtxnpb[i];
				stuple.jet_SecVtxTrkmass[i] = stuple_temp.jet_SecVtxTrkmass[i];
					
			}
		} // for
		loopupto--;
	}
	
	for (unsigned int i=0 ; i < stuple.jet_num; i++) {
		if (i+1 < stuple.jet_num) assert(stuple.jet_Pt[i] >= stuple.jet_Pt[i+1]);
	}
}


// -------------------JES UP JET -----------------------------------------------------------------
void NewJetList::AddUnusedUp(Stuple& stuple, std::vector<int>UsedPho, std::vector<int>UsedEle)
{
	Stuple stuple_temp = stuple;
	
	bool bNeedSort = false;
	
	// add the un-used photons to jet list
	if (UsedPho.size() <= stuple.pho_up_num) {			// do not check for size>0, I am adding the un-used one in this loop.
		//add these phos
			std::vector<int> PhoNoMatch;
			for (unsigned int i =0; i < stuple.pho_up_num; i++) {
				if (IsUsed(i, UsedPho)) continue;						//check if this is the photon that was picked as THE photon
				if (stuple.pho_up_matchJetIndex[i] >=0) {
					bool bFoundMatch = false;
					for (unsigned int j=0; j < stuple.jet_up_num; j++) {
						if (stuple.pho_up_matchJetIndex[i] != stuple.jet_up_Index[j]) continue;
						bFoundMatch = true;
						bNeedSort = true;
						stuple_temp.jet_up_E[j]  = stuple.pho_up_E[i] + stuple.jet_up_E[j];
						stuple_temp.jet_up_Px[j] = stuple.pho_up_Px[i] + stuple.jet_up_Px[j];
						stuple_temp.jet_up_Py[j] = stuple.pho_up_Py[i] + stuple.jet_up_Py[j];
						stuple_temp.jet_up_Pz[j] = stuple.pho_up_Pz[i] + stuple.jet_up_Pz[j];
						stuple_temp.jet_up_Pt[j]  = stuple.pho_up_Etc[i] + stuple.jet_up_Pt[j];
						TLorentzVector vJet;
						vJet.SetPxPyPzE(stuple_temp.jet_up_Px[j],stuple_temp.jet_up_Py[j],stuple_temp.jet_up_Pz[j],stuple_temp.jet_up_E[j]);
						stuple_temp.jet_up_Pt[j] = vJet.Pt();

						//following will work only if the EM obj hold the large fraction of energ from the cluster.
						// need to thinks of a solution to this.
						// this also effects when adding completely removed (fully macthed) objects
						stuple.jet_up_Index[stuple.jet_up_num] 	= stuple.pho_up_matchJetIndex[i];
						stuple.jet_up_DetEta[stuple.jet_up_num] 	= stuple.pho_up_DetEta[i];
						stuple.jet_up_DetPhi[stuple.jet_up_num] 	= stuple.pho_up_DetPhi[i];
						stuple.jet_up_HadEm[stuple.jet_up_num] 	= stuple.pho_up_HadEm[i];
						stuple.jet_up_Emfr[stuple.jet_up_num] 		= 0;
						stuple.jet_up_Ntowers[stuple.jet_up_num] 	= 0;
						stuple.jet_up_Ntracks[stuple.jet_up_num] 	= 0;
						stuple.jet_up_SeedIPhi[stuple.jet_up_num]	= 0;
						stuple.jet_up_SeedIEta[stuple.jet_up_num] = 0;
			
						break;
					} // for

					if (!bFoundMatch) PhoNoMatch.push_back(i);

				} else PhoNoMatch.push_back(i);
	
			}	// for

			if (PhoNoMatch.size()>0) { 
				AddUpPhoToUpJetList(stuple_temp, PhoNoMatch);
				bNeedSort = true;
			}
		
	} else {
		StdOut(__FILE__,__LINE__,3,"Number of used photons exceeds the number of photons in the list!");
		exit (1);
	}
	

	// add the electrons
	if (UsedEle.size() <= stuple.ele_up_num) {
		//add these phos
			std::vector<int> EleNoMatch;
			for (unsigned int i =0; i < stuple.ele_up_num; i++) {
				if (IsUsed(i, UsedEle)) continue;						//skip if this is the electron that was picked as THE electron
				if (stuple.ele_up_matchJetIndex[i] >=0) {
					bool bFoundMatch = false;
					for (unsigned int j=0; j < stuple.jet_up_num; j++) {
						if (stuple.ele_up_matchJetIndex[i] != stuple.jet_up_Index[j]) continue;
						bFoundMatch = true;
						bNeedSort = true;
						stuple_temp.jet_up_E[j]  = stuple.ele_up_E[i] + stuple.jet_up_E[j];
						stuple_temp.jet_up_Px[j] = stuple.ele_up_Px[i] + stuple.jet_up_Px[j];
						stuple_temp.jet_up_Py[j] = stuple.ele_up_Py[i] + stuple.jet_up_Py[j];
						stuple_temp.jet_up_Pz[j] = stuple.ele_up_Pz[i] + stuple.jet_up_Pz[j];
						stuple_temp.jet_up_Pt[j]  = stuple.ele_up_Etc[i] + stuple.jet_up_Pt[j];
						TLorentzVector vJet;
						vJet.SetPxPyPzE(stuple_temp.jet_up_Px[j],stuple_temp.jet_up_Py[j],stuple_temp.jet_up_Pz[j],stuple_temp.jet_up_E[j]);
						stuple_temp.jet_up_Pt[j] = vJet.Pt();

						//following will work only if the EM obj hold the large fraction of energ from the cluster.
						// need to thinks of a solution to this.
						// this also effects when adding completely removed (fully macthed) objects
						stuple.jet_up_Index[stuple.jet_up_num] 	= stuple.ele_up_matchJetIndex[i];
						stuple.jet_up_DetEta[stuple.jet_up_num] 	= stuple.ele_up_DetEta[i];
						stuple.jet_up_DetPhi[stuple.jet_up_num] 	= stuple.ele_up_DetPhi[i];
						stuple.jet_up_HadEm[stuple.jet_up_num] 	= stuple.ele_up_HadEm[i];
						stuple.jet_up_Emfr[stuple.jet_up_num] 		= 0;
						stuple.jet_up_Ntowers[stuple.jet_up_num] 	= 0;
						stuple.jet_up_Ntracks[stuple.jet_up_num] 	= 0;
						stuple.jet_up_SeedIPhi[stuple.jet_up_num]	= 0;
						stuple.jet_up_SeedIEta[stuple.jet_up_num] = 0;
			
						break;
					} // for

					if (!bFoundMatch) EleNoMatch.push_back(i);

				} else EleNoMatch.push_back(i);
	
			}	// for

			if (EleNoMatch.size()>0) { 
				AddUpEleToUpJetList(stuple_temp, EleNoMatch);
				bNeedSort = true;
			}
		
	} else {
		StdOut(__FILE__,__LINE__,3,"Number of used electron exceeds the number of electrons in the list!");
		exit (1);
	}
	
		
	if (bNeedSort) {
		ReorderUpJets(stuple_temp);
		stuple = stuple_temp;
	}

}

//-------------------------------------------------------------------
void NewJetList::AddUpPhoToUpJetList(Stuple& stuple, std::vector<int> PhoNoMatch)
{

		//since I am adding a em obj, it may not have all the information (or variables)
		// defined for a jet. So I'll set those to zero for the em obj being added.
	
		for (unsigned int i=0; i < PhoNoMatch.size(); i++) {
			if (stuple.jet_up_num < stuple.Njet) {
				TLorentzVector vJet;
				vJet.SetPxPyPzE(stuple.jet_up_Px[stuple.jet_up_num],
									 stuple.jet_up_Py[stuple.jet_up_num],
									 stuple.jet_up_Pz[stuple.jet_up_num],
									 stuple.jet_up_E[stuple.jet_up_num]);
				stuple.jet_up_Pt[stuple.jet_up_num] = vJet.Pt();
				stuple.jet_up_E[stuple.jet_up_num] = stuple.pho_up_E[PhoNoMatch[i]]; 
				stuple.jet_up_Px[stuple.jet_up_num] = stuple.pho_up_Px[PhoNoMatch[i]]; 
				stuple.jet_up_Py[stuple.jet_up_num] = stuple.pho_up_Py[PhoNoMatch[i]]; 
				stuple.jet_up_Pz[stuple.jet_up_num] = stuple.pho_up_Pz[PhoNoMatch[i]]; 
//				std::cout << "new jet Pt= " << stuple.jet_up_Pt[stuple.jet_up_num] <<std::endl;

				stuple.jet_up_Index[stuple.jet_up_num] 	= stuple.pho_up_matchJetIndex[i];
				stuple.jet_up_DetEta[stuple.jet_up_num] 	= stuple.pho_up_DetEta[i];
				stuple.jet_up_DetPhi[stuple.jet_up_num] 	= stuple.pho_up_DetPhi[i];
				stuple.jet_up_HadEm[stuple.jet_up_num] 	= stuple.pho_up_HadEm[i];
				stuple.jet_up_Emfr[stuple.jet_up_num] 		= 0;
				stuple.jet_up_Ntowers[stuple.jet_up_num] 	= 0;
				stuple.jet_up_Ntracks[stuple.jet_up_num] 	= 0;
				stuple.jet_up_SeedIPhi[stuple.jet_up_num]	= 0;
				stuple.jet_up_SeedIEta[stuple.jet_up_num] = 0;
				stuple.jet_up_num++;
			} else {
				std::cout <<__FILE__ <<"::" << __LINE__ << ":: Maximum number of jets " << stuple.Njet << " reached while adding stuff to jet list!" << std::endl;
				exit (1);
			}
		}

}

//-------------------------------------------------------------------
void NewJetList::AddUpEleToUpJetList(Stuple& stuple, std::vector<int> EleNoMatch)
{

		//since I am adding a em obj, it may not have all the information (or variables)
		// defined for a jet. So I'll set those to zero for the em obj being added.
	
		for (unsigned int i=0; i < EleNoMatch.size(); i++) {
			if (stuple.jet_up_num < stuple.Njet) {
				TLorentzVector vJet;
				vJet.SetPxPyPzE(stuple.jet_up_Px[stuple.jet_up_num],
									 stuple.jet_up_Py[stuple.jet_up_num],
									 stuple.jet_up_Pz[stuple.jet_up_num],
									 stuple.jet_up_E[stuple.jet_up_num]);
				stuple.jet_up_Pt[stuple.jet_up_num] = vJet.Pt();
				stuple.jet_up_E[stuple.jet_up_num] = stuple.ele_up_E[EleNoMatch[i]]; 
				stuple.jet_up_Px[stuple.jet_up_num] = stuple.ele_up_Px[EleNoMatch[i]]; 
				stuple.jet_up_Py[stuple.jet_up_num] = stuple.ele_up_Py[EleNoMatch[i]]; 
				stuple.jet_up_Pz[stuple.jet_up_num] = stuple.ele_up_Pz[EleNoMatch[i]]; 
//				std::cout << "new jet Pt= " << stuple.jet_up_Pt[stuple.jet_up_num] <<std::endl;

				stuple.jet_up_Index[stuple.jet_up_num] 	= stuple.ele_up_matchJetIndex[i];
				stuple.jet_up_DetEta[stuple.jet_up_num] 	= stuple.ele_up_DetEta[i];
				stuple.jet_up_DetPhi[stuple.jet_up_num] 	= stuple.ele_up_DetPhi[i];
				stuple.jet_up_HadEm[stuple.jet_up_num] 	= stuple.ele_up_HadEm[i];
				stuple.jet_up_Emfr[stuple.jet_up_num] 		= 0;
				stuple.jet_up_Ntowers[stuple.jet_up_num] 	= 0;
				stuple.jet_up_Ntracks[stuple.jet_up_num] 	= 0;
				stuple.jet_up_SeedIPhi[stuple.jet_up_num]	= 0;
				stuple.jet_up_SeedIEta[stuple.jet_up_num] = 0;
				stuple.jet_up_num++;
			} else {
				std::cout <<__FILE__ <<"::" << __LINE__ << ":: Maximum number of jets " << stuple.Njet << " reached while adding stuff to jet list!" << std::endl;
				exit (1);
			}
		}

}



//-------------------------------------------------------------------
void NewJetList::ReorderUpJets(Stuple& stuple)
{
	int loopupto = stuple.jet_up_num - 1;

	while (loopupto !=0) {
		for (int i=0; i<loopupto;i++) {
			if (stuple.jet_up_Pt[i] < stuple.jet_up_Pt[i+1]) {
				Stuple stuple_temp;

				stuple_temp.jet_up_Index[i] 	= stuple.jet_up_Index[i];
				stuple_temp.jet_up_Pt[i] 		= stuple.jet_up_Pt[i];
				stuple_temp.jet_up_E[i] 		= stuple.jet_up_E[i];
				stuple_temp.jet_up_Px[i] 		= stuple.jet_up_Px[i];
				stuple_temp.jet_up_Py[i] 		= stuple.jet_up_Py[i];
				stuple_temp.jet_up_Pz[i] 		= stuple.jet_up_Pz[i];
				stuple_temp.jet_up_DetEta[i] 	= stuple.jet_up_DetEta[i];
				stuple_temp.jet_up_DetPhi[i] 	= stuple.jet_up_DetPhi[i];
				stuple_temp.jet_up_HadEm[i] 	= stuple.jet_up_HadEm[i];
				stuple_temp.jet_up_Emfr[i] 	= stuple.jet_up_Emfr[i];
				stuple_temp.jet_up_Ntowers[i] 	= stuple.jet_up_Ntowers[i];
				stuple_temp.jet_up_Ntracks[i] 	= stuple.jet_up_Ntracks[i];
				stuple_temp.jet_up_SeedIPhi[i] 	= stuple.jet_up_SeedIPhi[i];
				stuple_temp.jet_up_SeedIEta[i] 	= stuple.jet_up_SeedIEta[i];
		

				stuple.jet_up_Index[i] 	= stuple.jet_up_Index[i+1];
				stuple.jet_up_Pt[i] 		= stuple.jet_up_Pt[i+1];
				stuple.jet_up_E[i] 		= stuple.jet_up_E[i+1];
				stuple.jet_up_Px[i] 		= stuple.jet_up_Px[i+1];
				stuple.jet_up_Py[i] 		= stuple.jet_up_Py[i+1];
				stuple.jet_up_Pz[i] 		= stuple.jet_up_Pz[i+1];
				stuple.jet_up_DetEta[i] = stuple.jet_up_DetEta[i+1];
				stuple.jet_up_DetPhi[i] = stuple.jet_up_DetPhi[i+1];
				stuple.jet_up_HadEm[i] 	= stuple.jet_up_HadEm[i+1];
				stuple.jet_up_Emfr[i] 	= stuple.jet_up_Emfr[i+1];
				stuple.jet_up_Ntowers[i] 	= stuple.jet_up_Ntowers[i+1];
				stuple.jet_up_Ntracks[i] 	= stuple.jet_up_Ntracks[i+1];
				stuple.jet_up_SeedIPhi[i] 	= stuple.jet_up_SeedIPhi[i+1];
				stuple.jet_up_SeedIEta[i] 	= stuple.jet_up_SeedIEta[i+1];


				stuple.jet_up_Index[i+1] 	= stuple_temp.jet_up_Index[i];
				stuple.jet_up_Pt[i+1] 		= stuple_temp.jet_up_Pt[i];
				stuple.jet_up_E[i+1] 		= stuple_temp.jet_up_E[i];
				stuple.jet_up_Px[i+1] 		= stuple_temp.jet_up_Px[i];
				stuple.jet_up_Py[i+1] 		= stuple_temp.jet_up_Py[i];
				stuple.jet_up_Pz[i+1] 		= stuple_temp.jet_up_Pz[i];
				stuple.jet_up_DetEta[i+1] 	= stuple_temp.jet_up_DetEta[i];
				stuple.jet_up_DetPhi[i+1] 	= stuple_temp.jet_up_DetPhi[i];
				stuple.jet_up_HadEm[i+1] 	= stuple_temp.jet_up_HadEm[i];
				stuple.jet_up_Emfr[i+1] 	= stuple_temp.jet_up_Emfr[i];
				stuple.jet_up_Ntowers[i+1] = stuple_temp.jet_up_Ntowers[i];
				stuple.jet_up_Ntracks[i+1] = stuple_temp.jet_up_Ntracks[i];
				stuple.jet_up_SeedIPhi[i+1] = stuple_temp.jet_up_SeedIPhi[i];
				stuple.jet_up_SeedIEta[i+1] = stuple_temp.jet_up_SeedIEta[i];
					
			}
		} // for
		loopupto--;
	}
	
	for (unsigned int i=0 ; i < stuple.jet_up_num; i++) {
		if (i+1 < stuple.jet_up_num) assert(stuple.jet_up_Pt[i] >= stuple.jet_up_Pt[i+1]);
	}
}



// -------------------JES DOWN JET -----------------------------------------------------------------
void NewJetList::AddUnusedDown(Stuple& stuple, std::vector<int>UsedPho, std::vector<int>UsedEle)
{
	Stuple stuple_temp = stuple;
	
	bool bNeedSort = false;
	
	// add the un-used photons to jet list
	if (UsedPho.size() <= stuple.pho_down_num) {			// do not check for size>0, I am adding the un-used one in this loop.
		//add these phos
			std::vector<int> PhoNoMatch;
			for (unsigned int i =0; i < stuple.pho_down_num; i++) {
				if (IsUsed(i, UsedPho)) continue;						//check if this is the photon that was picked as THE photon
				if (stuple.pho_down_matchJetIndex[i] >=0) {
					bool bFoundMatch = false;
					for (unsigned int j=0; j < stuple.jet_down_num; j++) {
						if (stuple.pho_down_matchJetIndex[i] != stuple.jet_down_Index[j]) continue;
						bFoundMatch = true;
						bNeedSort = true;
						stuple_temp.jet_down_E[j]  = stuple.pho_down_E[i] + stuple.jet_down_E[j];
						stuple_temp.jet_down_Px[j] = stuple.pho_down_Px[i] + stuple.jet_down_Px[j];
						stuple_temp.jet_down_Py[j] = stuple.pho_down_Py[i] + stuple.jet_down_Py[j];
						stuple_temp.jet_down_Pz[j] = stuple.pho_down_Pz[i] + stuple.jet_down_Pz[j];
						stuple_temp.jet_down_Pt[j]  = stuple.pho_down_Etc[i] + stuple.jet_down_Pt[j];
						TLorentzVector vJet;
						vJet.SetPxPyPzE(stuple_temp.jet_down_Px[j],stuple_temp.jet_down_Py[j],stuple_temp.jet_down_Pz[j],stuple_temp.jet_down_E[j]);
						stuple_temp.jet_down_Pt[j] = vJet.Pt();

						//following will work only if the EM obj hold the large fraction of energ from the cluster.
						// need to thinks of a solution to this.
						// this also effects when adding completely removed (fully macthed) objects
						stuple.jet_down_Index[stuple.jet_down_num] 	= stuple.pho_down_matchJetIndex[i];
						stuple.jet_down_DetEta[stuple.jet_down_num] 	= stuple.pho_down_DetEta[i];
						stuple.jet_down_DetPhi[stuple.jet_down_num] 	= stuple.pho_down_DetPhi[i];
						stuple.jet_down_HadEm[stuple.jet_down_num] 	= stuple.pho_down_HadEm[i];
						stuple.jet_down_Emfr[stuple.jet_down_num] 		= 0;
						stuple.jet_down_Ntowers[stuple.jet_down_num] 	= 0;
						stuple.jet_down_Ntracks[stuple.jet_down_num] 	= 0;
						stuple.jet_down_SeedIPhi[stuple.jet_down_num]	= 0;
						stuple.jet_down_SeedIEta[stuple.jet_down_num] = 0;
			
						break;
					} // for

					if (!bFoundMatch) PhoNoMatch.push_back(i);

				} else PhoNoMatch.push_back(i);
	
			}	// for

			if (PhoNoMatch.size()>0) { 
				AddDownPhoToDownJetList(stuple_temp, PhoNoMatch);
				bNeedSort = true;
			}
		
	} else {
		StdOut(__FILE__,__LINE__,3,"Number of used photons exceeds the number of photons in the list!");
		exit (1);
	}
	

	// add the electrons
	if (UsedEle.size() <= stuple.ele_down_num) {
		//add these phos
			std::vector<int> EleNoMatch;
			for (unsigned int i =0; i < stuple.ele_down_num; i++) {
				if (IsUsed(i, UsedEle)) continue;						//skip if this is the electron that was picked as THE electron
				if (stuple.ele_down_matchJetIndex[i] >=0) {
					bool bFoundMatch = false;
					for (unsigned int j=0; j < stuple.jet_down_num; j++) {
						if (stuple.ele_down_matchJetIndex[i] != stuple.jet_down_Index[j]) continue;
						bFoundMatch = true;
						bNeedSort = true;
						stuple_temp.jet_down_E[j]  = stuple.ele_down_E[i] + stuple.jet_down_E[j];
						stuple_temp.jet_down_Px[j] = stuple.ele_down_Px[i] + stuple.jet_down_Px[j];
						stuple_temp.jet_down_Py[j] = stuple.ele_down_Py[i] + stuple.jet_down_Py[j];
						stuple_temp.jet_down_Pz[j] = stuple.ele_down_Pz[i] + stuple.jet_down_Pz[j];
						stuple_temp.jet_down_Pt[j]  = stuple.ele_down_Etc[i] + stuple.jet_down_Pt[j];
						TLorentzVector vJet;
						vJet.SetPxPyPzE(stuple_temp.jet_down_Px[j],stuple_temp.jet_down_Py[j],stuple_temp.jet_down_Pz[j],stuple_temp.jet_down_E[j]);
						stuple_temp.jet_down_Pt[j] = vJet.Pt();

						//following will work only if the EM obj hold the large fraction of energ from the cluster.
						// need to thinks of a solution to this.
						// this also effects when adding completely removed (fully macthed) objects
						stuple.jet_down_Index[stuple.jet_down_num] 	= stuple.ele_down_matchJetIndex[i];
						stuple.jet_down_DetEta[stuple.jet_down_num] 	= stuple.ele_down_DetEta[i];
						stuple.jet_down_DetPhi[stuple.jet_down_num] 	= stuple.ele_down_DetPhi[i];
						stuple.jet_down_HadEm[stuple.jet_down_num] 	= stuple.ele_down_HadEm[i];
						stuple.jet_down_Emfr[stuple.jet_down_num] 		= 0;
						stuple.jet_down_Ntowers[stuple.jet_down_num] 	= 0;
						stuple.jet_down_Ntracks[stuple.jet_down_num] 	= 0;
						stuple.jet_down_SeedIPhi[stuple.jet_down_num]	= 0;
						stuple.jet_down_SeedIEta[stuple.jet_down_num] = 0;
			
						break;
					} // for

					if (!bFoundMatch) EleNoMatch.push_back(i);

				} else EleNoMatch.push_back(i);
	
			}	// for

			if (EleNoMatch.size()>0) { 
				AddDownEleToDownJetList(stuple_temp, EleNoMatch);
				bNeedSort = true;
			}
		
	} else {
		StdOut(__FILE__,__LINE__,3,"Number of used electron exceeds the number of electrons in the list!");
		exit (1);
	}
	
		
	if (bNeedSort) {
		ReorderDownJets(stuple_temp);
		stuple = stuple_temp;
	}

}

//-------------------------------------------------------------------
void NewJetList::AddDownPhoToDownJetList(Stuple& stuple, std::vector<int> PhoNoMatch)
{

		//since I am adding a em obj, it may not have all the information (or variables)
		// defined for a jet. So I'll set those to zero for the em obj being added.
	
		for (unsigned int i=0; i < PhoNoMatch.size(); i++) {
			if (stuple.jet_down_num < stuple.Njet) {
				TLorentzVector vJet;
				vJet.SetPxPyPzE(stuple.jet_down_Px[stuple.jet_down_num],
									 stuple.jet_down_Py[stuple.jet_down_num],
									 stuple.jet_down_Pz[stuple.jet_down_num],
									 stuple.jet_down_E[stuple.jet_down_num]);
				stuple.jet_down_Pt[stuple.jet_down_num] = vJet.Pt();
				stuple.jet_down_E[stuple.jet_down_num] = stuple.pho_down_E[PhoNoMatch[i]]; 
				stuple.jet_down_Px[stuple.jet_down_num] = stuple.pho_down_Px[PhoNoMatch[i]]; 
				stuple.jet_down_Py[stuple.jet_down_num] = stuple.pho_down_Py[PhoNoMatch[i]]; 
				stuple.jet_down_Pz[stuple.jet_down_num] = stuple.pho_down_Pz[PhoNoMatch[i]]; 
//				std::cout << "new jet Pt= " << stuple.jet_down_Pt[stuple.jet_down_num] <<std::endl;

				stuple.jet_down_Index[stuple.jet_down_num] 	= stuple.pho_down_matchJetIndex[i];
				stuple.jet_down_DetEta[stuple.jet_down_num] 	= stuple.pho_down_DetEta[i];
				stuple.jet_down_DetPhi[stuple.jet_down_num] 	= stuple.pho_down_DetPhi[i];
				stuple.jet_down_HadEm[stuple.jet_down_num] 	= stuple.pho_down_HadEm[i];
				stuple.jet_down_Emfr[stuple.jet_down_num] 		= 0;
				stuple.jet_down_Ntowers[stuple.jet_down_num] 	= 0;
				stuple.jet_down_Ntracks[stuple.jet_down_num] 	= 0;
				stuple.jet_down_SeedIPhi[stuple.jet_down_num]	= 0;
				stuple.jet_down_SeedIEta[stuple.jet_down_num] = 0;
				stuple.jet_down_num++;
			} else {
				std::cout <<__FILE__ <<"::" << __LINE__ << ":: Maximum number of jets " << stuple.Njet << " reached while adding stuff to jet list!" << std::endl;
				exit (1);
			}
		}

}

//-------------------------------------------------------------------
void NewJetList::AddDownEleToDownJetList(Stuple& stuple, std::vector<int> EleNoMatch)
{

		//since I am adding a em obj, it may not have all the information (or variables)
		// defined for a jet. So I'll set those to zero for the em obj being added.
	
		for (unsigned int i=0; i < EleNoMatch.size(); i++) {
			if (stuple.jet_down_num < stuple.Njet) {
				TLorentzVector vJet;
				vJet.SetPxPyPzE(stuple.jet_down_Px[stuple.jet_down_num],
									 stuple.jet_down_Py[stuple.jet_down_num],
									 stuple.jet_down_Pz[stuple.jet_down_num],
									 stuple.jet_down_E[stuple.jet_down_num]);
				stuple.jet_down_Pt[stuple.jet_down_num] = vJet.Pt();
				stuple.jet_down_E[stuple.jet_down_num] = stuple.ele_down_E[EleNoMatch[i]]; 
				stuple.jet_down_Px[stuple.jet_down_num] = stuple.ele_down_Px[EleNoMatch[i]]; 
				stuple.jet_down_Py[stuple.jet_down_num] = stuple.ele_down_Py[EleNoMatch[i]]; 
				stuple.jet_down_Pz[stuple.jet_down_num] = stuple.ele_down_Pz[EleNoMatch[i]]; 
//				std::cout << "new jet Pt= " << stuple.jet_down_Pt[stuple.jet_down_num] <<std::endl;

				stuple.jet_down_Index[stuple.jet_down_num] 	= stuple.ele_down_matchJetIndex[i];
				stuple.jet_down_DetEta[stuple.jet_down_num] 	= stuple.ele_down_DetEta[i];
				stuple.jet_down_DetPhi[stuple.jet_down_num] 	= stuple.ele_down_DetPhi[i];
				stuple.jet_down_HadEm[stuple.jet_down_num] 	= stuple.ele_down_HadEm[i];
				stuple.jet_down_Emfr[stuple.jet_down_num] 		= 0;
				stuple.jet_down_Ntowers[stuple.jet_down_num] 	= 0;
				stuple.jet_down_Ntracks[stuple.jet_down_num] 	= 0;
				stuple.jet_down_SeedIPhi[stuple.jet_down_num]	= 0;
				stuple.jet_down_SeedIEta[stuple.jet_down_num] = 0;
				stuple.jet_down_num++;
			} else {
				std::cout <<__FILE__ <<"::" << __LINE__ << ":: Maximum number of jets " << stuple.Njet << " reached while adding stuff to jet list!" << std::endl;
				exit (1);
			}
		}

}



//-------------------------------------------------------------------
void NewJetList::ReorderDownJets(Stuple& stuple)
{
	int loopupto = stuple.jet_down_num - 1;

	while (loopupto !=0) {
		for (int i=0; i<loopupto;i++) {
			if (stuple.jet_down_Pt[i] < stuple.jet_down_Pt[i+1]) {
				Stuple stuple_temp;

				stuple_temp.jet_down_Index[i] 	= stuple.jet_down_Index[i];
				stuple_temp.jet_down_Pt[i] 		= stuple.jet_down_Pt[i];
				stuple_temp.jet_down_E[i] 		= stuple.jet_down_E[i];
				stuple_temp.jet_down_Px[i] 		= stuple.jet_down_Px[i];
				stuple_temp.jet_down_Py[i] 		= stuple.jet_down_Py[i];
				stuple_temp.jet_down_Pz[i] 		= stuple.jet_down_Pz[i];
				stuple_temp.jet_down_DetEta[i] 	= stuple.jet_down_DetEta[i];
				stuple_temp.jet_down_DetPhi[i] 	= stuple.jet_down_DetPhi[i];
				stuple_temp.jet_down_HadEm[i] 	= stuple.jet_down_HadEm[i];
				stuple_temp.jet_down_Emfr[i] 	= stuple.jet_down_Emfr[i];
				stuple_temp.jet_down_Ntowers[i] 	= stuple.jet_down_Ntowers[i];
				stuple_temp.jet_down_Ntracks[i] 	= stuple.jet_down_Ntracks[i];
				stuple_temp.jet_down_SeedIPhi[i] 	= stuple.jet_down_SeedIPhi[i];
				stuple_temp.jet_down_SeedIEta[i] 	= stuple.jet_down_SeedIEta[i];
		

				stuple.jet_down_Index[i] 	= stuple.jet_down_Index[i+1];
				stuple.jet_down_Pt[i] 		= stuple.jet_down_Pt[i+1];
				stuple.jet_down_E[i] 		= stuple.jet_down_E[i+1];
				stuple.jet_down_Px[i] 		= stuple.jet_down_Px[i+1];
				stuple.jet_down_Py[i] 		= stuple.jet_down_Py[i+1];
				stuple.jet_down_Pz[i] 		= stuple.jet_down_Pz[i+1];
				stuple.jet_down_DetEta[i] = stuple.jet_down_DetEta[i+1];
				stuple.jet_down_DetPhi[i] = stuple.jet_down_DetPhi[i+1];
				stuple.jet_down_HadEm[i] 	= stuple.jet_down_HadEm[i+1];
				stuple.jet_down_Emfr[i] 	= stuple.jet_down_Emfr[i+1];
				stuple.jet_down_Ntowers[i] 	= stuple.jet_down_Ntowers[i+1];
				stuple.jet_down_Ntracks[i] 	= stuple.jet_down_Ntracks[i+1];
				stuple.jet_down_SeedIPhi[i] 	= stuple.jet_down_SeedIPhi[i+1];
				stuple.jet_down_SeedIEta[i] 	= stuple.jet_down_SeedIEta[i+1];


				stuple.jet_down_Index[i+1] 	= stuple_temp.jet_down_Index[i];
				stuple.jet_down_Pt[i+1] 		= stuple_temp.jet_down_Pt[i];
				stuple.jet_down_E[i+1] 		= stuple_temp.jet_down_E[i];
				stuple.jet_down_Px[i+1] 		= stuple_temp.jet_down_Px[i];
				stuple.jet_down_Py[i+1] 		= stuple_temp.jet_down_Py[i];
				stuple.jet_down_Pz[i+1] 		= stuple_temp.jet_down_Pz[i];
				stuple.jet_down_DetEta[i+1] 	= stuple_temp.jet_down_DetEta[i];
				stuple.jet_down_DetPhi[i+1] 	= stuple_temp.jet_down_DetPhi[i];
				stuple.jet_down_HadEm[i+1] 	= stuple_temp.jet_down_HadEm[i];
				stuple.jet_down_Emfr[i+1] 	= stuple_temp.jet_down_Emfr[i];
				stuple.jet_down_Ntowers[i+1] = stuple_temp.jet_down_Ntowers[i];
				stuple.jet_down_Ntracks[i+1] = stuple_temp.jet_down_Ntracks[i];
				stuple.jet_down_SeedIPhi[i+1] = stuple_temp.jet_down_SeedIPhi[i];
				stuple.jet_down_SeedIEta[i+1] = stuple_temp.jet_down_SeedIEta[i];
					
			}
		} // for
		loopupto--;
	}
	
	for (unsigned int i=0 ; i < stuple.jet_down_num; i++) {
		if (i+1 < stuple.jet_down_num) assert(stuple.jet_down_Pt[i] >= stuple.jet_down_Pt[i+1]);
	}
}



void NewJetList::AddUnused(PhotonList& phos, ElectronList& eles, JetList& jets,
							const std::vector<int>UsedPho, const std::vector<int>UsedEle,
							const float fMinJetEt, const float fMaxJetEta)
{
	bool bNeedSort = false;

	PhotonList tphos = phos;
	ElectronList teles = eles;
	JetList tjets = jets;
	
	assert ( tphos.pho_num == tphos.pho_Etc.size() && "NewJetList::AddUnused:: pho looping over numbers differ!");
	assert ( teles.ele_num == teles.ele_Etc.size() && "NewJetList::AddUnused:: ele looping over numbers differ!");
	assert ( tjets.jet_num == tjets.jet_Pt.size() && "NewJetList::AddUnused:: jet looping over numbers differ!");

	std::vector<unsigned int> PhoNoMatch;
	for (unsigned int i =0; i < tphos.pho_num; i++)
	{
		if (IsUsed(i, UsedPho)) continue;						//check if this is the photon that was picked as THE photon

		if (tphos.pho_matchJetIndex[i] >=0)
		{
			bool bFoundMatch = false;
			for (unsigned int j=0; j < jets.jet_num; j++)
			{
				if (tphos.pho_matchJetIndex[i] != jets.jet_Index[j]) continue;
				bFoundMatch = true;
				bNeedSort = true;

				//std::cout << __FILE__ << ":" << __LINE__ << " index pho, jet = " << tphos.pho_matchJetIndex[i] << ", " << jets.jet_Index[j] << std::endl;
				TLorentzVector pvec(phos.pho_Px[i], phos.pho_Py[i], phos.pho_Pz[i], phos.pho_E[i]);
				TLorentzVector jvec(jets.jet_Px[j], jets.jet_Py[j], jets.jet_Pz[j], jets.jet_E[j]);
				TLorentzVector sum = pvec + jvec;
				
				tjets.jet_E[j]  = phos.pho_E[i]  + jets.jet_E[j];			// i should remove this em objs after they are added to the jet!
				tjets.jet_Px[j] = phos.pho_Px[i] + jets.jet_Px[j];
				tjets.jet_Py[j] = phos.pho_Py[i] + jets.jet_Py[j];
				tjets.jet_Pz[j] = phos.pho_Pz[i] + jets.jet_Pz[j];
				TLorentzVector vJet;
				vJet.SetPxPyPzE(tjets.jet_Px[j], tjets.jet_Py[j], tjets.jet_Pz[j], tjets.jet_E[j]);
				tjets.jet_Pt[j] = vJet.Pt();

				//std::cout << "old/new E = " << tjets.jet_E[j] << ", " << sum.E() << std::endl;
				//std::cout << "old/new Px= " << tjets.jet_Px[j] << ", " << sum.Px() << std::endl;
				//std::cout << "old/new Py= " << tjets.jet_Py[j] << ", " << sum.Py() << std::endl;
				//std::cout << "old/new Pz= " << tjets.jet_Pz[j] << ", " << sum.Pz() << std::endl;
				//std::cout << "old/new Pt= " << tjets.jet_Pt[j] << ", " << sum.Pt() << std::endl;

				
				//following will work only if the EM obj hold the large fraction of energy from the cluster.
				// need to thinks of a solution to this.
				// this also effects when adding completely removed (fully macthed) objects
				//stuple.jet_Index[stuple.jet_num] 	= stuple.pho_matchJetIndex[i]; // this should be there, in this case
				tjets.jet_DetEta[j] 	= phos.pho_DetEta[i];	
				tjets.jet_DetPhi[j] 	= phos.pho_DetPhi[i];
				tjets.jet_HadEm[j] 	= phos.pho_HadEm[i];
				tjets.jet_Emfr[j] 		= 1.0/(1.0+phos.pho_HadEm[i]);
				
				/*These varaiables do not change when em obj is removed from a jet. so keep them as it is. -05-06-2010
				
				tjets.jet_Ntowers[j] 	= 0;
				tjets.jet_Ntracks[j] 	= 0;
				tjets.jet_SeedIPhi[j]	= 0;
				tjets.jet_SeedIEta[j] = 0;
				*/
				bNeedSort = true;

				break;
			} // for

			if (!bFoundMatch) PhoNoMatch.push_back(i);

		} else PhoNoMatch.push_back(i);

	}	// for

	// now add the em objs that were completely removed from jet list.
	if (PhoNoMatch.size()>0)
	{
		for (unsigned int i=0; i < PhoNoMatch.size(); i++)
		{
			if (phos.pho_Etc[PhoNoMatch[i]] > fMinJetEt && 
				fabs(phos.pho_DetEta.at(i)) < fMaxJetEta)			//I need jets above 15GeV. This is not required as all my EM objs are Et>30GeV
			{
				tjets.jet_num++;							// number of jets so we know how much to read from arrays
				tjets.jet_NJet15++;
				tjets.jet_Index.push_back(phos.pho_matchJetIndex[PhoNoMatch.at(i)]);				//original index in jet block
				TLorentzVector vJet;
				vJet.SetPxPyPzE(phos.pho_Px[PhoNoMatch.at(i)], phos.pho_Py[PhoNoMatch.at(i)], phos.pho_Pz[PhoNoMatch.at(i)], phos.pho_E[PhoNoMatch.at(i)]);
				tjets.jet_Pt.push_back(vJet.Pt());
				tjets.jet_E.push_back(phos.pho_E[PhoNoMatch.at(i)]);
				tjets.jet_Px.push_back(phos.pho_Px[PhoNoMatch.at(i)]);
				tjets.jet_Py.push_back(phos.pho_Py[PhoNoMatch.at(i)]);
				tjets.jet_Pz.push_back(phos.pho_Pz[PhoNoMatch.at(i)]);
				tjets.jet_DetEta.push_back(phos.pho_DetEta[PhoNoMatch.at(i)]);
				tjets.jet_DetPhi.push_back(phos.pho_DetPhi[PhoNoMatch.at(i)]);
				float hm=phos.pho_HadEm[PhoNoMatch.at(i)];
				tjets.jet_HadEm.push_back(hm);
				tjets.jet_Emfr.push_back(1/(hm+1.0));
				tjets.jet_Ntowers.push_back(-999999);
				tjets.jet_Ntracks.push_back(-999999);
				tjets.jet_SeedIPhi.push_back(-999999);
				tjets.jet_SeedIEta.push_back(-999999);
				tjets.jet_EmTime.push_back(phos.pho_EmTime[PhoNoMatch.at(i)]);
				tjets.jet_SecVtxTag.push_back(-999999);
				tjets.jet_SecVtxppb.push_back(-999999);
				tjets.jet_SecVtxnpb.push_back(-999999);
				tjets.jet_SecVtxTrkmass.push_back(-999999);

			}
		}

		bNeedSort = true;

		// now remove the photons added to jet list from photon list
		PhotonList newpholist;

		for (unsigned int i=0; i < PhoNoMatch.size(); i++)
		{
			for (unsigned int j=0; j < phos.pho_num; j++)
			{
				if (j == PhoNoMatch[i]) continue;			//true=this has been added to the jet list.

				newpholist.pho_num++;
				if (phos.pho_TightId[j] == 0) newpholist.pho_Ntight++; 
				if (phos.pho_LooseId[j] == 0) newpholist.pho_Nloose++;

				newpholist.pho_Index.push_back(phos.pho_Index[j]);
				newpholist.pho_PhoBlockIndex.push_back(phos.pho_PhoBlockIndex[j]);
				newpholist.pho_Etc.push_back(phos.pho_Etc[j]);
				newpholist.pho_E.push_back(phos.pho_E[j]);
				newpholist.pho_Px.push_back(phos.pho_Px[j]);
				newpholist.pho_Py.push_back(phos.pho_Py[j]);
				newpholist.pho_Pz.push_back(phos.pho_Pz[j]);
				newpholist.pho_Detector.push_back(phos.pho_Detector[j]);
				newpholist.pho_DetEta.push_back(phos.pho_DetEta[j]);
				newpholist.pho_DetPhi.push_back(phos.pho_DetPhi[j]);
				newpholist.pho_XCes.push_back(phos.pho_XCes[j]);
				newpholist.pho_ZCes.push_back(phos.pho_ZCes[j]);
				newpholist.pho_HadEm.push_back(phos.pho_HadEm[j]);
				newpholist.pho_Chi2Mean.push_back(phos.pho_Chi2Mean[j]);
				newpholist.pho_N3d.push_back(phos.pho_N3d[j]);
				newpholist.pho_Iso4.push_back(phos.pho_Iso4[j]);
				newpholist.pho_TrkPt.push_back(phos.pho_TrkPt[j]);	
				newpholist.pho_TrkIso.push_back(phos.pho_TrkIso[j]);
				newpholist.pho_CesWireE2.push_back(phos.pho_CesWireE2[j]);
				newpholist.pho_CesStripE2.push_back(phos.pho_CesStripE2[j]);
				newpholist.pho_PhiWedge.push_back(phos.pho_PhiWedge[j]);
				newpholist.pho_NMuonStubs.push_back(phos.pho_NMuonStubs[j]);
				newpholist.pho_EmTime.push_back(phos.pho_EmTime[j]);
				newpholist.pho_TightId.push_back(phos.pho_TightId[j]);
				newpholist.pho_LooseId.push_back(phos.pho_LooseId[j]);
				newpholist.pho_PhoenixId.push_back(phos.pho_PhoenixId[j]);
				newpholist.pho_Halo_seedWedge.push_back(phos.pho_Halo_seedWedge[j]);
				newpholist.pho_Halo_eastNhad.push_back(phos.pho_Halo_eastNhad[j]);
				newpholist.pho_Halo_westNhad.push_back(phos.pho_Halo_westNhad[j]);
				newpholist.pho_matchJetIndex.push_back(phos.pho_matchJetIndex[j]);

			}
		}

		//phos.Clear();			//don't do this yet. I should change whole thing to use the original em ind.
		//phos = newpholist;		// so when I remove and compact, the array index is change and things messes up!
	}	// if PhoNoMatch.size>0




	std::vector<unsigned int> EleNoMatch;
	for (unsigned int i =0; i < teles.ele_num; i++)
	{
		if (IsUsed(i, UsedEle)) continue;						//check if this is the photon that was picked as THE photon
		if (teles.ele_matchJetIndex[i] >=0)
		{
			bool bFoundMatch = false;
			for (unsigned int j=0; j < jets.jet_num; j++)
			{
				if (teles.ele_matchJetIndex[i] != jets.jet_Index[j]) continue;
				bFoundMatch = true;
				bNeedSort = true;
				tjets.jet_E[j]  = eles.ele_E[i] + jets.jet_E[j];
				tjets.jet_Px[j] = eles.ele_Px[i] + jets.jet_Px[j];
				tjets.jet_Py[j] = eles.ele_Py[i] + jets.jet_Py[j];
				tjets.jet_Pz[j] = eles.ele_Pz[i] + jets.jet_Pz[j];
				TLorentzVector vJet;
				vJet.SetPxPyPzE(jets.jet_Px[j], jets.jet_Py[j], jets.jet_Pz[j], jets.jet_E[j]);
				jets.jet_Pt[j] = vJet.Pt();

				//following will work only if the EM obj hold the large fraction of energy from the cluster.
				// need to thinks of a solution to this.
				// this also effects when adding completely removed (fully macthed) objects
				//stuple.jet_Index[stuple.jet_num] 	= stuple.pho_matchJetIndex[i]; // this should be there, in this case
				tjets.jet_DetEta[j] 	= eles.ele_DetEta[i];	
				tjets.jet_DetPhi[j] 	= eles.ele_DetPhi[i];
				tjets.jet_HadEm[j] 	= eles.ele_HadEm[i];
				tjets.jet_Emfr[j] 		= 1.0/(1.0+eles.ele_HadEm[i]);
				/*These varaiables do not change when em obj is removed from a jet. so keep them as it is. -05-06-2010
				tjets.jet_Ntowers[j] 	= 0;
				tjets.jet_Ntracks[j] 	= 0;
				tjets.jet_SeedIPhi[j]	= 0;
				tjets.jet_SeedIEta[j] = 0;
				*/
				bNeedSort = true;

				break;
			} // for

			if (!bFoundMatch) EleNoMatch.push_back(i);

		} else EleNoMatch.push_back(i);

	}	// for

	// now add the em objs that were completely removed from jet list.
	if (EleNoMatch.size()>0)
	{
		for (unsigned int i=0; i < EleNoMatch.size(); i++)
		{
			if (eles.ele_Etc[EleNoMatch[i]] > 15.0)			//I need jets above 15GeV. This is not required as all my EM objs are Et>30GeV
			{
				tjets.jet_num++;							// number of jets so we know how much to read from arrays
				tjets.jet_NJet15++;
				tjets.jet_Index.push_back(eles.ele_matchJetIndex[EleNoMatch.at(i)]);				//original index in jet block
				TLorentzVector vJet;
				vJet.SetPxPyPzE(eles.ele_Px[EleNoMatch.at(i)], eles.ele_Py[EleNoMatch.at(i)], eles.ele_Pz[EleNoMatch.at(i)], eles.ele_E[EleNoMatch.at(i)]);
				tjets.jet_Pt.push_back(vJet.Pt());
				tjets.jet_E.push_back(eles.ele_E[EleNoMatch.at(i)]);
				tjets.jet_Px.push_back(eles.ele_Px[EleNoMatch.at(i)]);
				tjets.jet_Py.push_back(eles.ele_Py[EleNoMatch.at(i)]);
				tjets.jet_Pz.push_back(eles.ele_Pz[EleNoMatch.at(i)]);
				tjets.jet_DetEta.push_back(eles.ele_DetEta[EleNoMatch.at(i)]);
				tjets.jet_DetPhi.push_back(eles.ele_DetPhi[EleNoMatch.at(i)]);
				float hm=eles.ele_HadEm[EleNoMatch.at(i)];
				tjets.jet_HadEm.push_back(hm);
				tjets.jet_Emfr.push_back(1/(hm+1.0));
				tjets.jet_Ntowers.push_back(-999999);
				tjets.jet_Ntracks.push_back(-999999);
				tjets.jet_SeedIPhi.push_back(-999999);
				tjets.jet_SeedIEta.push_back(-999999);
				tjets.jet_EmTime.push_back(eles.ele_EmTime[EleNoMatch.at(i)]);
				tjets.jet_SecVtxTag.push_back(-999999);
				tjets.jet_SecVtxppb.push_back(-999999);
				tjets.jet_SecVtxnpb.push_back(-999999);
				tjets.jet_SecVtxTrkmass.push_back(-999999);
			}
		}

		bNeedSort = true;

		// now remove the electrons added to jet list from photon list
		ElectronList newelelist;

		for (unsigned int i=0; i < EleNoMatch.size(); i++)
		{
			for (unsigned int j=0; j < eles.ele_num; j++)
			{
				if (j != EleNoMatch[i]) continue;			//true=this has not been added to the jet list.

				newelelist.ele_num++;
				if (eles.ele_TightId[j] == 0) newelelist.ele_Ntight++; 
				if (eles.ele_LooseId[j] == 0) newelelist.ele_Nloose++;
				newelelist.ele_Index.push_back(teles.ele_Index[j]);
				newelelist.ele_PhoBlockIndex.push_back(teles.ele_PhoBlockIndex[j]);
				newelelist.ele_EleBlockIndex.push_back(teles.ele_EleBlockIndex[j]);
				newelelist.ele_Etc.push_back(teles.ele_Etc[j]);
				newelelist.ele_E.push_back(teles.ele_E[j]);
				newelelist.ele_Px.push_back(teles.ele_Px[j]);
				newelelist.ele_Py.push_back(teles.ele_Py[j]);
				newelelist.ele_Pz.push_back(teles.ele_Pz[j]);
				newelelist.ele_Detector.push_back(teles.ele_Detector[j]);
				newelelist.ele_DetEta.push_back(teles.ele_DetEta[j]);
				newelelist.ele_DetPhi.push_back(teles.ele_DetPhi[j]);
				newelelist.ele_XCes.push_back(teles.ele_XCes[j]);
				newelelist.ele_ZCes.push_back(teles.ele_ZCes[j]);
				newelelist.ele_HadEm.push_back(teles.ele_HadEm[j]);
				newelelist.ele_Chi2Mean.push_back(teles.ele_Chi2Mean[j]);
				newelelist.ele_N3d.push_back(teles.ele_N3d[j]);
				newelelist.ele_Iso4.push_back(teles.ele_Iso4[j]);
				newelelist.ele_TrkIso.push_back(teles.ele_TrkIso[j]);		
				newelelist.ele_CesWireE2.push_back(teles.ele_CesWireE2[j]);
				newelelist.ele_CesStripE2.push_back(teles.ele_CesStripE2[j]);
				newelelist.ele_PhiWedge.push_back(teles.ele_PhiWedge[j]);
				newelelist.ele_NMuonStubs.push_back(teles.ele_NMuonStubs[j]);
				newelelist.ele_EmTime.push_back(teles.ele_EmTime[j]);
				newelelist.ele_PhoenixId.push_back(teles.ele_PhoenixId[j]);
				newelelist.ele_Halo_seedWedge.push_back(teles.ele_Halo_seedWedge[j]);
				newelelist.ele_Halo_eastNhad.push_back(teles.ele_Halo_eastNhad[j]);
				newelelist.ele_Halo_westNhad.push_back(teles.ele_Halo_westNhad[j]);
				newelelist.ele_matchJetIndex.push_back(teles.ele_matchJetIndex[j]);
				newelelist.ele_Ntracks.push_back(teles.ele_Ntracks[j]);	
				newelelist.ele_Emfr.push_back(teles.ele_Emfr[j]);
				newelelist.ele_EoverP.push_back(teles.ele_EoverP[j]);
				newelelist.ele_TrackPt.push_back(teles.ele_TrackPt[j]);
				newelelist.ele_TrackBcPt.push_back(teles.ele_TrackBcPt[j]);
				newelelist.ele_TrackPhi.push_back(teles.ele_TrackPhi[j]);
				newelelist.ele_Nssl.push_back(teles.ele_Nssl[j]);
				newelelist.ele_Nasl.push_back(teles.ele_Nasl[j]);
				newelelist.ele_TightId.push_back(teles.ele_TightId[j]);
				newelelist.ele_LooseId.push_back(teles.ele_LooseId[j]);
				newelelist.ele_ConversionId.push_back(teles.ele_ConversionId[j]);
			}
		}

		//eles.Clear();			//don't do this yet. I should change whole thing to use the original em ind.
		//eles = newelelist;		// so when I remove and compact, the array index is change and things messes up!
	}	// if EleNoMatch.size>0


	if (bNeedSort)
	{
		ReorderJets(tjets);
		jets.Clear();
		jets = tjets;
	}
}

void NewJetList::ReorderJets(JetList& jets)
{
	//when there are no jets to begin with
	if (jets.jet_num == 0) return;

	int loopupto = jets.jet_num - 1;

	while (loopupto !=0) {
		for (int i=0; i<loopupto;i++) {
			if (jets.jet_Pt[i] < jets.jet_Pt[i+1]) {

				int jet_Index 		= jets.jet_Index[i];
				float jet_Pt 		= jets.jet_Pt[i];
				float jet_E 		= jets.jet_E[i];
				float jet_Px 		= jets.jet_Px[i];
				float jet_Py 		= jets.jet_Py[i];
				float jet_Pz 		= jets.jet_Pz[i];
				float jet_DetEta 	= jets.jet_DetEta[i];
				float jet_DetPhi 	= jets.jet_DetPhi[i];
				float jet_HadEm 	= jets.jet_HadEm[i];
				float jet_Emfr 	= jets.jet_Emfr[i];
				int jet_Ntowers 	= jets.jet_Ntowers[i];
				int jet_Ntracks 	= jets.jet_Ntracks[i];
				int jet_SeedIPhi 	= jets.jet_SeedIPhi[i];
				int jet_SeedIEta 	= jets.jet_SeedIEta[i];
				float jet_EmTime = jets.jet_EmTime[i];
				int   jet_SecVtxTag = jets.jet_SecVtxTag[i];
				float jet_SecVtxppb = jets.jet_SecVtxppb[i];
				float jet_SecVtxnpb = jets.jet_SecVtxnpb[i];
				float jet_SecVtxTrkmass = jets.jet_SecVtxTrkmass[i];
		
				jets.jet_Index[i] 	= jets.jet_Index[i+1];
				jets.jet_Pt[i] 		= jets.jet_Pt[i+1];
				jets.jet_E[i] 			= jets.jet_E[i+1];
				jets.jet_Px[i] 		= jets.jet_Px[i+1];
				jets.jet_Py[i] 		= jets.jet_Py[i+1];
				jets.jet_Pz[i] 		= jets.jet_Pz[i+1];
				jets.jet_DetEta[i] 	= jets.jet_DetEta[i+1];
				jets.jet_DetPhi[i] 	= jets.jet_DetPhi[i+1];
				jets.jet_HadEm[i] 	= jets.jet_HadEm[i+1];
				jets.jet_Emfr[i] 		= jets.jet_Emfr[i+1];
				jets.jet_Ntowers[i] 	= jets.jet_Ntowers[i+1];
				jets.jet_Ntracks[i] 	= jets.jet_Ntracks[i+1];
				jets.jet_SeedIPhi[i] = jets.jet_SeedIPhi[i+1];
				jets.jet_SeedIEta[i] = jets.jet_SeedIEta[i+1];
				jets.jet_EmTime[i]   = jets.jet_EmTime[i+1];
				jets.jet_SecVtxTag[i] = jets.jet_SecVtxTag[i+1];
				jets.jet_SecVtxppb[i] = jets.jet_SecVtxppb[i+1];
				jets.jet_SecVtxnpb[i] = jets.jet_SecVtxnpb[i+1];
				jets.jet_SecVtxTrkmass[i] = jets.jet_SecVtxTrkmass[i+1];

				jets.jet_Index[i+1] 	= jet_Index;
				jets.jet_Pt[i+1] 		= jet_Pt;
				jets.jet_E[i+1] 		= jet_E;
				jets.jet_Px[i+1] 		= jet_Px;
				jets.jet_Py[i+1] 		= jet_Py;
				jets.jet_Pz[i+1] 		= jet_Pz;
				jets.jet_DetEta[i+1] = jet_DetEta;
				jets.jet_DetPhi[i+1] = jet_DetPhi;
				jets.jet_HadEm[i+1] 	= jet_HadEm;
				jets.jet_Emfr[i+1] 	= jet_Emfr;
				jets.jet_Ntowers[i+1] = jet_Ntowers;
				jets.jet_Ntracks[i+1] = jet_Ntracks;
				jets.jet_SeedIPhi[i+1] = jet_SeedIPhi;
				jets.jet_SeedIEta[i+1] = jet_SeedIEta;
				jets.jet_EmTime[i+1]    = jet_EmTime;
				jets.jet_SecVtxTag[i+1] = jet_SecVtxTag;
				jets.jet_SecVtxppb[i+1] = jet_SecVtxppb;
				jets.jet_SecVtxnpb[i+1] = jet_SecVtxnpb;
				jets.jet_SecVtxTrkmass[i+1] = jet_SecVtxTrkmass;
			}
		} // for
		loopupto--;
	}
	
	for (unsigned int i=0 ; i < jets.jet_num; i++) {
		if (i+1 < jets.jet_num) assert(jets.jet_Pt[i] >= jets.jet_Pt[i+1]);
	}

	
}

