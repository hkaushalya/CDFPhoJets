//this is to annotate the root files created in jobs
//with a READ ME canvas
#include <iostream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TDatime.h"

//void Annotate()
void Format(const std::string readme, std::vector<std::string>& lines)
{
	//break long lines into 40 column text
	const int max_len = 40;
	const int strlen = readme.length();
	if (strlen<=max_len) lines.push_back(readme);
	else 
	{
		bool iterate = true;
		size_t i = 0;
		while (iterate)
		{
			int j = max_len;
			if (i+max_len>strlen) { j = strlen - i+1; iterate = false; }
			const std::string substr = readme.substr(i,j);
			std::cout << substr << std::endl;
			lines.push_back(substr);
			i +=max_len;
		}
		TDatime date;
		std::string sdate(date.AsString());
		lines.push_back(sdate);

	}

}

void Annotate(const std::string file="test.root", const std::string readme="test read me")
{
	if (file.length() < 1 || readme.length()<1)
	{
		std::cout << "file or 'readme' text not found! returning.." << std::endl;
		return;
	}
	TFile f(file.c_str(),"UPDATE");
	if (f.IsZombie())
	{
		std::cout << "File not found! returning.." << std::endl;
		return;
	}
	
	std::vector<std::string> lines;
	Format(readme, lines);
	const std::string canvas_name("README"), tp_name("README");

	//check and see if there is a README canvas
	TPaveText *tp = 0;
	TCanvas *c1 = 0;
/*	TCanvas *c1 = dynamic_cast<TCanvas*> (f.Get(canvas_name.c_str()));
	if (c1 != NULL)
	{
		c1->ls();
		tp = dynamic_cast<TPaveText*> (c1->GetPrimitive(tp_name.c_str()));
		if (tp != NULL)
		{
			tp->Print();
		}
	}
*/

	if (c1 == 0)
	{
		c1 = new TCanvas();
		c1->SetTitle(canvas_name.c_str());
		c1->SetName(canvas_name.c_str());
	}
	if (tp == 0)
	{
		tp = new TPaveText(0.1,0.1,0.9,0.9,"NDC TR");
		tp->SetLabel(tp_name.c_str());
		tp->SetName(tp_name.c_str());
	}	

	for (unsigned i = 0; i < lines.size(); ++i)
	{
		tp->AddText(lines.at(i).c_str());
	}

	c1->Draw();
	tp->Draw();

	c1->Write();
	f.ls();
	f.Close();

}
