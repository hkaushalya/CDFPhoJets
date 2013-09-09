#include "TText.h"
#include <iostream>
#include "TCanvas.h"
#include <string>
#include "TPaveText.h"

void TText()
{
	std::string str("msg text fdsf safh akjdfas laf hads flahsd fjhslf as flafh daslj hsdf jkahfsd ");

	TText *t = new TText(0.2,0.2,str.c_str());

	t->Draw();


}
