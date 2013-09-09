{

	int pmean = 5;
	int x = 3;
	double p = TMath::Poisson(x,pmean);

  //double f = pow(x, mean)*exp(-x)/TMath::Factorial(int(mean+0.001));
  double f = pow(pmean, x)*exp(-pmean)/TMath::Factorial(x);
	cout << "p = " << p << endl;
	cout << "f = " << f << endl;


	double gmean = 4;
	double gsigma = 2;
  //double g = 1/sqrt(2*pi)/par[1]*exp(-pow(x[0]-par[0],2)/2.0/par[1]/par[1]);
  double g = 1/sqrt(2*TMath::Pi())/gsigma*exp(-pow(x[0]-gmean,2)/2.0/gsigma/gsigma);
  double g1 = TMath::Gaus(x, gmean, gsigma,1);
	cout << "g = " << g << endl;
	cout << "g1 = " << g1 << endl;
	
}
