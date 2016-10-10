double delPhi(double phi1, double phi2)
{
	double dPhi = phi1 - phi2;
	if (dPhi>M_PI) dPhi = dPhi - 2*M_PI;
	if (dPhi<-M_PI) dPhi = dPhi + 2*M_PI;
	return dPhi;
}

double delEta(double eta1, double eta2)
{
	double dEta = eta1 - eta2;
	return dEta;
}
    
double delR(double dPhi, double dEta)
{
  double dR = sqrt(dPhi*dPhi + dEta*dEta);
  return dR;
}