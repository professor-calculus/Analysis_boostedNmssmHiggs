#include <math.h>

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

double invMass_v1(double E1, double E2, double px1, double px2, double py1, double py2, double pz1, double pz2)
{
	double invmass = sqrt( E1*E2 - (px1*px2 + py1*py2 + pz1*pz2) );
	return invmass;
}

// adds two pt vectors together in the (2d) transverse plane
// output: first element is magnitude, second element is phi
std::vector<double> addTwoPtVectors(double pt1, double phi1, double pt2, double phi2)
{
	std::vector<double> output;
	double xpt = pt1*cos(phi1) + pt2*cos(phi2);
	double ypt = pt1*sin(phi1) + pt2*sin(phi2);
	output.push_back(sqrt(pow(xpt,2)+pow(ypt,2)));
	double phiOut = atan(ypt/xpt);
    if ( xpt < 0 && ypt > 0 ) phiOut = phiOut + M_PI;
    if ( xpt < 0 && ypt < 0 ) phiOut = phiOut - M_PI;
	output.push_back(phiOut);
	return output;
}