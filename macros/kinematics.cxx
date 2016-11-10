#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>

void kinematics(){

	/*
	We calculate kintematic properties of the following interaction:
	pp -> squark + squark

	squark -> quark + NLSP;
	NLSP -> LSP + higgs;
	higgs -> bb;

	interesting properties:

	quark: energy/momentum, transverse momentum
	NLSP: energy, momentum, transverse momentum

	LSP:
	in rest frame of NLSP - energy, momentum
	in lab frame - energy, momentum, transverse momentum

	Higgs:
	in rest frame of NLSP - energy, momentum
	in lab frame - energy, momentum, transverse momentum

	Splitting angle of the bbar pair
	(if LSP pt small, it will carry the NLSP momentum)
	*/


	// --- USER INPUT --- //
	// set the four mass parameters
	double mass_squark = 1000.0; // in GeV
	double mass_higgs = 83.0;
	double mass_delta = 1.0; // = mass_NLSP - mass_higgs - mass_LSP
	double mass_ratio = 0.93; // = mass_higgs / mass_NLSP (must be less than one) 
	// ------------------ //
	std::cout << std::endl;
	std::cout << "USER INPUTS" << std::endl;
	std::cout << "Mass Squark: " << mass_squark << std::endl;
	std::cout << "Mass Higgs: " << mass_higgs << std::endl;
	std::cout << "Mass Splitting: " << mass_delta << std::endl;
	std::cout << "Higgs Mass / NLSP mass: " << mass_ratio << std::endl;
	std::cout << std::endl;


	// calculate the LSP and NLSP masses from input
	double mass_NLSP = mass_higgs / mass_ratio;
	double mass_LSP = mass_NLSP - mass_higgs - mass_delta;
	std::cout << "CALCULATIONS" << std::endl;
	std::cout << "Mass NLSP: " << mass_NLSP << std::endl;
	std::cout << "Mass LSP: " << mass_LSP << std::endl;
	std::cout << std::endl;


	// squark -> quark + NLSP
	double energy_quark = (mass_squark * mass_squark - mass_NLSP * mass_NLSP) / (2 * mass_squark);
	double momentum_quark = energy_quark;
	double energy_NLSP = mass_squark - energy_quark;
	double momentum_NLSP = momentum_quark;
	std::cout << "Energy Quark: " << energy_quark << std::endl;
	std::cout << "Momentum Quark: " << momentum_quark << std::endl;
	std::cout << "Energy NLSP: " << energy_NLSP << std::endl;
	std::cout << "Momentum NLSP: " << momentum_NLSP << std::endl;
	std::cout << std::endl;

	// TODO TRANSVERSE MOMENTUM THINGZZ





	// NLSP -> LSP + higgs
	double boost_gamma = energy_NLSP / mass_NLSP;
	double boost_betaGamma = momentum_NLSP / mass_NLSP;

	// TODO: process this angle in a smarter way...eg tfunction with a loop or something...
	double boost_angleLSP = 0; // angle between boost axis and the momentum of the LSP
	double boost_angleHiggs = M_PI - boost_angleLSP;

	// LSP: rest frame
	double energy_LSP_restFrame = (mass_NLSP * mass_NLSP + mass_LSP * mass_LSP - mass_higgs * mass_higgs) / (2 * mass_NLSP);
	double momentum_LSP_restFrame = sqrt(energy_LSP_restFrame * energy_LSP_restFrame - mass_LSP * mass_LSP);

	// higgs: rest frame
	double energy_higgs_restFrame = mass_NLSP - energy_LSP_restFrame;
	double momentum_higgs_restFrame = momentum_LSP_restFrame;

	std::cout << "rest frame" << std::endl;
	std::cout << "Energy LSP: " << energy_LSP_restFrame << std::endl;
	std::cout << "Momentum LSP: " << momentum_LSP_restFrame << std::endl;
	std::cout << "Energy Higgs: " << energy_higgs_restFrame << std::endl;
	std::cout << "Momentum Higgs: " << momentum_higgs_restFrame << std::endl;
	std::cout << std::endl;
	std::cout << "boost to lab frame params" << std::endl;
	std::cout << "Gamma: " << boost_gamma << std::endl;
	std::cout << "betaGamma: " << boost_betaGamma << std::endl;
	std::cout << std::endl;

	// LSP: lab frame
	double energy_LSP = boost_gamma * energy_LSP_restFrame + boost_betaGamma * momentum_LSP_restFrame * cos(boost_angleLSP);
	double momentum_LSP = sqrt(energy_LSP * energy_LSP - mass_LSP * mass_LSP);

	// higgs: lab frame
	double energy_higgs = boost_gamma * energy_higgs_restFrame + boost_betaGamma * momentum_higgs_restFrame * cos(boost_angleHiggs);
	double momentum_higgs = sqrt(energy_higgs * energy_higgs - mass_higgs * mass_higgs);

	std::cout << "lab frame" << std::endl;
	std::cout << "Energy LSP: " << energy_LSP << std::endl;
	std::cout << "Momentum LSP: " << momentum_LSP << std::endl;
	std::cout << "Energy higgs: " << energy_higgs << std::endl;
	std::cout << "Momentum higgs: " << momentum_higgs << std::endl;
	std::cout << std::endl;
	// TODOtrnasverse momentum thingzzz...random dist




	// higgs -> bbar
	double seperation_bb = 2 * mass_higgs / momentum_higgs; //TODO transverse momentum
	std::cout << "Average dR(h->bb): " << seperation_bb << std::endl;


































}