
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "evolution.hpp"
#include "qcd.hpp"
#include "transport.hpp"




//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                         EVOLUTION EQUATIONS                      ::
//                                                                  ::
//     Derivatives of hydrodynamic quantities with respect to       ::
//	   longitudinal proper time for bjorken flow. 		            ::
//                                                                  ::
//        dTtt_dtau	   dTtx_dtau	dTty_dtau    dTtn_dtau          ::
//																	::
//		  dpi_dtau                                                  ::
//																	::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


double dTtt_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pi, double Pi, double ut, double ux, double uy, double un, double e, double p, double tau)
{

	//double T = effectiveTemperature(e);
	//double piNS = 4.0 * (e+p) / (3.0*T*tau) * shearViscosityToEntropyDensity(T);
	//double PiNS = - (e+p) / (tau*T) * bulkViscosityToEntropyDensity(T);
	//double Tnn = (e+p+PiNS)*un*un + (p+PiNS-piNS)/(tau*tau);


	double Tnn = (e+p+Pi)*un*un + (p+Pi-pi)/(tau*tau);

	double Tttdot = - (Ttt + tau*tau*Tnn) / tau;

	return Tttdot;
}


double dTtx_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pi, double Pi, double ut, double ux, double uy, double un, double e, double p, double tau)
{
	return - Ttx / tau;
}


double dTty_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pi, double Pi, double ut, double ux, double uy, double un, double e, double p, double tau)
{
	return - Tty / tau;
}


double dTtn_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pi, double Pi, double ut, double ux, double uy, double un, double e, double p, double tau)
{
	return - Ttn / tau;
}


// time derivative of shear component pi = -tau^2 * pinn
double dpi_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pi, double Pi, double ut, double ux, double uy, double un, double e, double p, double tau)
{
	// lattice qcd
	double T = effectiveTemperature(e); 		 // temperature
	double s = (e+p)/T; 						 // entropy density

	double etas = shearViscosityToEntropyDensity(T);  // specific shear viscosity
	double eta = s * etas;                        	  // shear viscosity

	// Model 1: fixed mass and m/T << 1
	double taupiInv = 0.2 * T / etas;                // shear relaxation rate
	double taupipi = 10.0 / 7.0;
	double deltapipi = 4.0 / 3.0;
	double lambdapiPi = 1.2;


	// Model 2: quasiparticle model
	// double taupiInv = beta_shear(T)/eta;
	// double taupipi = tauSS(T);
	// double deltapipi = deltaSS(T);
	// double lambdapiPi = lambdaSB(T);


	// shear stress relaxation equation
	double pidot = taupiInv * (-pi + 4.0*eta/3.0/tau) - (taupipi/3.0 + deltapipi)*pi/tau + 2.0*lambdapiPi*Pi/3.0/tau;

	return pidot;
}


// time derivative of bulk pressure
double dPi_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pi, double Pi, double ut, double ux, double uy, double un, double e, double p, double tau)
{
	// lattice qcd
	double T = effectiveTemperature(e); 		 // temperature
	double s = (e+p)/T; 						 // entropy density
	double cs2 = speedOfSoundSquared(e); 	     // speed of sound squared
	double b2 = (1.0/3.0 - cs2);

	if(b2 == 0)
	{
		return 0.0; // conformal eos
	}

	double etas = shearViscosityToEntropyDensity(T);  // specific shear viscosity
	double eta = s * etas;

	double zetas = bulkViscosityToEntropyDensity(T);      // specific bulk viscosity
	double zeta = s * zetas;						      // bulk viscosity

	// Model 1: fixed mass and m/T << 1
	double tauPiInv = 15.0 * b2 * b2 * T / zetas;      // bulk relaxation rate
	double deltaPiPi = 2.0 / 3.0;
	double lambdaPipi = 1.6 * b2;

	//double taupiInv = beta_shear(T)/eta;

	// Model 2: quasiparticle model
	// double tauPiInv = beta_bulk(T)/zeta;
	// double deltaPiPi = deltaBB(T);
	// double lambdaPipi = lambdaBS(T);

	// bulk relaxation equation
	double Pidot = -tauPiInv * (Pi + zeta/tau) + (-deltaPiPi*Pi + lambdaPipi*pi) / tau;

	//double Pidot = -taupiInv * Pi - beta_bulk(T)/tau + (-deltaPiPi*Pi + lambdaPipi*pi) / tau;

	return Pidot;
}




















