
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

#if (KINETIC == 1)
	double taupiInv = beta_shear(T)/eta;
	double taupipi = tauSS(T);
	double deltapipi = deltaSS(T);
	double lambdapiPi = lambdaSB(T);
#else
	double taupiInv = 0.2 * T / etas;
	double taupipi = 10.0 / 7.0;
	double deltapipi = 4.0 / 3.0;
	double lambdapiPi = 1.2;
#endif

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

#if (KINETIC == 1)
	double tauPiInv = beta_bulk(T)/zeta;
	double deltaPiPi = deltaBB(T);
	double lambdaPipi = lambdaBS(T);
#else
	double tauPiInv = 15.0 * b2 * b2 * T / zetas;
	double deltaPiPi = 2.0 / 3.0;
	double lambdaPipi = 1.6 * b2;
#endif

	double Pidot = -tauPiInv * (Pi + zeta/tau) + (-deltaPiPi*Pi + lambdaPipi*pi) / tau;

	return Pidot;
}




















