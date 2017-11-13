#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string>
#include <string.h>
#include <iostream>
#include <iomanip>
using namespace std;
#include <sstream>
#include <fstream>
#include "qcd.hpp"
#include "evolution.hpp"
#include "inferredvariables.hpp"
//#include <gsl/gsl_sf.h>

#define GEV_TO_INVERSE_FM 5.067731


int main()
{
	// Bjorken flow

	// input parameters
	const double T0 = 0.6 * GEV_TO_INVERSE_FM;  // initial temperature in fm^-1
	const double tau0 = 0.25;					// initial time in fm
	const double tauf = 20.0;					// final time in fm



	// initial flow profile
	double ut = 1.0;
	double ux = 0.0;
	double uy = 0.0;
	double un = 0.0;


	// initial energy density and pressure  (units = [fm^-4])
	const double e0 = equilibriumEnergyDensity(T0);
	double e = e0;
	const double p0 = equilibriumPressure(e0);
	double p = p0;

	double s0 = (e0+p0)/T0; 						    // initial entropy density
	double zetas0 = bulkViscosityToEntropyDensity(T0);  // initial specific bulk viscosity


	// initial Tt\mu components (units = [fm^-4])
	double Ttt = e0;
	double Ttx = 0.0;
	double Tty = 0.0;
	double Ttn = 0.0;


	// initial shear stress: pi = - tau^2 * pinn (units = [fm^-4])
	double pi0 = 4.0 * s0 / (3.0 * tau0) * ETAS; // (Navier Stokes)
	double pi = 0.0*pi0;


	// initial bulk pressure (units = [fm^-4])
	double Pi0 = - s0 * zetas0 / tau0;  // (Navier Stokes)
	double Pi = 0.0*Pi0; // set to zero for now


	// intermediate and end values for heun's rule
	double Ttt_mid, Ttt_end;
	double Ttx_mid, Ttx_end;
	double Tty_mid, Tty_end;
	double Ttn_mid, Ttn_end;
	double pi_mid, pi_end;
	double Pi_mid, Pi_end;




	// initial time variable set to tau0, number of steps, stepsize
	double tau = tau0;
	const double dtau = 0.001;
	const int n = floor((tauf - tau0) / dtau);
	const int timesteps_per_write = 10;

	// Data files for plots
	ofstream eplot, piplot, bulkplot, plptplot;

	eplot.open("eplot.dat", ios::out);
	piplot.open("piplot.dat", ios::out);
	bulkplot.open("bulkplot.dat", ios::out);
	plptplot.open("plptplot.dat", ios::out);

	eplot << "tau [fm]" << "\t\t" << "e/e0" << endl << setprecision(6) << tau << "\t\t" << 1.0 << endl;
	piplot << "tau [fm]" << "\t\t" << "pi/p" << endl << setprecision(6) << tau << "\t\t" << 0.0 << endl;
	bulkplot << "tau [fm]" << "\t\t" << "Pi [Gev/fm^3]" << endl << setprecision(6) << tau << "\t\t" << Pi / GEV_TO_INVERSE_FM << endl;
	plptplot << "tau [fm]" << "\t\t" << "PL/PT" << endl << setprecision(6) << tau << "\t\t" << (p + Pi - pi) / (p + Pi + 0.5*pi) << endl;



	// start evolution
	for(int i = 0; i < n; i++)
	{
		// compute intermediate values with Euler step
		Ttt_mid = Ttt + dtau * dTtt_dtau(Ttt, Ttx, Tty, Ttn, pi, Pi, ut, ux, uy, un, e, p, tau);
		Ttx_mid = Ttx + dtau * dTtx_dtau(Ttt, Ttx, Tty, Ttn, pi, Pi, ut, ux, uy, un, e, p, tau);
		Tty_mid = Tty + dtau * dTty_dtau(Ttt, Ttx, Tty, Ttn, pi, Pi, ut, ux, uy, un, e, p, tau);
		Ttn_mid = Ttn + dtau * dTtn_dtau(Ttt, Ttx, Tty, Ttn, pi, Pi, ut, ux, uy, un, e, p, tau);
		pi_mid = pi + dtau * dpi_dtau(Ttt, Ttx, Tty, Ttn, pi, Pi, ut, ux, uy, un, e, p, tau);
		Pi_mid = Pi + dtau * dPi_dtau(Ttt, Ttx, Tty, Ttn, pi, Pi, ut, ux, uy, un, e, p, tau);


		// find intermediate inferred variables
		get_inferred_variables(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pi_mid, Pi_mid, &ut, &ux, &uy, &un, &e, &p, tau + dtau);


		// add Euler step with respect to the intermediate value
		Ttt_end = Ttt_mid + dtau * dTtt_dtau(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pi_mid, Pi_mid, ut, ux, uy, un, e, p, tau + dtau);
		Ttx_end = Ttx_mid + dtau * dTtx_dtau(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pi_mid, Pi_mid, ut, ux, uy, un, e, p, tau + dtau);
		Tty_end = Tty_mid + dtau * dTty_dtau(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pi_mid, Pi_mid, ut, ux, uy, un, e, p, tau + dtau);
		Ttn_end = Ttn_mid + dtau * dTty_dtau(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pi_mid, Pi_mid, ut, ux, uy, un, e, p, tau + dtau);
		pi_end = pi_mid + dtau * dpi_dtau(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pi_mid, Pi_mid, ut, ux, uy, un, e, p, tau + dtau);
		Pi_end = Pi_mid + dtau * dPi_dtau(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pi_mid, Pi_mid, ut, ux, uy, un, e, p, tau + dtau);

		// increase time step
		tau += dtau;


		// Heun's Rule (average initial and end)
		// updated variables at new time step
		Ttt = 0.5 * (Ttt + Ttt_end);
		Ttx = 0.5 * (Ttx + Ttx_end);
		Tty = 0.5 * (Tty + Tty_end);
		Ttn = 0.5 * (Ttn + Ttn_end);
		pi = 0.5 * (pi + pi_end);
		Pi = 0.5 * (Pi + Pi_end);


		// find inferred variables at new time step
		get_inferred_variables(Ttt, Ttx, Tty, Ttn, pi, Pi, &ut, &ux, &uy, &un, &e, &p, tau);


		// write updated energy density to file
		if((i+1)%timesteps_per_write == 0)
		{
			eplot << setprecision(6) << tau << "\t\t" << e / e0 << "\t\t" << endl;
			piplot << setprecision(6) << tau << "\t\t" << pi / p0 << "\t\t" << endl;
			bulkplot << setprecision(6) << tau << "\t\t" << Pi / GEV_TO_INVERSE_FM << "\t\t" << endl;
			plptplot << setprecision(6) << tau << "\t\t" << (p + Pi - pi) / (p + Pi + 0.5*pi) << "\t\t" << endl;
		}
	}



	// close plot data files
	eplot.close();
	piplot.close();
	bulkplot.close();
	plptplot.close();



	// free memory
	printf("Freeing memory...");
	printf("done\n\n");

	return 0;
}






