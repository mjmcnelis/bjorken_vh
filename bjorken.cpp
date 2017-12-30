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
#include "transport.hpp"
//#include <gsl/gsl_sf.h>

#define GEV_TO_INVERSE_FM 5.067731


int main()
{
	// Bjorken flow


	// input parameters
	const double T0 = 0.6 * GEV_TO_INVERSE_FM;  // initial temperature in fm^-1
	const double tau0 = 0.25;					// initial time in fm
	const double tauf = 50.0;					// final time in fm


	// initialize temperature
	double T = T0;


	// initial flow profile
	double ut = 1.0;
	double ux = 0.0;
	double uy = 0.0;
	double un = 0.0;


	// initialize energy density and pressure  (units = [fm^-4])
	const double e0 = equilibriumEnergyDensity(T0);
	double e = e0;
	const double p0 = equilibriumPressure(e0);
	double p = p0;


	double cs2 = speedOfSoundSquared(e0);


	double s0 = (e0+p0)/T0; 						    // initial entropy density
	double zetas0 = bulkViscosityToEntropyDensity(T0);  // initial specific bulk viscosity
	double etas0 = shearViscosityToEntropyDensity(T0);  // initial specific shear viscosity


	double taupi = (s0*etas0) / beta_shear(T0);             // quasiparticle model
	double taubulk = (s0*zetas0) / beta_bulk(T0);

	cout << sqrt(1.5)*taupi/tau0 << endl;
	cout << taubulk/tau0 << endl;

	//double taupi = (s0*etas0) / (0.2*(e0+p0));            // m/T << 1 fixed mass model
	//double taubulk = (s0*zetas0) / (15.0*pow(1.0/3.0-cs2,2)*(e0+p0));


	// initial Tt\mu components (units = [fm^-4])
	double Ttt = e0;
	double Ttx = 0.0;
	double Tty = 0.0;
	double Ttn = 0.0;

	// initialize shear stress: pi = - tau^2 * pinn (units = [fm^-4])
	double pi0 = 4.0 * s0 / (3.0 * tau0) * etas0; // (Navier Stokes)
	double pi = 0.0*pi0;

	// initialize bulk pressure (units = [fm^-4])
	double Pi0 = - s0 * zetas0 / tau0;  // (Navier Stokes)
	double Pi = 0.0*Pi0; // set to zero for now


// Glasma initial conditions:

	double PL = 0.00149925 * p;
	double PT = 1.49925 * p;
	pi = 2.0*(PT-PL)/3.0;
	Pi = (2.0*PT/3.0 + PL/3.0 - p);


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
	ofstream RpiInvplot, RbulkInvplot;
	ofstream taupiplot, taubulkplot;

	eplot.open("eplot_vh.dat", ios::out);
	piplot.open("piplot_vh.dat", ios::out);
	bulkplot.open("bulkplot_vh.dat", ios::out);
	plptplot.open("plptplot_vh.dat", ios::out);

	RpiInvplot.open("RpiInvplot_vh.dat", ios::out);
	RbulkInvplot.open("RbulkInvplot_vh.dat", ios::out);

	taupiplot.open("taupiplot_vh.dat", ios::out);
	taubulkplot.open("taubulkplot_vh.dat", ios::out);


	eplot << "tau [fm]" << "\t\t" << "e/e0" << endl << setprecision(6) << tau << "\t\t" << e/e0 << endl;
	piplot << "tau [fm]" << "\t\t" << "pi [fm^-4]" << endl << setprecision(6) << tau << "\t\t" << pi << endl;
	bulkplot << "tau [fm]" << "\t\t" << "Pi [fm^-4]" << endl << setprecision(6) << tau << "\t\t" << Pi << endl;
	plptplot << "tau [fm]" << "\t\t" << "PL/PT" << endl << setprecision(6) << tau << "\t\t" << (p + Pi - pi) / (p + Pi + 0.5*pi) << endl;

	RpiInvplot << "tau [fm]" << "\t\t" << "R_pi^-1" << endl << setprecision(6) << tau << "\t\t" << sqrt(1.5) * pi / p << endl;
	RbulkInvplot << "tau [fm]" << "\t\t" << "R_Pi^-1" << endl << setprecision(6) << tau << "\t\t" << Pi / p << endl;

	taupiplot << "tau [fm]" << "\t\t" << "tau_pi" << endl << setprecision(6) << tau << "\t\t" << taupi << endl;
	taubulkplot << "tau [fm]" << "\t\t" << "tau_Pi" << endl << setprecision(6) << tau << "\t\t" << taubulk << endl;



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

		T = effectiveTemperature(e);
		cs2 = speedOfSoundSquared(e);

		//piNS = 4.0 * (e+p) / (3.0*T*tau) * shearViscosityToEntropyDensity(T);
		//PiNS = - (e+p) / (tau*T) * bulkViscosityToEntropyDensity(T);


		// quasiparticle model
		taupi = (e+p) * shearViscosityToEntropyDensity(T) / (T*beta_shear(T));
		taubulk = (e+p) * bulkViscosityToEntropyDensity(T) / (T*beta_bulk(T));


		// m/T << 1 model
		//taupi = 5.0 * shearViscosityToEntropyDensity(T) / T;
		//taubulk = bulkViscosityToEntropyDensity(T) / (15.0*T*pow(1.0/3.0-cs2,2));


		// write updated energy density to file
		if((i+1)%timesteps_per_write == 0)
		{
			eplot << setprecision(6) << tau << "\t\t" << e/e0 << "\t\t" << endl;
			piplot << setprecision(6) << tau << "\t\t" << pi << "\t\t" << endl;
			bulkplot << setprecision(6) << tau << "\t\t" << Pi << "\t\t" << endl;
			plptplot << setprecision(6) << tau << "\t\t" << (p + Pi - pi) / (p + Pi + 0.5*pi) << "\t\t" << endl;

			RpiInvplot << setprecision(6) << tau << "\t\t" << sqrt(1.5) * pi / p << "\t\t" << endl;
			RbulkInvplot << setprecision(6) << tau << "\t\t" << Pi / p << "\t\t" << endl;

			taupiplot << setprecision(6) << tau << "\t\t" << taupi << "\t\t" << endl;
			taubulkplot << setprecision(6) << tau << "\t\t" << taubulk  << "\t\t" << endl;

		}
	}



	// close plot data files

	eplot.close();
	piplot.close();
	bulkplot.close();
	plptplot.close();

	RpiInvplot.close();
	RbulkInvplot.close();

	taupiplot.close();
	taubulkplot.close();




	// free memory
	//printf("Freeing memory...");
	printf("Done\n\n");

	return 0;
}






