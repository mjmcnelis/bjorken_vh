
#include <stdlib.h>

#ifndef EVOLUTION_H

#define EVOLUTION_H

#define ETAS 0.2;  // specific shear viscosity


double dTtt_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pi, double Pi, double ut, double ux, double uy, double un, double e, double p, double tau);

double dTtx_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pi, double Pi, double ut, double ux, double uy, double un, double e, double p, double tau);

double dTty_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pi, double Pi, double ut, double ux, double uy, double un, double e, double p, double tau);

double dTtn_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pi, double Pi, double ut, double ux, double uy, double un, double e, double p, double tau);

double dpi_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pi, double Pi, double ut, double ux, double uy, double un, double e, double p, double tau);

double dPi_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pi, double Pi, double ut, double ux, double uy, double un, double e, double p, double tau);


#endif