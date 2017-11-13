
#include <stdlib.h>

#ifndef QCD_H

#define QCD_H

//#define CONFORMAL_EOS // if using P_eq = e_eq / 3 or m = 0

// ideal gas of massless quarks and gluons
//#define EOS_FACTOR 15.6269   // Nc=3, Nf=3
#define EOS_FACTOR 13.8997   // Nc=3, Nf=2.5


double equilibriumPressure(double e);

double speedOfSoundSquared(double e);

double effectiveTemperature(double e);

double equilibriumEnergyDensity(double T);

double derivativeEnergyDensityWithRespectToTemperature(double T);

double bulkViscosityToEntropyDensity(double T);

double z_qcd(double T);

double zdmdT(double T);

//double zdmdT(double T);

double z_Quasiparticle(double T);
double mdmdT_Quasiparticle(double T);


#endif