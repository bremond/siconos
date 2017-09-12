/* -*- c -*- */
#include <math.h>
#define NBDirsMax 6


void Fnat_dv(double mu, double un, double und[NBDirsMax], double ut1, double 
        ut1d[NBDirsMax], double ut2, double ut2d[NBDirsMax], double rn, double
        rnd[NBDirsMax], double rt1, double rt1d[NBDirsMax], double rt2, double
        rt2d[NBDirsMax], double *result, double (*resultd)[6], int 
             nbdirs);
