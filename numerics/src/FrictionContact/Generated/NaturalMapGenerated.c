#include "NaturalMapGenerated.h"
#include "assert.h"
#include "op3x3.h"

void fc3d_NaturalMapFunctionGenerated(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *f,
  double *A,
  double *B)
{
  double result[3];
  double resultd[3*NBDirsMax];

  assert(reaction);
  assert(velocity);
  assert(rho);

  SET3(reaction);
  SET3(velocity);
  SET3(rho);

  double und[6] =  { 1., 0., 0., 0., 0., 0. };
  double ut1d[6] = { 0., 1., 0., 0., 0., 0. };
  double ut2d[6] = { 0., 0., 1., 0., 0., 0. };

  double rnd[6] =  { 0., 0., 0., 1., 0., 0. };
  double rt1d[6] = { 0., 0., 0., 0., 1., 0. };
  double rt2d[6] = { 0., 0., 0., 0., 0., 1. };

  Fnat_dv(mu, *velocity0, und, *velocity1, ut1d, *velocity2, ut2d,
          *reaction0, rnd, *reaction1, rt1d, *reaction2, rt2d,
          result, (double(*)[6]) resultd, 6);

  if (f && A && B)
  {
      cpy3(result,     f);

      for(unsigned i=0; i<3; ++i)
      {
        for(unsigned j=0; j<3; ++j)
        {
          A[i + 3*j] = resultd[6*i + j];
          B[i + 3*j] = resultd[6*i + j + 3];
        }
      }
  }
  else
  {
    if (f)
    {
      cpy3(result, f);
    }

    if (A && B)
    {
      for(unsigned i=0; i<3; ++i)
      {
        for(unsigned j=0; j<3; ++j)
        {
          A[i + 3*j] = resultd[6*i + j];
          B[i + 3*j] = resultd[6*i + j + 3];
        }
      }
    }
  }
}
