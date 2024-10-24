/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include "numerics_verbose.h"
#include "op3x3.h"
#include "SparseBlockMatrix.h"
#include "fc3d_Solvers.h"
#include "FrictionContactProblem.h"
#include "fc3d_compute_error.h"
#include "FischerBurmeisterGenerated.h"
#include "fc3d_nonsmooth_Newton_solvers.h"
#include "fc3d_nonsmooth_Newton_FischerBurmeister.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "Friction_cst.h"
#include "SiconosLapack.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

void fc3d_FischerBurmeisterFunction(
  unsigned int problemSize,
  FischerBurmeisterFun3x3Ptr computeACFun3x3,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *result,
  double *A,
  double *B)
{
  assert(reaction);
  assert(velocity);
  assert(rho);
  assert(mu);

  assert(problemSize / 3 > 0);
  assert(problemSize % 3 == 0);

  unsigned int i;
  for (i = 0; i < problemSize; i += 3)
  {

    computeACFun3x3(reaction, velocity, *mu, rho, result, A, B);

    reaction += 3;
    velocity += 3;
    mu++;
    rho += 3;

    if (result)
      result += 3;

    if (A)
      A += 9;

    if (B)
      B += 9;

  }

}


int fc3d_nonsmooth_Newton_FischerBurmeister_compute_error(
    FrictionContactProblem* problem,
    double *z , double *w, double tolerance,
    SolverOptions * options, double * error)
{

  double *A = NULL;
  double *B = NULL;

  unsigned int problemSize = 3 * problem->numberOfContacts;

  double *rho = (double*) malloc(problemSize*sizeof(double));
  double *F = (double *) malloc(problemSize*sizeof(double));

  FischerBurmeisterFun3x3Ptr computeACFun3x3;

  switch (options->iparam[10])
  {
  case 0:
  {

    computeACFun3x3 = &fc3d_FischerBurmeisterFunctionGenerated;
    break;
  }
  }

  fc3d_FischerBurmeisterFunction(
    problemSize,
    computeACFun3x3,
    z, w,
    problem->mu, rho,
    F, A, B);

  *error=0.;
  for(unsigned int i=0; i<problemSize;
      i+=3)
  {
    *error += sqrt(F[i]*F[i] + F[i+1]*F[i+1] + F[i+2]*F[i+2]);
  }

  *error /= (problem->numberOfContacts + 1);

  free(F);
  free(rho);

  if (*error > tolerance)
  {
    if (verbose > 1)
      printf(" Numerics - fc3d_compute_error: error = %g > tolerance = %g.\n",
             *error, tolerance);
    return 1;
  }
  else
  {
    return 0;
  }
}

int fc3d_nonsmooth_Newton_FischerBurmeister_setDefaultSolverOptions(
  SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the default solver options for the NSN_FB Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_NSN_FB;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 14;
  options->dSize = 14;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  options->iparam[0] = 200;
  options->iparam[1] = 1;
  options->iparam[3] = 100000; /* nzmax*/
  options->iparam[5] = 1;
  options->iparam[7] = 1;      /* erritermax */
  options->dparam[0] = 1e-3;
  options->dparam[3] = 1;      /* default rho */

  options->iparam[8] = -1;     /* mpi com fortran */
  options->iparam[10] = 0;
  options->iparam[11] = 0;     /* 0 GoldsteinPrice line search, 1 FBLSA */
  options->iparam[12] = 100;   /* max iter line search */

#ifdef WITH_MUMPS
  options->iparam[13] = 1;
#else
  options->iparam[13] = 0;     /* Linear solver used at each Newton iteration. 0: cs_lusol, 1 mumps */
#endif

  options->internalSolvers = NULL;

#ifdef HAVE_MPI
  options->solverData = MPI_COMM_NULL;
#endif

  return 0;
}


typedef struct
{
  FischerBurmeisterFun3x3Ptr computeACFun3x3;
} FischerBurmeisterParams;

void nonsmoothEqnFischerBurmeisterFun(void* arg,
                                      unsigned int problemSize,
                                      double* reaction,
                                      double* velocity,
                                      double* mu,
                                      double* rho,
                                      double* result,
                                      double* A,
                                      double* B);
void nonsmoothEqnFischerBurmeisterFun(void* arg,
                                      unsigned int problemSize,
                                      double* reaction,
                                      double* velocity,
                                      double* mu,
                                      double* rho,
                                      double* result,
                                      double* A,
                                      double* B)
{
  FischerBurmeisterParams* acparams_p = (FischerBurmeisterParams *) arg;

  fc3d_FischerBurmeisterFunction(problemSize,
                                         acparams_p->computeACFun3x3,
                                         reaction,
                                         velocity,
                                         mu,
                                         rho,
                                         result,
                                         A,
                                         B);
}




void fc3d_nonsmooth_Newton_FischerBurmeister(
  FrictionContactProblem* problem,
  double *reaction,
  double *velocity,
  int *info,
  SolverOptions *options)
{
  assert(problem);
  assert(reaction);
  assert(velocity);
  assert(info);
  assert(options);

  assert(problem->dimension == 3);

  assert(options->iparam);
  assert(options->dparam);

  assert(problem->q);
  assert(problem->mu);
  assert(problem->M);

  assert(!options->iparam[4]); // only host

  FischerBurmeisterParams acparams;

  switch (options->iparam[10])
  {
  case 0:
  {
    acparams.computeACFun3x3 = &fc3d_FischerBurmeisterFunctionGenerated;
    break;
  }
  }

  fc3d_nonsmooth_Newton_solvers equation;

  equation.problem = problem;
  equation.data = (void *) &acparams;
  equation.function = &nonsmoothEqnFischerBurmeisterFun;

  fc3d_nonsmooth_Newton_solvers_solve(&equation, reaction, velocity, info,
                                   options);

}
