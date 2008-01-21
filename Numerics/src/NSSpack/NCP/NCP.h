/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

#ifndef NCP_H
#define NCP_H

/*!\file NCP.h
  Functions related to NCP formulation and solvers.
  *
  *
  * It provides routines to compute Fischer-Burmeister function and its jacobian,
  * written as a NCP-function:
  *
  *  \f[
  *    \phi(z,F(z)) = \sqrt( z^2 + F(z)^2) - z - F(z)
  *  \f]
  *
  * For details see the paper of Kanzow and Kleinmichel, "A New Class of Semismooth Newton-type Methods for Nonlinear
  * Complementarity Problems", Computational Optimization and Applications 11, 227-251 (1998).
  *
  * The notations below are more or less those of this paper.
  *
  * Functions:
  *
  * phi_FB(int size, double* z, double* F, double* phiVector)
  *
  * jacobianPhi_FB(int size, double* z, double* F, double* jacobianF, double* phiVector)
  *
  * \author Houari Khenous, Franck Perignon last modification (13/12/2007)
  *
  */

#include "SparseBlockMatrix.h"
#ifdef __cplusplus
extern "C" {
#endif

  /** NCP Fischer Burmeister function, \f$ \phi(z,F(z)) \f$
      \param size of vector z
      \param vector z
      \param vector F(z)
      \param vector \f$ \phi(z,F(z)) \f$, in-out arg.
  */
  void phi_FB(int, double*, double*, double*);

  /** Jacobian of NCP Fischer Burmeister function, \f$ \nabla_z \phi(z,F(z)) \f$
      \param size of vector z
      \param vector z
      \param vector F(z)
      \param \f$ \nabla_z F(z) \f$
      \param \f$ \nabla_z \phi(z,F(z)) \f$, in-out arg.
  */
  void jacobianPhi_FB(int, double*, double*, double*, double*);

  /**
   * This function checks the validity of the vector z as a solution \n
   * of the LCP : \n
   * \f$
   *    0 \le z \perp Mz + q \ge 0
   * \f$
   * \author Houari Khenous
   \warning temporary function - To be reviewed
  */
  void NCP_compute_error(int n, double *vec , double *q , double *z , int verbose, double *w, double *err);

  /** This function adapts the NCP_compute_error routine for M saved as a SparseBlockStructuredMatrix.
   * TEMPORARY FUNCTION, used to compare pfc3D and FrictionContact3D functions.
   * \author Franck Perignon
   */
  void NCP_block_compute_error(int n, SparseBlockStructuredMatrix *M , double *q , double *z , int verbose, double *w, double *err);

#ifdef __cplusplus
}
#endif

#endif
