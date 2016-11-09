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
#ifndef HARMONICBC_HPP
#define HARMONICBC_HPP


#include "BoundaryCondition.hpp"

/** \class HarmonicBC
 *  \brief This class models a simple harmonic boundary conditions for
 *   prescribing the velocities in a Dynamical System. A simple
 *   boundary condition is considered to fix a component \f$ j \f$ of
 * the velocity vector, i.e., \f$ v_j(t) = a +  b cos( \omega t+ \phi)\f$ where \f$
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.2.0.
 *  \date (Re-Creation) November 2016
 *
 */
class HarmonicBC : public  BoundaryCondition
{
public:

  /** \fn HarmonicBC(SP::UnsignedIntVector  newVelocityIndices);
   *  \brief Basic constructor
   *  \param newVelocityIndices the indices of the velocity subjected to prescribed velocities
   */

  HarmonicBC(SP::UnsignedIntVector newVelocityIndices,
             double a, double b,
             double omega, double phi) ;


  /** destructor */
  virtual ~HarmonicBC();

  /** default function to compute the precribed velocities
   *  \param  time : the current time
   */
  virtual void computePrescribedVelocity(double time);

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(HarmonicBC);

  /** protected default constructor */
  HarmonicBC(): BoundaryCondition() {};

  /** Constant additive term of the prescribed velocity  */
  double _a;
  /** Constant multiplicative term of the prescribed velocity  */
  double _b;
  /** Constant frequency  */
  double _omega;
  /** Constant phase  */
  double _phi;

  

  //   /*Link to the precribed DynamicalSystem*/
  //   SP::DynamicalSystem _DS;
};

TYPEDEF_SPTR(HarmonicBC)
#endif // HARMONICBC_HPP