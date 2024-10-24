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
/*! \file EqualityConditionNSL.hpp

*/
#ifndef EQUALITYCONDITIONNSLAW_H
#define EQUALITYCONDITIONNSLAW_H

#include "NonSmoothLaw.hpp"


/** Equality NonSmoothLaw
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 27, 2004
 *
 *
 **/
class EqualityConditionNSL : public NonSmoothLaw
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(EqualityConditionNSL);

  /** default constructor
   */
  EqualityConditionNSL() {};

public:
  /** basic constructor
  *  \param size of the non smooth law
  */
  EqualityConditionNSL(unsigned int size);

  /** Destructor */
  ~EqualityConditionNSL();


  /** print the data to the screen
  */
  inline void display()const {};

  /** Visitors hook
   */
  ACCEPT_STD_VISITORS();
};

TYPEDEF_SPTR(EqualityConditionNSL)

#endif // EQUALITYCONDITIONNSLAW_H
