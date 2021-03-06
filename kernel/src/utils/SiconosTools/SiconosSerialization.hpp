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

/*! \file SiconosSerialization.hpp
  serialization for Siconos
*/

#ifndef SiconosSerialization_hpp
#define SiconosSerialization_hpp

/** install serialization hooks. Must be used inside a protected zone
    of class definition
    \parameter a class name
 */
namespace boost
{
namespace serialization
{
class access;
}
}

#define ACCEPT_SERIALIZATION(CLASS)                             \
  typedef void serializable;                                    \
  template<typename Archive>                                    \
  friend void siconos_io(Archive&, CLASS&, const unsigned int); \
  friend class boost::serialization::access

#endif
