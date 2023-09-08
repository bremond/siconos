#pragma once

#include "siconos/utils/environment.hpp"
#include "siconos/storage/storage.hpp"
#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/algebra/numerics.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/model/nslaws.hpp"
#include "siconos/simul/topology.hpp"
#include "siconos/simul/interaction.hpp"
#include "siconos/simul/one_step_nonsmooth_problem.hpp"
#include "siconos/simul/one_step_integrator.hpp"
#include "siconos/simul/time_discretization.hpp"
#include "siconos/simul/time_stepping.hpp"
#include <charconv>
#include <tuple>
#include <functional>
#include <typeinfo>
#if defined(__clang__)
#include <boost/hana/experimental/type_name.hpp>
#endif
namespace siconos
{
}

