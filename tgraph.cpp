#include <benchmark/benchmark.h>

#include <array>
#include <variant>
#include <random>
#include <functional>
#include <algorithm>
#include <iostream>
#include <memory>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>

#include "SiconosKernel.hpp"

static const size_t SIZE = 100000;

std::random_device dev;
std::mt19937 rng(dev());
std::uniform_int_distribution<std::mt19937::result_type> random_type(0,3);
std::uniform_int_distribution<std::mt19937::result_type> random_indice(0,SIZE-1);

template <typename T, size_t N>
std::array<T, N> random_vector_stack(std::function<T(decltype(random_type(rng)))> fun) {
  std::array<T, N> a;
  std::generate(a.begin(), a.end(), [&] {
    return fun(random_type(rng));
  });
  return a;
}

template <typename T, size_t N>
auto random_vector_heap(std::function<T(decltype(random_type(rng)))> fun) {
  std::shared_ptr<std::vector<T> > a(new std::vector<T>());
  a->resize(N);
  std::generate(a->begin(), a->end(), [&] {
      return fun(random_type(rng));
    });
  return a;
}


template<typename T>
auto random_pitem(const T& array)
{
  return array[random_indice(rng)];
}

template<typename T>
auto random_item(const T& array)
{
  return &array[random_indice(rng)];
}

SP::LagrangianLinearTIDS newBall()
{
  SP::SiconosVector q(new SiconosVector(3));
  SP::SiconosVector v(new SiconosVector(3));
  SP::SiconosMatrix mass(new SimpleMatrix(3,3));
  SP::LagrangianLinearTIDS ball(new LagrangianLinearTIDS(q, v, mass));
  return ball;
};

struct FlatDS
{
  size_t id;
  std::array<double, 3> position;
  std::array<double, 3> velocity;
  std::array<double, 9> mass;
  // padding...
  bool operator < (const FlatDS& ds) const { return id < ds.id; };
};

struct SiconosGraphBench : public benchmark::Fixture
{
  typedef SiconosGraph < int, int,
                         boost::no_property, boost::no_property,
                         boost::no_property > IntGraph;
  typedef SiconosGraph < std::shared_ptr<FlatDS>, int,
                         boost::no_property, boost::no_property,
                         boost::no_property > PFlatDSGraph;
  typedef SiconosGraph < FlatDS, int,
                         boost::no_property, boost::no_property,
                         boost::no_property > FlatDSGraph;

  DynamicalSystemsGraph graph;
  PFlatDSGraph pflatds_graph;
  std::shared_ptr<FlatDSGraph> flatds_graph = std::shared_ptr<FlatDSGraph>(new FlatDSGraph{});
  IntGraph int_graph;

  void SetUp(const ::benchmark::State& state)
  {
    for(unsigned int i=0; i<SIZE; ++i)
    {
      graph.add_vertex(newBall());
      pflatds_graph.add_vertex(std::shared_ptr<FlatDS>(new FlatDS{}));
      flatds_graph->add_vertex(FlatDS{i});
      int_graph.add_vertex(i);
    }

  }

  void TearDown(const ::benchmark::State& state)
  {
  }
};

BENCHMARK_F(SiconosGraphBench, DynamicalSystemAccess)(benchmark::State& state)
{
  DynamicalSystemsGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      SP::DynamicalSystem ds = graph.bundle(*vi);
      SP::SiconosVector velocity = std::static_pointer_cast<LagrangianLinearTIDS>(ds)->velocity();
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
    }
  }
}

BENCHMARK_F(SiconosGraphBench, PFlatDSAccess)(benchmark::State& state)
{
  PFlatDSGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = pflatds_graph.vertices(); vi!=vend; ++vi)
    {
      auto ds = pflatds_graph.bundle(*vi);
      auto velocity = ds->velocity;
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
    }
  }
}

BENCHMARK_F(SiconosGraphBench, FlatDSAccess)(benchmark::State& state)
{
  FlatDSGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = flatds_graph->vertices(); vi!=vend; ++vi)
    {
      auto ds = flatds_graph->bundle(*vi);
      auto velocity = ds.velocity;
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
    }
  }
}

BENCHMARK_F(SiconosGraphBench, IntAccess)(benchmark::State& state)
{
  IntGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = int_graph.vertices(); vi!=vend; ++vi)
    {
      int res = int_graph.bundle(*vi);
      benchmark::DoNotOptimize(res);
    }
  }
}

BENCHMARK_MAIN();
