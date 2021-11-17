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

#include "XSiconosGraph.hpp"
#include "SiconosKernel.hpp"


static const size_t SIZE = 100000;

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
  bool operator < (const FlatDS& ds) const { return id < ds.id; };
  bool operator == (const FlatDS& ds) const { return id == ds.id; };
};

template<typename G>
typename G::vertex_t new_vertex(const G& graph, size_t i) {};

template<typename G>
struct SiconosGraphBench : public benchmark::Fixture
{
  G graph;

  void SetUp(const ::benchmark::State& state)
  {
    for(unsigned int i=0; i<SIZE; ++i)
    {
      graph.add_vertex(new_vertex(graph, i));
    }
  }

  void TearDown(const ::benchmark::State& state)
  {
    graph.clear();
  }
};


// struct SiconosGraphBench : public benchmark::Fixture
// {
//   typedef SiconosGraph < size_t, size_t,
//                          FlatDS, boost::no_property,
//                          boost::no_property > IntGraph;

//   typedef XSiconosGraph < size_t, size_t,
//                          FlatDS, boost::no_property,
//                          boost::no_property > XIntGraph;
//   typedef SiconosGraph < std::shared_ptr<FlatDS>, int,
//                          boost::no_property, boost::no_property,
//                          boost::no_property > PFlatDSGraph;
//   typedef SiconosGraph < FlatDS, int,
//                          boost::no_property, boost::no_property,
//                          boost::no_property > FlatDSGraph;

//   typedef XSiconosGraph < FlatDS, int,
//                           boost::no_property, boost::no_property,
//                           boost::no_property > XFlatDSGraph;

//   DynamicalSystemsGraph graph;
//   PFlatDSGraph pflatds_graph;
//   std::shared_ptr<FlatDSGraph> flatds_graph = std::shared_ptr<FlatDSGraph>(new FlatDSGraph{});
//   std::shared_ptr<XFlatDSGraph> xflatds_graph = std::shared_ptr<XFlatDSGraph>(new XFlatDSGraph{});
//   IntGraph int_graph;
//   XIntGraph xint_graph;

//   void SetUp(const ::benchmark::State& state)
//   {
//     for(unsigned int i=0; i<SIZE; ++i)
//     {
//       graph.add_vertex(newBall());
//       pflatds_graph.add_vertex(std::shared_ptr<FlatDS>(new FlatDS{}));
//       flatds_graph->add_vertex(FlatDS{i});
//       xflatds_graph->add_vertex(FlatDS{i});
//       int_graph.add_vertex(i);
//       xint_graph.add_vertex(i);
//     }

//   }

//   void TearDown(const ::benchmark::State& state)
//   {
//   }
// };

template<>
DynamicalSystemsGraph::vertex_t new_vertex(const DynamicalSystemsGraph& g, size_t i)
{
  return newBall();
}

BENCHMARK_TEMPLATE_F(SiconosGraphBench, DynamicalSystemAccess, DynamicalSystemsGraph)(benchmark::State& state)
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

typedef SiconosGraph < std::shared_ptr<FlatDS>, int,
                       boost::no_property, boost::no_property,
                       boost::no_property > PFlatDSGraph;

template<>
PFlatDSGraph::vertex_t new_vertex(const PFlatDSGraph& g, size_t i)
{
  return std::shared_ptr<FlatDS>(new FlatDS{i});
};

BENCHMARK_TEMPLATE_F(SiconosGraphBench, PFlatDSAccess, PFlatDSGraph)(benchmark::State& state)
{
  PFlatDSGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      auto ds = graph.bundle(*vi);
      auto velocity = ds->velocity;
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
    }
  }
}

typedef SiconosGraph < FlatDS, int,
                       boost::no_property, boost::no_property,
                       boost::no_property > FlatDSGraph;

template<>
FlatDSGraph::vertex_t new_vertex(const FlatDSGraph& g, size_t i)
{
  return FlatDS{i};
};

BENCHMARK_TEMPLATE_F(SiconosGraphBench, FlatDSAccess, FlatDSGraph)(benchmark::State& state)
{
  FlatDSGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      auto ds = graph.bundle(*vi);
      auto velocity = ds.velocity;
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
    }
  }
}

typedef XSiconosGraph < FlatDS, int,
                        boost::no_property, boost::no_property,
                        boost::no_property > XFlatDSGraph;

template<>
XFlatDSGraph::vertex_t new_vertex(const XFlatDSGraph& g, size_t i)
{
  return FlatDS{i};
};

BENCHMARK_TEMPLATE_F(SiconosGraphBench, XFlatDSAccess, XFlatDSGraph)(benchmark::State& state)
{
  XFlatDSGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      auto ds = graph.bundle(*vi);
      auto velocity = ds.velocity;
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
    }
  }
}

typedef SiconosGraph < size_t, size_t,
                          FlatDS, boost::no_property,
                          boost::no_property > IntGraph;

template<>
IntGraph::vertex_t new_vertex(const IntGraph& g, size_t i)
{
  return i;
};


BENCHMARK_TEMPLATE_F(SiconosGraphBench, IntAccess, IntGraph)(benchmark::State& state)
{
  IntGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      int res = graph.bundle(*vi);
      benchmark::DoNotOptimize(res);
    }
  }
}


BENCHMARK_TEMPLATE_F(SiconosGraphBench, IntPropertyAccess, IntGraph)(benchmark::State& state)
{
  IntGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      auto ds = graph.properties(*vi);
      auto velocity = ds.velocity;
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
    }
  }
}

typedef XSiconosGraph < size_t, size_t,
                       FlatDS, boost::no_property,
                       boost::no_property > XIntGraph;

template<>
XIntGraph::vertex_t new_vertex(const XIntGraph& g, size_t i)
{
  return i;
};


BENCHMARK_TEMPLATE_F(SiconosGraphBench, XIntAccess, XIntGraph)(benchmark::State& state)
{
  XIntGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      int res = graph.bundle(*vi);
      benchmark::DoNotOptimize(res);
    }
  }
}

BENCHMARK_TEMPLATE_F(SiconosGraphBench, XIntPropertyAccess, XIntGraph)(benchmark::State& state)
{
  XIntGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      auto ds = graph.properties(*vi);
      auto velocity = ds.velocity;
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
    }
  }
}

BENCHMARK_TEMPLATE_F(SiconosGraphBench, IntAddRemove, IntGraph)(benchmark::State& state)
{
  std::random_device dev;
  std::mt19937 rng(dev());

  IntGraph::VIterator vi,vend;

  for (auto _ : state)
  {
//    for(size_t i=0; i<graph.size()/10; ++i)
    {
      std::uniform_int_distribution<std::mt19937::result_type> random_indice(0,graph.size()-1);
      graph.update_vertices_indices();

      size_t indx = random_indice(rng);
      if (graph.is_vertex(indx))
      {
        graph.remove_vertex(indx);
        graph.add_vertex(new_vertex(graph, indx));
      }
    }

    // for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    // {
    //   auto ds = graph.properties(*vi);
    //   auto velocity = ds.velocity;
    //   benchmark::DoNotOptimize(ds);
    //   benchmark::DoNotOptimize(velocity);
    // }
  }
}

BENCHMARK_TEMPLATE_F(SiconosGraphBench, XIntAddRemove, XIntGraph)(benchmark::State& state)
{
  std::random_device dev;
  std::mt19937 rng(dev());

  XIntGraph::VIterator vi,vend;

  for (auto _ : state)
  {

//    for(size_t i=0; i<graph.size()/10; ++i)
    {
      std::uniform_int_distribution<std::mt19937::result_type> random_indice(0,graph.size()-1);
//    xint_graph.update_vertices_indices();

      size_t indx = random_indice(rng);
      if (graph.is_vertex(indx))
      {
        graph.remove_vertex(indx);
        graph.add_vertex(new_vertex(graph, indx));
      }
    }

    // for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    // {
    //   auto ds = graph.properties(*vi);
    //   auto velocity = ds.velocity;
    //   benchmark::DoNotOptimize(ds);
    //   benchmark::DoNotOptimize(velocity);
    // }
  }
}

BENCHMARK_MAIN();
