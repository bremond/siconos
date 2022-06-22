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

#include <boost/pool/object_pool.hpp>

#include "XSiconosGraph.hpp"
#include "XSiconosVector.hpp"
#include "SiconosGraph.hpp"
#include "SiconosKernel.hpp"

static const size_t SIZE = 10000;

SP::LagrangianLinearTIDS newBall()
{
  SP::SiconosVector q(new SiconosVector(3));
  SP::SiconosVector v(new SiconosVector(3));
  SP::SiconosMatrix mass(new SimpleMatrix(3,3));
  SP::LagrangianLinearTIDS ball(new LagrangianLinearTIDS(q, v, mass));
  return ball;
};

struct PseudoDS
{
  SP::SiconosVector _q;
  SP::SiconosVector _v;
  SP::SiconosMatrix _mass;
  PseudoDS(SP::SiconosVector q, SP::SiconosVector v, SP::SiconosMatrix m) :
    _q(q), _v(v), _mass(m) {};

  SP::SiconosVector velocity() { return _v; };
};

TYPEDEF_SPTR(PseudoDS)

std::shared_ptr<PseudoDS> newPseudoBall()
{
  SP::SiconosVector q(new SiconosVector(3));
  SP::SiconosVector v(new SiconosVector(3));
  SP::SiconosMatrix mass(new SimpleMatrix(3,3));
  std::shared_ptr<PseudoDS> ball(new PseudoDS(q, v, mass));
  return ball;
};


struct XPseudoDS
{
  SP::XSiconosVector _q;
  SP::XSiconosVector _v;
  SP::SiconosMatrix _mass;
  XPseudoDS(SP::XSiconosVector q, SP::XSiconosVector v, SP::SiconosMatrix m) :
    _q(q), _v(v), _mass(m) {};

  SP::XSiconosVector velocity() { return _v; };
};

TYPEDEF_SPTR(XPseudoDS)

std::shared_ptr<XPseudoDS> newXPseudoBall()
{
  SP::XSiconosVector q(new XSiconosVector(3));
  SP::XSiconosVector v(new XSiconosVector(3));
  SP::SiconosMatrix mass(new SimpleMatrix(3,3));
  std::shared_ptr<XPseudoDS> ball(new XPseudoDS(q, v, mass));
  return ball;
};

TYPEDEF_SPTR(XDenseVect)

struct XDPseudoDS
{
  XDenseVect _q;
  XDenseVect _v;
  SP::SiconosMatrix _mass;
  XDPseudoDS(XDenseVect q, XDenseVect v, SP::SiconosMatrix m) :
    _q(q), _v(v), _mass(m) {};

  XDenseVect& velocity() { return _v; };

};

TYPEDEF_SPTR(XDPseudoDS)

std::vector<SP::XDenseVect> XDenseVectMem;


SP::XDPseudoDS newXDPseudoBall()
{
  SP::XDenseVect q(new XDenseVect(3));
  SP::XDenseVect v(new XDenseVect(3));
  XDenseVectMem.push_back(q);
  XDenseVectMem.push_back(v);
  SP::SiconosMatrix mass(new SimpleMatrix(3,3));
  std::shared_ptr<XDPseudoDS> ball(new XDPseudoDS(*q, *v, mass));
  return ball;
};

struct newPoolXPseudoBall
{
  boost::object_pool<XSiconosVector> xsiconos_vector_pool;
  boost::object_pool<XPseudoDS> xpseudods_pool;

  SP::XPseudoDS operator() ()
  {
    SP::XSiconosVector q(xsiconos_vector_pool.construct(3));
    SP::XSiconosVector v(xsiconos_vector_pool.construct(3));
    SP::SiconosMatrix mass(new SimpleMatrix(3,3));

    SP::XPseudoDS ball(xpseudods_pool.construct(q,v,mass));
    return ball;
  }

  ~newPoolXPseudoBall()
  {
  }
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
typename G::vertex_t new_vertex(G& graph, size_t i) {};

template<typename G>
struct SiconosGraphBench : public benchmark::Fixture
{
  G graph;
  std::vector<typename G::vertex_t> systems;
  std::vector<double> values;

  void SetUp(const ::benchmark::State& state)
  {
    for(unsigned int i=0; i<SIZE; ++i)
    {
      auto nv = new_vertex(graph, i);
      graph.add_vertex(nv);
      systems.push_back(nv);
      values.push_back(double(i));
    }
  }

  void TearDown(const ::benchmark::State& state)
  {
    graph.clear();
  }
};


#include <iostream>
#include <memory_resource>
#include <chrono>
#include <functional>

class debug_resource : public std::pmr::memory_resource {
public:
  explicit debug_resource(std::string name,
                          std::pmr::memory_resource* up = std::pmr::new_delete_resource())
    : _name{ std::move(name) }, _upstream{ up }
  { }

  void* do_allocate(size_t bytes, size_t alignment) override {
    std::cout << _name << " do_allocate(): " << bytes << '\n';
    void* ret = _upstream->allocate(bytes, alignment);
                  return ret;
  }
  void do_deallocate(void* ptr, size_t bytes, size_t alignment) override {
    std::cout << _name << " do_deallocate(): " << bytes << '\n';
    _upstream->deallocate(ptr, bytes, alignment);
  }
  bool do_is_equal(const std::pmr::memory_resource& other) const noexcept override {
    return this == &other;
  }

private:
  std::string _name;
  std::pmr::memory_resource* _upstream;
};
const auto BufferPrinter = [](std::string_view buf, std::string_view title)
{
  std::cout << title << ":\n";
  for (auto& ch : buf) {
    std::cout << (ch >= ' ' ? ch : '#');
  }
  std::cout << '\n';
};



static void pmrVector(benchmark::State& state)
{
  constexpr size_t BUF_SIZE = 1000000000;
  std::pmr::pool_options options;
  options.max_blocks_per_chunk = 40;
  options.largest_required_pool_block = 640;

  //    alignas(8) std::array<char,BUF_SIZE> buffer; // a small buffer on the stack

  //std::cout <<options.largest_required_pool_block << std::endl;
  char* buffer = new char[BUF_SIZE];
  std::pmr::monotonic_buffer_resource pool{buffer, BUF_SIZE};
  std::pmr::synchronized_pool_resource mem (options,&pool);
  for (auto _ : state)
  {
    {
      std::pmr::vector<char> vec{ &mem };
      for(int i = 0; i != 100000000; ++i)
      {

        vec.emplace_back('a');
      }
      benchmark::DoNotOptimize(vec);


    }
    mem.release();
    pool.release();
  }
  delete[] buffer;
}
static void stdVector(benchmark::State& state)
{
  for (auto _ : state)
  {
    std::vector<char> vec{};
    for(int i = 0; i != SIZE; ++i)
    {

      vec.emplace_back('a');
    }

  }
}

static void boostPoolVector(benchmark::State& state)
{

  for (auto _ : state)
  {
    std::vector<char, boost::pool_allocator<char>> vec{};
    {
      for(int i = 0; i != 100000000; ++i)
      {
        vec.emplace_back('a');
      }
      benchmark::DoNotOptimize(vec);
    }
  }
}

void rm(auto& vec, const auto& i)
{
//  std::remove(vec.begin(), vec.end(), *(vec.begin()
  *(vec.begin()+i) = std::move(vec.back());
  vec.pop_back();
//  std::sort(vec.begin(), vec.end());
//  std::swap(vec.begin()+i, vec.end(), vec[i]);
}

BENCHMARK_TEMPLATE_F(SiconosGraphBench, StdVectorRm2, DynamicalSystemsGraph)(benchmark::State& state)
{
  double k = 0;
  for (auto _ : state)
  {
    values.erase(values.begin()+values.size()/2);
    values.push_back(k++);
    benchmark::DoNotOptimize(values);
  }
}


BENCHMARK_TEMPLATE_F(SiconosGraphBench, StdVectorRm1, DynamicalSystemsGraph)(benchmark::State& state)
{
  double k = 0;
  for (auto _ : state)
  {
    rm(values, values.size()/2);
    values.push_back(k++);
    benchmark::DoNotOptimize(values);
  }
}



BENCHMARK(pmrVector);
BENCHMARK(stdVector);
BENCHMARK(boostPoolVector);

void update_ds_vertices_indices(auto& graph)
{
  DynamicalSystemsGraph::VIterator vi, viend;
  size_t i;
  for (std::tie(vi, viend) = boost::vertices(graph.g), i = 0;
       vi != viend; ++vi, ++i)
  {
    graph.bundle(*vi)->setNumber(i);
    graph.index(*vi) = graph.bundle(*vi)->number();
  }
}


BENCHMARK_TEMPLATE_F(SiconosGraphBench, DynamicalSystemDescriptor1, DynamicalSystemsGraph)(benchmark::State& state)
{
  std::random_device dev;
  std::mt19937 rng(dev());
  DynamicalSystemsGraph::VIterator vi,vend;
  std::uniform_int_distribution<std::mt19937::result_type> random_indice(0, graph.size()-1);
  update_ds_vertices_indices(graph);
  for (auto _ : state)
  {
    SP::DynamicalSystem ds = graph.bundle(boost::vertex(random_indice(rng), graph.g));

    DynamicalSystemsGraph::VDescriptor vd = graph.descriptor(ds);
    benchmark::DoNotOptimize(vd);
  }
}

BENCHMARK_TEMPLATE_F(SiconosGraphBench, DynamicalSystemDescriptor2, DynamicalSystemsGraph)(benchmark::State& state)
{
  std::random_device dev;
  std::mt19937 rng(dev());
  DynamicalSystemsGraph::VIterator vi,vend;
  std::uniform_int_distribution<std::mt19937::result_type> random_indice(0, graph.size()-1);
  update_ds_vertices_indices(graph);
  for (auto _ : state)
  {
    SP::DynamicalSystem ds = graph.bundle(boost::vertex(random_indice(rng), graph.g));
    DynamicalSystemsGraph::VDescriptor vd = boost::vertex(ds->number(), graph.g);
    benchmark::DoNotOptimize(vd);
  }
}


template<>
DynamicalSystemsGraph::vertex_t new_vertex(DynamicalSystemsGraph& g, size_t i)
{
  SP::DynamicalSystem nb = newBall();
  return nb;
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

class PseudoDSGraph :
  public SiconosGraph < std::shared_ptr<PseudoDS>,
                        std::shared_ptr<Interaction>,
                        DynamicalSystemProperties, InteractionProperties,
                        GraphProperties >
{
};



template<>
PseudoDSGraph::vertex_t new_vertex(PseudoDSGraph& g, size_t i)
{
  std::shared_ptr<PseudoDS> nb = newPseudoBall();
  return nb;
}

BENCHMARK_TEMPLATE_F(SiconosGraphBench, PseudoDSAccess, PseudoDSGraph)(benchmark::State& state)
{
  PseudoDSGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      std::shared_ptr<PseudoDS> ds = graph.bundle(*vi);
      SP::SiconosVector velocity = std::static_pointer_cast<PseudoDS>(ds)->velocity();
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
    }
  }
}


class XPseudoDSGraph :
  public SiconosGraph < std::shared_ptr<XPseudoDS>,
                        std::shared_ptr<Interaction>,
                        DynamicalSystemProperties, InteractionProperties,
                        GraphProperties >
{
public:
//  newPoolXPseudoBall xpseudo_ball_maker;
};



template<>
XPseudoDSGraph::vertex_t new_vertex(XPseudoDSGraph& g, size_t i)
{
  std::shared_ptr<XPseudoDS> nb = newXPseudoBall();//g.xpseudo_ball_maker();
  return nb;
}

BENCHMARK_TEMPLATE_F(SiconosGraphBench, XPseudoDSAccess, XPseudoDSGraph)(benchmark::State& state)
{
  XPseudoDSGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      std::shared_ptr<XPseudoDS> ds = graph.bundle(*vi);
      SP::XSiconosVector velocity = std::static_pointer_cast<XPseudoDS>(ds)->velocity();
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
    }
  }
}



class XDPseudoDSGraph :
  public SiconosGraph < std::shared_ptr<XDPseudoDS>,
                        std::shared_ptr<Interaction>,
                        DynamicalSystemProperties, InteractionProperties,
                        GraphProperties >
{
public:
//  newPoolXPseudoBall xpseudo_ball_maker;
};



template<>
XDPseudoDSGraph::vertex_t new_vertex(XDPseudoDSGraph& g, size_t i)
{
  std::shared_ptr<XDPseudoDS> nb = newXDPseudoBall();//g.xpseudo_ball_maker();
  return nb;
}

BENCHMARK_TEMPLATE_F(SiconosGraphBench, XDPseudoDSAccess, XDPseudoDSGraph)(benchmark::State& state)
{
  XDPseudoDSGraph::VIterator vi,vend;

  for (auto _ : state)
  {
    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      std::shared_ptr<XDPseudoDS> ds = graph.bundle(*vi);
      XDenseVect& velocity = std::static_pointer_cast<XDPseudoDS>(ds)->velocity();
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
    }
  }
}




class XDynamicalSystemsGraph :
  public XSiconosGraph < std::shared_ptr<DynamicalSystem>,
                         std::shared_ptr<Interaction>,
                         DynamicalSystemProperties, InteractionProperties,
                         GraphProperties >
{
};

template<>
XDynamicalSystemsGraph::vertex_t new_vertex(XDynamicalSystemsGraph& g, size_t i)
{
  SP::DynamicalSystem nb = newBall();
  return nb;
}


BENCHMARK_TEMPLATE_F(SiconosGraphBench, XDynamicalSystemAccess, XDynamicalSystemsGraph)(benchmark::State& state)
{
  XDynamicalSystemsGraph::VIterator vi,vend;

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
PFlatDSGraph::vertex_t new_vertex(PFlatDSGraph& g, size_t i)
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
FlatDSGraph::vertex_t new_vertex(FlatDSGraph& g, size_t i)
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
XFlatDSGraph::vertex_t new_vertex(XFlatDSGraph& g, size_t i)
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
IntGraph::vertex_t new_vertex(IntGraph& g, size_t i)
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
XIntGraph::vertex_t new_vertex(XIntGraph& g, size_t i)
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

BENCHMARK_TEMPLATE_F(SiconosGraphBench, DynamicalSystemAddRemove, DynamicalSystemsGraph)(benchmark::State& state)
{
  std::random_device dev;
  std::mt19937 rng(dev());

  DynamicalSystemsGraph::VIterator vi,vend;

  for (auto _ : state)
  {

    for(size_t i=0; i<graph.size()/10; ++i)
    {
      std::uniform_int_distribution<std::mt19937::result_type> random_indice(0, graph.size()-1);
      graph.update_vertices_indices();

      SP::DynamicalSystem ds = graph.bundle(boost::vertex(random_indice(rng), graph.g));

//      std::cout << "remove? " << ds->number() <<std::endl;
      if (graph.is_vertex(ds))
      {

        graph.remove_vertex(ds);

        SP::DynamicalSystem nds = newBall();
        graph.add_vertex(ds);
      }
    }

    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      SP::DynamicalSystem ds = graph.bundle(*vi);
      SP::SiconosVector velocity = std::static_pointer_cast<LagrangianLinearTIDS>(ds)->velocity();
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
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



BENCHMARK_TEMPLATE_F(SiconosGraphBench, XDynamicalSystemAddRemove, XDynamicalSystemsGraph)(benchmark::State& state)
{
  std::random_device dev;
  std::mt19937 rng(dev());

  XDynamicalSystemsGraph::VIterator vi,vend;

  for (auto _ : state)
  {

    for(size_t i=0; i<graph.size()/10; ++i)
    {
      std::uniform_int_distribution<std::mt19937::result_type> random_indice(0,graph.size()-1);
      //graph.update_vertices_indices();
      SP::DynamicalSystem ds = graph.bundle(boost::vertex(random_indice(rng), graph.g));

//      std::cout << "xremove? " << ds->number() <<std::endl;

//      if (graph.is_vertex(ds))
      {
        graph.remove_vertex(ds);

        SP::DynamicalSystem nds = newBall();
        graph.add_vertex(nds);
      }
    }

    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      SP::DynamicalSystem ds = graph.bundle(*vi);
      SP::SiconosVector velocity = std::static_pointer_cast<LagrangianLinearTIDS>(ds)->velocity();
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
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


BENCHMARK_TEMPLATE_F(SiconosGraphBench, PseudoDSAddRemove, PseudoDSGraph)(benchmark::State& state)
{
  std::random_device dev;
  std::mt19937 rng(dev());

  PseudoDSGraph::VIterator vi,vend;

  for (auto _ : state)
  {

//    for(size_t i=0; i<graph.size()/10; ++i)
    {
      std::uniform_int_distribution<std::mt19937::result_type> random_indice(0, graph.size()-1);
      graph.update_vertices_indices();

      SP::PseudoDS ds = graph.bundle(boost::vertex(random_indice(rng), graph.g));

      if (graph.is_vertex(ds))
      {
        graph.remove_vertex(ds);

        SP::PseudoDS nds = newPseudoBall();
        graph.add_vertex(nds);
      }
    }

    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      SP::PseudoDS ds = graph.bundle(*vi);
      SP::SiconosVector velocity = std::static_pointer_cast<PseudoDS>(ds)->velocity();
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
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

BENCHMARK_TEMPLATE_F(SiconosGraphBench, XPseudoDSAddRemove, XPseudoDSGraph)(benchmark::State& state)
{
  std::random_device dev;
  std::mt19937 rng(dev());

  XPseudoDSGraph::VIterator vi,vend;

  for (auto _ : state)
  {

//    for(size_t i=0; i<graph.size()/10; ++i)
    {
      std::uniform_int_distribution<std::mt19937::result_type> random_indice(0, graph.size()-1);
      graph.update_vertices_indices();

      SP::XPseudoDS ds = graph.bundle(boost::vertex(random_indice(rng), graph.g));

//      std::cout << "remove? " << ds->number() <<std::endl;
      if (graph.is_vertex(ds))
      {

        graph.remove_vertex(ds);

        SP::XPseudoDS nds = newXPseudoBall();
        graph.add_vertex(nds);
      }
    }

    for(std::tie(vi,vend) = graph.vertices(); vi!=vend; ++vi)
    {
      SP::XPseudoDS ds = graph.bundle(*vi);
      SP::XSiconosVector velocity = std::static_pointer_cast<XPseudoDS>(ds)->velocity();
      benchmark::DoNotOptimize(ds);
      benchmark::DoNotOptimize(velocity);
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
