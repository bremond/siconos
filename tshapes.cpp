#include <benchmark/benchmark.h>

#include <array>
#include <variant>
#include <random>
#include <functional>
#include <algorithm>
#include <iostream>
#include <memory>


static const size_t SIZE = 100000;

typedef float scal;

struct VShape
{
  virtual const scal& param(size_t) const noexcept = 0;
  virtual ~VShape(){};
  virtual const size_t size() const noexcept = 0;
};

template<size_t N>
struct VPoly : VShape
{
  size_t _size = N;

  std::shared_ptr<std::array<scal, N>> _data;

  VPoly<N>() : _data(new std::array<scal, N>()) {};

  const scal& param(size_t i) const noexcept override
  {
    return (*_data)[i];
  }

  const size_t size() const noexcept override { return _size; };
};

template<size_t N>
struct VSPoly : VShape
{
  size_t _size = N;

  std::array<scal, N> _data;

  const scal& param(size_t i) const noexcept override
  {
    return _data[i];
  }

  const size_t size() const noexcept override { return _size; };
};


template<typename T, typename DT>
struct Shape {

  std::shared_ptr<DT> _data;

  Shape() : _data(new DT()) {};

  ~Shape() { };

  const scal& param(size_t i) const noexcept
  {
    return (*_data)[i];
  }

};

struct Disk : Shape<Disk, std::array<scal, 1> >
{
  auto radius()
  {
    return (*(this->_data))[0];
  }
};

template <size_t N>
struct Poly : Shape<Poly<N>, std::array<scal, N> >
{
  size_t size=N;
  auto points()
  {
    return *(this->_data);
  }
};


std::random_device dev;
std::mt19937 rng(dev());
std::uniform_int_distribution<std::mt19937::result_type> random_pick(0,3);
std::uniform_int_distribution<std::mt19937::result_type> random_indice(0,SIZE-1);

template <typename T, size_t N>
std::array<T, N> random_vector_stack(std::function<T(decltype(random_pick(rng)))> func) {
  std::array<T, N> a;
  std::generate(a.begin(), a.end(), [&] {
    return func(random_pick(rng));
  });
  return a;
}

template <typename T, size_t N>
auto random_vector_heap(std::function<T(decltype(random_pick(rng)))> func) {
  std::shared_ptr<std::vector<T> > a(new std::vector<T>());
  a->resize(N);
  std::generate(a->begin(), a->end(), [&] {
    return func(random_pick(rng));
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

struct MyFixture : public benchmark::Fixture
{
  typedef std::variant<Poly<1>,Poly<2>,Poly<3>,Poly<4> > vtype;

  std::shared_ptr<std::vector<vtype>> shapes_with_variant =
    random_vector_heap<vtype, SIZE>([&] (auto r) -> vtype
    {
      switch(r) {
      case 0: return Poly<1>{};
      case 1: return Poly<2>{};
      case 2: return Poly<3>{};
      case 3: return Poly<4>{};
      default: exit(1);
      }
    });;

  std::shared_ptr<std::vector<VShape*>> shapes_with_inheritance =
    random_vector_heap<VShape*, SIZE>([&] (auto r) -> VShape*
    {
      switch(r) {
      case 0: return new VPoly<1>{};
      case 1: return new VPoly<2>{};
      case 2: return new VPoly<3>{};
      case 3: return new VPoly<4>{};
      default: exit(1);
      }
    });


  void SetUp(const ::benchmark::State& state)
  {
  }

  void TearDown(const ::benchmark::State& state)
  {
  }
};

BENCHMARK_F(MyFixture, Inheritance)(benchmark::State& state)
{
  VShape* shape = nullptr;

  scal res;

  for (auto _ : state)
  {
    shape = random_pitem(*(this->shapes_with_inheritance));

    res = shape->param(0);
//    std::cout << shape->size() << std::endl;

    benchmark::DoNotOptimize(shape);
    benchmark::DoNotOptimize(res);
  }
}


BENCHMARK_F(MyFixture,Variant)(benchmark::State& state)
{

  const vtype* shape = nullptr;

  scal res;

  for (auto _ : state)
  {
    shape = random_item(*(this->shapes_with_variant));

    switch(shape->index())
    {
    case 0:
      res = std::get<Poly<1>>(*shape).param(0);
//      std::cout << std::get<Poly<1>>(*shape).size << std::endl;
      break;
    case 1:
      res = std::get<Poly<2>>(*shape).param(0);
//      std::cout << std::get<Poly<2>>(*shape).size << std::endl;
      break;
    case 2:
      res = std::get<Poly<3>>(*shape).param(0);
//      std::cout << std::get<Poly<3>>(*shape).size << std::endl;
      break;
    case 3:
      res = std::get<Poly<4>>(*shape).param(0);
//      std::cout << std::get<Poly<4>>(*shape).size << std::endl;
      break;
    }

    benchmark::DoNotOptimize(shape);
    benchmark::DoNotOptimize(res);

  }
}



BENCHMARK_MAIN();
