#include <boost/numeric/ublas/vector.hpp>
#include <boost/pool/object_pool.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <SiconosVector.hpp>
#include <memory_resource>

#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/std/vector.hpp>


template<size_t MBPC, size_t LRPB, size_t S>
struct data_container
{
  std::pmr::pool_options options = {
    .max_blocks_per_chunk = MBPC,
    .largest_required_pool_block = LRPB};
  char *buffer = new char[S];
  std::pmr::monotonic_buffer_resource pool = {buffer, S};
  std::pmr::synchronized_pool_resource mem = std::pmr::synchronized_pool_resource(options,&pool);

  ~data_container()
  {
    mem.release();
    pool.release();
    delete[] buffer;
  }
};


struct X
{
  std::vector<double, boost::pool_allocator<double>> x;
  std::vector<double, boost::pool_allocator<double>> v;
  double miamiam;
};

struct Y
{
  std::vector<double, boost::pool_allocator<double>> x;
  std::vector<double, boost::pool_allocator<double>> v;
  int ok;
  double okok;
  double z;
};



struct Xp
{
  std::vector<double> x;
  std::vector<double> v;
  double miamiam;

};


struct Yp
{
  std::vector<double> x;
  std::vector<double> v;
  int ok;
  double okok;
  double z;
};

int main()
{
  std::vector<X> Xs;
  std::vector<Y> Ys;

  std::vector<Xp> Xps;
  std::vector<Yp> Yps;

  for(size_t i = 0; i<10000; ++i)
  {
    Xs.push_back(X{{1.,2.,3.},{4.,5.,6.}}); Ys.push_back(Y{{1.,2.,3.},{4.,5.,6.}});
    Xps.push_back(Xp{{1.,2.,3.},{4.,5.,6.}}); Yps.push_back(Yp{{1.,2.,3.},{4.,5.,6.}});
  }

  std::cout << &Xs[99].x[0] - &Xs[98].x[0] << std::endl;

  std::cout << &Xps[99].x[0] - &Xps[98].x[0] << std::endl;

  data_container<40, 640, 1000> dc;

  std::pmr::vector<double> v { &dc.mem };

  std::vector<double, std::pmr::polymorphic_allocator<double>> w;

  w = v;

  typedef ublas::vector<double, std::vector<double,std::pmr::polymorphic_allocator<double> >> DenseVect;

  DenseVect uv;

}
