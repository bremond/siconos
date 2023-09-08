#pragma once
#include <petscmat.h>

#include "siconos/storage/traits/traits.hpp"

namespace siconos {
namespace petsc {

struct base_mat {
  using petsc_base_mat_t = void;
};

template <typename T>
struct mat : base_mat {
  using value_type = T;
  using petsc_mat_t = void;

  static constexpr auto value_nb_rows = traits::ncols(T{});
  static constexpr auto value_nb_cols = traits::nrows(T{});

  // Constructors and destructor
  mat() : mat(nullptr) {}
  mat& operator=(mat&&) = delete;
  mat(const mat&) = delete;
  mat(MPI_Comm comm, PetscInt m, PetscInt n) : comm(comm)
  {
    create(comm, m, n);
  }

  mat(mat&& other) { swap(other); }
  ~mat() { destroy(); }

  // Assignment operator
  mat& operator=(mat other)
  {
    swap(other);
    return *this;
  }

  // Member functions
  void create(MPI_Comm comm, PetscInt m, PetscInt n)
  {
    destroy();
    MatCreate(comm, &mat);
    MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, m, n);
    MatSetUp(mat);
  }

  void destroy()
  {
    if (mat) {
      MatDestroy(&mat);
      mat = nullptr;
    }
  }

  void swap(PetscMatrix& other) { std::swap(mat, other.mat); }

  void set(PetscInt i, PetscInt j, PetscScalar value)
  {
    MatSetValue(mat, i, j, value, INSERT_VALUES);
  }

  void assemble()
  {
    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
  }

  void view() const { MatView(mat, PETSC_VIEWER_STDOUT_WORLD); }

  Mat mat;
  MPI_Comm comm;
};
void resize(match::any_mat auto& amat)
{
  m.destroy();
    m.create(
}
}  // namespace petsc
}  // namespace siconos
