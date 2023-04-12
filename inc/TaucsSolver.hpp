
//=============================================================================
//
//  CLASS TaucsSolver: 
//  direct solver for symmetric positive definite sparse systems
//
//=============================================================================


#ifndef TAUCS_SOLVER_HH
#define TAUCS_SOLVER_HH


//== INCLUDES =================================================================

#include <OpenMesh/Core/System/config.hh> // to define OM_INCLUDE_TEMPLATES...
#include <vector>
#include <set>
#include <taucs.h>


//== CLASS DEFINITION =========================================================


class TaucsSolver
{
public:
   
  TaucsSolver();
  ~TaucsSolver();

  void begin_row();
  void add_value(int _i, double _val);
  void end_row();

  bool factorize(bool _use_supernodal=true);

  bool solve(std::vector<double>& _b, std::vector<double>& _x);

  template <class Vec>
  bool vec_solve(std::vector<Vec>& _b, std::vector<Vec>& _x);
  

private:

  void delete_matrices();


private:

  taucs_ccs_matrix           A, *PAP, *L;
  void                       *SL;
  std::vector<double>        values;
  std::vector<int>           colptr;
  std::vector<int>           rowind;
  int                        n_rows;
  int                        *perm, *invperm;
  bool                       supernodal_;
};


//=============================================================================
#if defined(OM_INCLUDE_TEMPLATES) && !defined(TAUCS_SOLVER_C)
#  define MB_TAUCS_SOLVER_TEMPLATES
#  include "TaucsSolver.cpp"
#endif
//=============================================================================
#endif // TAUCS_SOLVER_HH defined
//=============================================================================
