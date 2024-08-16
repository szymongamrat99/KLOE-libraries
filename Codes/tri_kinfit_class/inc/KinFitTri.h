#ifndef KINFITTRI_H
#define KINFITTRI_H

#include "KinFit.h"
#include "ConstraintsTri.h"

namespace KLOE
{

  class KinFitTri : public KinFit
  {
  private:
    ConstraintsTri constraints;
  public:
    KinFitTri(Int_t M_init);
    ~KinFitTri();
  };
  
  KinFitTri::KinFitTri(Int_t M_init) : constraints()
  {
    
    _chosen.resize(_M);
  }
  
  KinFitTri::~KinFitTri()
  {
  }
  
  
} // namespace KLOE


#endif // !KINFITTRI_H