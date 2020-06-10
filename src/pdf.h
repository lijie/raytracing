#ifndef __PDF_H__
#define __PDF_H__

#include "orthonormal_basis.h"
#include "rand.h"
#include "vec3.h"

class Pdf {
 public:
  virtual double Value(const Vec3& dir) const = 0;
  virtual Vec3 Generate() const = 0;
};

class CosinePdf : public Pdf {
 public:
  CosinePdf(const Vec3& n) { ob_.BuildFromW(n); }
  virtual double Value(const Vec3& dir) const {
    auto cosine = dot(unit_vector(dir), ob_.w());
    return (cosine < 0) ? 0 : cosine / M_PI;
  }

  virtual Vec3 Generate() const { return unit_vector(ob_.Local(RandomCosineDirection())); }

 private:
  OrthonormalBasis ob_;
};

#endif  // __PDF_H__
