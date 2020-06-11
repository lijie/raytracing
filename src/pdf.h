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

  virtual Vec3 Generate() const {
    return unit_vector(ob_.Local(RandomCosineDirection()));
  }

 private:
  OrthonormalBasis ob_;
};

class Hitable;
class ShapePdf : public Pdf {
 public:
  ShapePdf(Hitable* shape, const Vec3& origin)
      : shape_(shape), origin_(origin) {}

  double Value(const Vec3& dir) const override;
  Vec3 Generate() const override;

 private:
  Hitable* shape_;
  Vec3 origin_;
};

class MixturePdf : public Pdf {
 public:
  MixturePdf(Pdf* p1, Pdf* p2, double mix) : p1_(p1), p2_(p2), mix_(mix) {}

  double Value(const Vec3& dir) const override;
  Vec3 Generate() const override;

 private:
  Pdf* p1_;
  Pdf* p2_;
  double mix_;
};

#endif  // __PDF_H__
