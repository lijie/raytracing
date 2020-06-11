#include "pdf.h"

#include "hitable.h"

double ShapePdf::Value(const Vec3& dir) const {
  return shape_->SurfaceRandomPdf(origin_, dir);
}

Vec3 ShapePdf::Generate() const { return shape_->SurfaceRandomPoint(origin_); }

double MixturePdf::Value(const Vec3& dir) const {
  return p1_->Value(dir) * mix_ + p2_->Value(dir) * (1 - mix_);
}

Vec3 MixturePdf::Generate() const {
  return (RandDouble() < mix_) ? p1_->Generate() : p2_->Generate();
}
