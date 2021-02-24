#include <IMP/foo/MyRestraint.h>
#include <IMP/core/XYZ.h>

IMPFOO_BEGIN_NAMESPACE

namespace {
  class ScoreCalculator {
    core::XYZ xyz_;
    double k_, k2_;
  public:
    ScoreCalculator(core::XYZ xyz, double k) : xyz_(xyz), k_(k) {}

    double get_score() { return .5 * k2_ * square(xyz_.get_z()); }
  };
}

MyRestraint::MyRestraint(Model *m, ParticleIndex p, double k)
    : Restraint(m, "MyRestraint%1%"), p_(p), k_(k) {}

void MyRestraint::do_add_score_and_derivatives(ScoreAccumulator sa) const {
  core::XYZ d(get_model(), p_);
  ScoreCalculator *calc = new ScoreCalculator(d, k_);
  double score = calc->get_score();
  if (sa.get_derivative_accumulator()) {
    double deriv = k_ * d.get_z();
    d.add_to_derivative(2, deriv, *sa.get_derivative_accumulator());
  }
  sa.add_score(score);
}

ModelObjectsTemp MyRestraint::do_get_inputs() const {
  return ModelObjectsTemp(1, get_model()->get_particle(p_));
}

IMPFOO_END_NAMESPACE
