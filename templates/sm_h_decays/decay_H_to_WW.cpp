// special case for H -> W+ W-
// TODO: implement higher order corrections
template <>
double CLASSNAME::get_partial_width<H, W, conj<W>::type>(
   const context_base& context, typename field_indices<H>::type const& indexIn,
   typename field_indices<conj<W>::type>::type const& indexOut1,
   typename field_indices<W>::type const& indexOut2) const
{

   const double mHOS = context.physical_mass<H>(indexIn);
   const double mWOS = context.physical_mass<W>(indexOut1);
   const double mW = context.mass<W>(indexOut1);
   const double x = Sqr(mWOS / mHOS);
   double res;
   // mH < mW
   // 4-body decay not implemented for a moment
   if (x > 1.0) {
      return 0.0;
   }
   // mW < mH < 2mW
   // three-body decays
   if (4 * x > 1.0) {
      // const auto vd = MODELPARAMETER(vd);
      // const auto vu = MODELPARAMETER(vu);
      // TODO: specify the vev correctly
      const auto vev = 246.0; // sqrt(Sqr(vd) + Sqr(vu));

      res = 3.0 / (128 * std::pow(Pi, 3)) * mHOS / Sqr(vev) * RT(x);
   // mH > 2mZ
   // two-body decay
   } else {

      const double flux = 1. / (2 * mHOS);
      // phase space without symmetry factor
      const double ps = 1. / (8. * Pi) * std::sqrt(KallenLambda(mHOS*mHOS, mWOS*mWOS, mWOS*mWOS))/(mHOS*mHOS);

      // matrix element squared
      const auto mat_elem = calculate_amplitude<H, W, conj<W>::type>(
         context, indexIn, indexOut1, indexOut2);
      const auto mat_elem_sq =mat_elem.square();

      // flux * phase space factor * matrix element squared
      return flux * ps * mat_elem_sq;
   }
   const auto indices = concatenate(indexOut2, indexOut1, indexIn);
   const auto ghWW =
      Vertex<conj<W>::type, W, H>::evaluate(indices, context).value() * std::pow(mWOS/mW, 2);
   return res * std::norm(ghWW);
}
