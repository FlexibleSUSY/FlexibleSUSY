// special case for H -> W+ W-
// TODO: implement higher order corrections
template <>
double CLASSNAME::get_partial_width<H, conj<W>::type, W>(
   const context_base& context, typename field_indices<H>::type const& indexIn,
   typename field_indices<conj<W>::type>::type const& indexOut1,
   typename field_indices<W>::type const& indexOut2) const
{

   const double mHOS = context.physical_mass<H>(indexIn);
   const double mWOS = context.physical_mass<W>(indexOut1);
   const double mW = context.mass<W>(indexOut1);
   const double x = Sqr(mWOS / mHOS);
   double res;

   // 4-body decay for mH < mW not implemented for a moment
   if (x > 1.0) {
      res = 0.0;
   }

   // three-body decays form mW < mH < w mW
   else if (4 * x > 1.0) {
      const auto vev = context.model.VEV();

      res = 3.0 / (128 * Power3(Pi)) * mHOS / Sqr(vev) * RT(x);

      const auto indices = concatenate(indexOut2, indexOut1, indexIn);
      const auto ghWW =
         Vertex<conj<W>::type, W, H>::evaluate(indices, context).value() * Sqr(mWOS/mW);
      res *= std::norm(ghWW);

   // two-body decay for mH > 2 mW
   } else {

      const double flux = 1. / (2 * mHOS);
      // phase space without symmetry factor
      const double ps = 1. / (8. * Pi) * std::sqrt(KallenLambda(mHOS*mHOS, mWOS*mWOS, mWOS*mWOS))/(mHOS*mHOS);

      // matrix element squared
      const auto mat_elem = calculate_amplitude<H, conj<W>::type, W>(
         context, indexIn, indexOut1, indexOut2);
      const auto mat_elem_sq =mat_elem.square();

      // flux * phase space factor * matrix element squared
      res = flux * ps * mat_elem_sq;
   }

   if (res < 0) {
      throw std::runtime_error("Error in H->WW. Partial width < 0.");
   } else {
      return res;
   }
}
