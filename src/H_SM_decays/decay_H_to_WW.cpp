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
   const double x = Sqr(mWOS/mHOS);
   double res;

   // 4-body decay for mH < mW not implemented for a moment
   if (x > 1.0) {
      const std::string index_as_string = indexIn.size() == 0 ? "" : "(" + std::to_string(indexIn[0]) + ")";
      WARNING("H" + index_as_string + "->W+W- decays: double off-shell decays currently not implemented.");
      res = 0.0;
   }
   // three-body decays form mW < mH < 2*mW
   else if (4 * x > 1.0) {

      res = 1./(768.*Power3(Pi)) * 1./mHOS * RT(x)/x;

      const auto indices = concatenate(indexOut2, indexOut1, indexIn);
      const auto ghWW =
         Vertex<conj<W>::type, W, H>::evaluate(indices, context).value();

      // absolute value of baru d W+ vertex (no CKM and no PL projector)
      const double g2 = context.model.get_g2();
      const double gWud = g2*M_SQRT1_2;

      res *= std::norm(ghWW*gWud);

      // multiply by number of final states
      constexpr double NLF = 3;  // number of lepton flavours
      constexpr double Nc = 3;   // number of colors
      constexpr double NQF = 2;  // number of quark flavours
      res *= NLF + Nc*NQF;

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

   return res;
}
