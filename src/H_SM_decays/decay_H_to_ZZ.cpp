// special case for Higgs -> Z Z
// TODO: implement higher order corrections

template <>
double CLASSNAME::get_partial_width<H,Z,Z>(
   const context_base& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<Z>::type const& indexOut1,
   typename field_indices<Z>::type const& indexOut2
   ) const
{

   const double mHOS = context.physical_mass<H>(indexIn);
   // There might be large differences between mZ from mass block
   // and one from slha input, especially in the decoupling limit
   // so we use the latter one. There might be a problem with
   // models where Z mixes with something else.
   // const double mZOS = context.physical_mass<Z>(indexOut1);
   const double mZOS = qedqcd.displayPoleMZ();
   const double x = Sqr(mZOS/mHOS);
   double res;

   // mH < mZ
   // 4-body decay not implemented for a moment
   if (x > 1.0) {
      const std::string index_as_string = indexIn.size() == 0 ? "" : "(" + std::to_string(indexIn[0]) + ")";
      WARNING("H" + index_as_string + "->ZZ decays: double off-shell decays currently not implemented.");
      return 0.0;
   // mZ < mH < 2*mZ
   // three-body decay
   } else if(4.0*x > 1.0) {
      const double sw2 = Sqr(std::sin(context.model.ThetaW()));
      const double deltaV = 7.0/12.0 - 10.0/9.0*sw2 + 40.0/27.0*Sqr(sw2);

      res = 3./(512.*Power3(Pi)) * 1./mHOS * deltaV * RT(x)/x;

      const auto indices = concatenate(indexOut2, indexOut1, indexIn);
      const auto ghZZ =
         Vertex<Z, Z, H>::evaluate(indices, context).value();

      const double g2 = context.model.get_g2();

      res *= std::norm(ghZZ*g2)/(1-sw2);
   // mH > 2mZ
   // two-body decay
   } else {

      const double flux = 1. / (2 * mHOS);
      // phase space without symmetry factor
      const double ps = 1. / (8. * Pi) * std::sqrt(KallenLambda(1., Sqr(mZOS/mHOS), Sqr(mZOS/mHOS)));

      // phase space symmetry factor
      const double ps_symmetry = 1. / 2.;

      // matrix element squared
      const auto mat_elem = calculate_amplitude<H, Z, Z>(
         context, indexIn, indexOut1, indexOut2);
      const auto mat_elem_sq = mat_elem.square();

      // flux * phase space factor * matrix element squared
      res = flux * ps * ps_symmetry * mat_elem_sq;
   }

   return res;
}
