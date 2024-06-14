template <>
double CLASSNAME::get_partial_width<Hp, Higgs, WpBoson>(
   const context_base& context,
   typename field_indices<Hp>::type const& indexIn,
   typename field_indices<Higgs>::type const& indexOut1,
   typename field_indices<WpBoson>::type const& indexOut2) const
{

   const double mHpOS = context.physical_mass<Hp>(indexIn);
   const double mHOS = context.physical_mass<Higgs>(indexOut1);
   const double mWOS = context.physical_mass<WpBoson>(indexOut2);

   if (mHpOS >= mHOS + mWOS) {
      const double xH = Sqr(mHOS/mHpOS);
      const double xW = Sqr(mWOS/mHpOS);

      const auto indices = concatenate(indexOut1, indexIn, indexOut2);
      const auto gHpHWm = Vertex<Higgs, Hp, WmBoson>::evaluate(indices, context).value(0, 1);

      return 1./(16.*Pi) * Cube(std::sqrt(KallenLambda(1., xH, xW))*mHpOS)/Sqr(mWOS) * std::norm(gHpHWm);
   }
   else if (mHpOS - mHOS > 6) {

      const double omega   = Sqr(mHOS /mWOS);
      const double omegapm = Sqr(mHpOS/mWOS);

      return 0.;
   }
   else {
      return 0.;
   }
}
