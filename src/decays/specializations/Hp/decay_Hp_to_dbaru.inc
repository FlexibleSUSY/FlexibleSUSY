template <>
double CLASSNAME::get_partial_width<Hp, bar<d>::type, u>(
   const context_base& context,
   typename field_indices<Hp>::type const& indexIn,
   typename field_indices<bar<d>>::type const& indexOut1,
   typename field_indices<u>::type const& indexOut2)
{
   const double mHp = context.physical_mass<Hp>(indexIn);
   const double md = context.physical_mass<d>(indexOut1);
   const double mu = context.physical_mass<u>(indexOut2);

   // phase space without symmetry factor
   const double ps = 1./(8.*Pi) * std::sqrt(KallenLambda(1., Sqr(mu/mHp), Sqr(md/mHp)));

   static constexpr double color_factor = 3;

   // matrix element squared
   const auto mat_elem_sq = amplitude_squared<Hp, bar<d>::type, u>(
      context, indexIn, indexOut1, indexOut2);

   // flux * phase space factor * symmetry factor * color factor * |matrix element|^2
   double result = 0.5/mHp * ps * color_factor * mat_elem_sq;

   // higher order corrections
   if (flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) {
      // H+ -> t bbar
      if(indexOut1.at(0) == 2 && indexOut2.at(0) == 2) {
         const double z = Sqr(mu/mHp);
         static constexpr double CF = 4./3.;
         // @todo: missing MSbar -> OS conversion factor?
         const double corr = 1. + CF*get_alphas(context)/Pi*Delta_Hp(z);
         // Delta_Hp has logarithmic singularity for z->0
         if (corr > 0) {
            result *= corr;
         }
      }
   }

   return result;
}

