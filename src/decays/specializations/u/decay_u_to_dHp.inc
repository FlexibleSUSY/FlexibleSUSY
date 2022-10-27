template <>
double CLASSNAME::get_partial_width<u, d, Hp>(
   const context_base& context,
   typename field_indices<u>::type const& indexIn,
   typename field_indices<d>::type const& indexOut1,
   typename field_indices<Hp>::type const& indexOut2)
{
   const double mu = context.physical_mass<u>(indexIn);
   const double md = context.physical_mass<d>(indexOut1);
   const double mHp = context.physical_mass<Hp>(indexOut2);


   const double ps = 1./(8.*Pi) * std::sqrt(KallenLambda(1., Sqr(md/mu), Sqr(mHp/mu)));

   // matrix element squared
   const auto mat_elem_sq = amplitude_squared<u, d, Hp>(
      context, indexIn, indexOut1, indexOut2);

   // flux * phase space factor * symmetry factor * color factor * |matrix element|^2
   double result = 0.5/mu * ps * mat_elem_sq;

   // higher order corrections
   if (flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) {
      // t -> b H^+
      if(indexIn.at(0) == 2 && indexOut1.at(0) == 2) {
         const double z = Sqr(mHp/mu);
         static constexpr double CF = 4./3.;
         // @todo: missing MSbar -> OS conversion term?
         result *= 1. + CF*get_alphas(context)/Pi*Delta_t(z);
      }
   }

   return result;
}

