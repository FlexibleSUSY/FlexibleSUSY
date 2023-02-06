template <>
double CLASSNAME::get_partial_width<UpTypeQuark, DownTypeQuark, WpBoson>(
   const context_base& context,
   typename field_indices<UpTypeQuark>::type const& indexIn,
   typename field_indices<DownTypeQuark>::type const& indexOut1,
   typename field_indices<WpBoson>::type const& indexOut2)
{
   const double mu = context.physical_mass<UpTypeQuark>(indexIn);
   const double md = context.physical_mass<DownTypeQuark>(indexOut1);
   const double mWp = context.physical_mass<WpBoson>(indexOut2);


   const double ps = 1./(8.*Pi) * std::sqrt(KallenLambda(1., Sqr(md/mu), Sqr(mWp/mu)));

   // matrix element squared
   const auto mat_elem_sq = amplitude_squared<UpTypeQuark, DownTypeQuark, WpBoson>(
      context, indexIn, indexOut1, indexOut2);

   // flux * phase space factor * |matrix element|^2
   double result = 0.5/mu * ps * mat_elem_sq;

   // higher order corrections
   if (flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) {
      // t -> b W^+
      if(indexIn.at(0) == 2 && indexOut1.at(0) == 2) {
         const double y = Sqr(mWp/mu);
         static constexpr double CF = 4./3.;
         result *= 1. + get_alphas(context)/(2.*Pi)*CF*Delta_tW(y);
      }
   }

   return result;
}

