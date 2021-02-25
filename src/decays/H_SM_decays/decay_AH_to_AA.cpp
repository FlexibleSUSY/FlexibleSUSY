template <>
double CLASSNAME::get_partial_width<AH, A, A>(
      const context_base& context,
      const typename field_indices<AH>::type& in_idx,
      const typename field_indices<A>::type& out1_idx,
      const typename field_indices<A>::type& out2_idx) const
{
   const auto amp = calculate_amplitude<AH, A, A>(context, in_idx, out1_idx, out2_idx);
   const double mH = context.physical_mass<AH>(in_idx);
   constexpr double ps {1./(8.*Pi)};
   constexpr double ps_symmetry {1./2.};
   const double flux = 0.5/mH;
   auto res = flux * ps * ps_symmetry * amp.square();

   // use alpha_em in the Thomson limit
   const double alpha_em_0 = physical_input.get(Physical_input::alpha_em_0);
   const double alpha_em = get_alpha(context);
   res *= Sqr(alpha_em_0/alpha_em);

   return res;
}
