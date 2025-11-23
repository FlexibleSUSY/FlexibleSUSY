template <>
double CLASSNAME::get_partial_width<Higgs, Higgs, Higgs>(
      const context_base& context,
      const typename field_indices<Higgs>::type& in_idx,
      const typename field_indices<Higgs>::type& out1_idx,
      const typename field_indices<Higgs>::type& out2_idx)
{
   const auto amp = calculate_amplitude_tree<Higgs, Higgs, Higgs>(context, in_idx, out1_idx, out2_idx);
   const double mHIn = context.physical_mass<Higgs>(in_idx);
   const double mHOut1 = context.physical_mass<Higgs>(out1_idx);
   const double mHOut2 = context.physical_mass<Higgs>(out2_idx);

   const auto xOS1 = Sqr(mHOut1/mHIn);
   const auto xOS2 = Sqr(mHOut2/mHIn);
   const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., xOS1, xOS2));
   const double flux = 0.5/mHIn;
   const double ps_symmetry = boost::range::equal(out1_idx, out2_idx) ? 0.5 : 1;
   auto res = flux * ps * ps_symmetry * amp.square();

   const auto amp_1l = calculate_amplitude_1l<Higgs, Higgs, Higgs>(context, in_idx, out1_idx, out2_idx);

   return res;
}
