Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   0                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   0                    # EFT matching scale (0 = SUSY scale)
   20   2                    # EFT loop order for upwards matching
   21   1                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   1                    # calculate BSM pole masses
   24   123111321            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   1                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   1                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   1                    # Higgs 3-loop corrections O(alpha_t^3)
   30   1                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
   31   -1                   # loop library (0 = softsusy)
   32   0                    # calculate decays (0 = no, 1 = yes)
Block FlexibleSUSYInput
    0   0.00729735           # alpha_em(0)
    1   125.09               # Mh pole
Block SMINPUTS               # Standard Model inputs
    1   1.279400100000000e+02   # alpha^{-1}(mZ)^MSbar
    2   1.166378700000000e-05   # G_Fermi
    3   1.182558783500000e-01   # alpha_s(mZ)^MSbar
    4   9.118760000000000e+02   # mZ(pole)
    5   4.180000000000000e+00   # mb(mb)^MSbar
    6   1.734791074910000e+02   # mtop(pole)
    7   1.776820000000000e+00   # mtau(pole)
    8   0.000000000000000e+00   # mnu3(pole)
    11  5.109989280000000e-04   # melectron(pole)
    12  0.000000000000000e+00   # mnu1(pole)
    13  1.056583720000000e-01   # mmuon(pole)
    14  0.000000000000000e+00   # mnu2(pole)
    21  4.800000000000000e-03   # md(2 GeV)^MSbar
    22  2.300000000000000e-03   # mu(2 GeV)^MSbar
    23  9.500000000000000e-02   # ms(2 GeV)^MSbar
    24  1.275000000000000e+00   # mc(mc)^MSbar
Block MINPAR                 # Input parameters
    1   9130.61751352        # m0
    2   3034.06628664        # m12
    3   52.9427817096        # TanBeta
    4   1                    # SignMu
    5   -45.7405445434       # Azero
