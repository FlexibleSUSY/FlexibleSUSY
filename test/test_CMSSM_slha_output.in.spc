Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   0                    # calculate SM pole masses
    4   1                    # pole mass loop order
    5   1                    # EWSB loop order
    6   2                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O(alpha_t^2 + alpha_t alpha_b + alpha_b^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   31   -1                   # loop library (0 = softsusy)
   32   0                    # calculate decays (0 = no, 1 = yes)
Block SOFTSUSY               # SOFTSUSY specific inputs
    1   1.000000000e-04      # tolerance
    2   2                    # up-quark mixing (=1) or down (=2)
    3   0                    # printout
    5   1                    # 2-loop running
    7   1                    # Higgs and EWSB loop order
Block MODSEL                 # Select model
    6    0                   # flavour violation
    1    1                   # mSUGRA
#   12    9.118760000e+01     # output scale (MZ)
Block SMINPUTS               # Standard Model inputs
    1   1.279340000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166370000e-05      # G_Fermi
    3   1.176000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.200000000e+00      # mb(mb) SM MSbar
    6   1.733000000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV)
   22   2.400000000e-03      # mu(2 GeV)
   23   1.040000000e-01      # ms(2 GeV)
   24   1.270000000e+00      # mc(2 GeV)
Block MINPAR                 # Input parameters
    1   125.                 # m0
    2   500.                 # m12
    3   10.                  # TanBeta
    4   1                    # sign(Mu)
    5   0                    # Azero
