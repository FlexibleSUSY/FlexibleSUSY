Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   1                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O(alpha_t^2 + alpha_t alpha_b + alpha_b^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
Block SMINPUTS               # Standard Model inputs
    1   1.279440000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166380000e-05      # G_Fermi
    3   1.184000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.425               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block EXTPAR                 # Input parameters
    0   2020                 # MSUSY
    1   2000                 # M1(MSUSY)
    2   2000                 # M2(MSUSY)
    3   2000                 # M3(MSUSY)
    4   2000                 # Mu(MSUSY)
    5   2000                 # mA(MSUSY)
    6   173.34               # MEWSB
    7   200                  # At(MSUSY)
   25   10                   # TanBeta(MSUSY)
  100   2                    # LambdaLoopOrder
  101   1                    # 2-loop at*as
  102   1                    # 2-loop ab*as
  103   1                    # 2-loop at*ab
  104   1                    # 2-loop atau*atau
  105   1                    # 2-loop at*at
  106   1                    # 1-loop ab
  107   1                    # 1-loop atau
Block MSQ2IN
  1  1     4e6               # mq2(1,1)
  2  2     4e6               # mq2(2,2)
  3  3     4e6               # mq2(3,3)
Block MSE2IN
  1  1     4e6               # me2(1,1)
  2  2     4e6               # me2(2,2)
  3  3     4e6               # me2(3,3)
Block MSL2IN
  1  1     4e6               # ml2(1,1)
  2  2     4e6               # ml2(2,2)
  3  3     4e6               # ml2(3,3)
Block MSU2IN
  1  1     4e6               # mu2(1,1)
  2  2     4e6               # mu2(2,2)
  3  3     4e6               # mu2(3,3)
Block MSD2IN
  1  1     4e6               # md2(1,1)
  2  2     4e6               # md2(2,2)
  3  3     4e6               # md2(3,3)
