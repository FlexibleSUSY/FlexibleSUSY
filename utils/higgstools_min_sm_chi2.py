import scipy.optimize

import Higgs.predictions as HP
import Higgs.bounds as HB
import Higgs.signals as HS

bounds = HB.Bounds("/run/media/scratch/Pobrane/hbdataset-master") # load HB dataset
signals = HS.Signals("/run/media/scratch/Pobrane/hsdataset-main") # load HS dataset

pred = HP.Predictions() # create the model predictions
cpls = HP.smLikeEffCouplings

h = pred.addParticle(HP.NeutralScalar("h", "even"))

def f(m):
    h.setMass(m)
    h.setMassUnc(0.0)
    HP.effectiveCouplingInput(h, cpls, HP.ReferenceModel.SMHiggsInterp, calcggH=False, calcHgamgam=False)
    return signals(pred)

res = scipy.optimize.shgo(f, [(120, 130)], n=200000, options={"f_tol": 1e-7, "workers": -1})

print(f'Min. SM chi2 is {res.fun} and happens for mh = {res.x[0]} GeV')
