import dadi, numpy, scipy
from optparse import OptionParser

def TwoPopAdmix(params,ns,pts):
    driftC,driftA,Tadmix,rate,innerdriftY,innerdriftZ = params    
    Nhum = 1
    Narch = ( (driftC+innerdriftY) * 2 * Nhum) / (2 * driftA)
    Tpre = innerdriftY - Tadmix
    earliest = max([innerdriftY,innerdriftZ,driftC,driftA])*10
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi,xx,(earliest-driftC),nu=Nhum)
    phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
    phi = dadi.Integration.two_pops(phi,xx,driftC,Narch,Nhum,m12=0,m21=0)
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx,phi)
    phi = dadi.Integration.three_pops(phi,xx,Tpre,Narch,Nhum,Nhum)
    phi = dadi.PhiManip.phi_3D_admix_1_and_3_into_2(phi,rate,0,xx,xx,xx)
    phi = dadi.Integration.three_pops(phi,xx,Tadmix,Narch,Nhum,Nhum)
    fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx,xx))
    fs = numpy.array(fs)
#    fs = fs[1:ns[0]]
    return(fs)


parser = OptionParser("$prog [options]")
parser.add_option("-c", "--driftC", dest="driftC", help="Drift C", default=None, type="float")
parser.add_option("-a", "--driftA", dest="driftA", help="Drift A", default=None, type="float")
parser.add_option("-x", "--admixrate", dest="admixrate", help="admixrate", default=None, type="float")
parser.add_option("-t", "--admixtime", dest="admixtime", help="admixtime", default=None, type="float")
parser.add_option("-m", "--nC", dest="nC", help="Number of samples from C", default=None, type="int")
parser.add_option("-b", "--nB", dest="nB", help="Number of samples from B", default=None, type="int")
parser.add_option("-n", "--nA", dest="nA", help="Number of samples from A", default=None, type="int")
parser.add_option("-y", "--innerdriftY",dest="innerdriftY",help="Inner drift Y", default=None, type="float")
parser.add_option("-z", "--innerdriftZ", dest="innerdriftZ", help="Inner drift Z", default=None, type="float")
(options,args) = parser.parse_args()


driftC = options.driftC
driftA = options.driftA
admixrate = options.admixrate
admixtime = options.admixtime
innerdriftY = options.innerdriftY
innerdriftZ = options.innerdriftZ
nC = options.nC
nB = options.nB
nA = options.nA

daditable = TwoPopAdmix([driftC,driftA,admixtime,admixrate,innerdriftY,innerdriftZ],(nA,nC,nB),20)

#print daditable

daditable.tolist()

for i in daditable:
    for j in i:
        for h in j:
            print h
