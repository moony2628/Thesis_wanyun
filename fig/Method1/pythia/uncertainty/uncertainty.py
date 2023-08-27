from ROOT import *
import numpy as np

    
var = "ntrk"  #change the var name according to the inputvar you want to read
mc = "pythia_SF"   #by setting it as "SF" or "MC", it will automatically making scale factor plots or MC closure plots
inputvar = var  #by setting it as bdt (or ntrk,width,c1..), it will read the corresponding histogram, but remember to change the TLine range according to X-axis of different variable, one can check it by browsing the histograms in root file.

if var == "ntrk":
    doreweight = "Quark"  
if var == "bdt":
    doreweight = 0

global rebins
rebins = 0
def rebin(input_hist):
    global rebins
    hist = input_hist.Clone()
    a = np.array([])
    for j in range(1,hist.GetNbinsX()+2):
        if var == "ntrk":
            if j == 1 or (j >= 6 and j<=30 and (j+5)%4 == 0) or j==31 or j ==41 or j ==61: 
                a=np.append(a,hist.GetBinLowEdge(j))
        if var == "bdt":
            if j == 1 or (j >= 13 and j<=49 and (j-1)%4 == 0) or j==61: 
                a=np.append(a,hist.GetBinLowEdge(j))
    rebins = rebins+1
    hist = hist.Rebin(len(a)-1,str(rebins),a)
    return(hist)

def myText(x,y,text, color = 1):
    l = TLatex()
    l.SetTextSize(0.025)
    l.SetNDC()
    l.SetTextColor(color)
    l.DrawLatex(x,y,text)
    pass

def matrixmethod(hq,lq,hg,lg):
    Quark_f = hq.Clone("")
    Gluon_f = hg.Clone("")

    hquark = hq.Clone("")
    hgluon = hg.Clone("")
    lquark = lq.Clone("")
    lgluon = lg.Clone("")
    lower_mc = lquark.Clone()
    higher_mc = hquark.Clone()
    lower_mc.Add(lgluon)
    higher_mc.Add(hgluon)

    fq1 = 0
    cq1 = 0
    fg1 = 0
    cg1 = 0

    for n in range(1,hq.GetNbinsX()+1):
        fq1 += hquark.GetBinContent(n)
        cq1 += lquark.GetBinContent(n)
        fg1 += hgluon.GetBinContent(n)
        cg1 += lgluon.GetBinContent(n)

    fq_f = fq1/(fq1+fg1)
    cq_f = cq1/(cq1+cg1)
    fg_f = 1.-fq_f
    cg_f = 1.-cq_f

    #print(fq_f,cq_f,fg_f,cg_f)

    if(hquark.Integral() != 0):
        hquark.Scale(1./hquark.Integral())
    if(lquark.Integral() != 0):
        lquark.Scale(1./lquark.Integral())
    if(hgluon.Integral() != 0):
        hgluon.Scale(1./hgluon.Integral())
    if(lgluon.Integral() != 0):
        lgluon.Scale(1./lgluon.Integral())
    if(higher_mc.Integral() != 0):
        higher_mc.Scale(1./higher_mc.Integral())
    if(lower_mc.Integral() != 0):
        lower_mc.Scale(1./lower_mc.Integral())
        
    if (doreweight=="Quark"):
            for i in range(1,hquark.GetNbinsX()+1):
                    if (lquark.GetBinContent(i) > 0 and lgluon.GetBinContent(i) > 0):
                            #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                            factor_gluon = hgluon.GetBinContent(i)/lgluon.GetBinContent(i)
                            factor_quark = hquark.GetBinContent(i)/lquark.GetBinContent(i)
                            lquark.SetBinContent(i,lquark.GetBinContent(i)*factor_quark)
                            lgluon.SetBinContent(i,lgluon.GetBinContent(i)*factor_quark)
                            lower_mc.SetBinContent(i,lower_mc.GetBinContent(i)*factor_quark)


        
    if (doreweight=="Gluon"):
            for i in range(1,higher_quark.GetNbinsX()+1):
                    if (lower_quark.GetBinContent(i) > 0 and lower_gluon.GetBinContent(i) > 0):
                            #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                            factor_gluon = higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                            factor_quark = higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i)
                            lquark.SetBinContent(i,lower_quark.GetBinContent(i)*factor_gluon)
                            lgluon.SetBinContent(i,lower_gluon.GetBinContent(i)*factor_gluon)
                            lower_mc.SetBinContent(i,lower_mc.GetBinContent(i)*factor_gluon)


    
    if(hquark.Integral() != 0):
        hquark.Scale(1./hquark.Integral())
    if(lquark.Integral() != 0):
        lquark.Scale(1./lquark.Integral())
    if(hgluon.Integral() != 0):
        hgluon.Scale(1./hgluon.Integral())
    if(lgluon.Integral() != 0):
        lgluon.Scale(1./lgluon.Integral())
    if(lower_mc.Integral() != 0):
        lower_mc.Scale(1./lower_mc.Integral())
    if(higher_mc.Integral() != 0):
        higher_mc.Scale(1./higher_mc.Integral())        
    higher_f = higher_mc.Clone()
    lower_f = lower_mc.Clone()



    for n in range(1,hq.GetNbinsX()+1):
        F_f = higher_f.GetBinContent(n)
        C_f = lower_f.GetBinContent(n)
        if((cg_f*fq_f-fg_f*cq_f) != 0 ):
            Q_f = -(C_f*fg_f-F_f*cg_f)/(cg_f*fq_f-fg_f*cq_f)
            G_f = (C_f*fq_f-F_f*cq_f)/(cg_f*fq_f-fg_f*cq_f)
            #print("Q",Q_f,hquark.GetBinContent(n))
            #print("G",G_f,hgluon.GetBinContent(n))
            Quark_f.SetBinContent(n,Q_f)
            Gluon_f.SetBinContent(n,G_f)

    return Quark_f,Gluon_f


def matrixmethod_data(hq,lq,hg,lg,hd,ld):
    Quark_f = hq.Clone("")
    Gluon_f = hg.Clone("")

    hquark = hq.Clone("")
    hgluon = hg.Clone("")
    lquark = lq.Clone("")
    lgluon = lg.Clone("")
    hdata = hd.Clone("")
    ldata = ld.Clone("")

    fq1 = 0
    cq1 = 0
    fg1 = 0
    cg1 = 0

    for n in range(1,hq.GetNbinsX()+1):
        fq1 += hquark.GetBinContent(n)
        cq1 += lquark.GetBinContent(n)
        fg1 += hgluon.GetBinContent(n)
        cg1 += lgluon.GetBinContent(n)

    fq_f = fq1/(fq1+fg1)
    cq_f = cq1/(cq1+cg1)
    fg_f = 1.-fq_f
    cg_f = 1.-cq_f

    #print(fq_f,cq_f,fg_f,cg_f)

    if(hquark.Integral() != 0):
        hquark.Scale(1./hquark.Integral())
    if(lquark.Integral() != 0):
        lquark.Scale(1./lquark.Integral())
    if(hgluon.Integral() != 0):
        hgluon.Scale(1./hgluon.Integral())
    if(lgluon.Integral() != 0):
        lgluon.Scale(1./lgluon.Integral())
    if(hdata.Integral() != 0):
        hdata.Scale(1./hdata.Integral())
    if(ldata.Integral() != 0):
        ldata.Scale(1./ldata.Integral())
    if (doreweight=="Quark"):
            for i in range(1,hquark.GetNbinsX()+1):
                    if (lquark.GetBinContent(i) > 0 and lgluon.GetBinContent(i) > 0):
                            #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                            factor_gluon = hgluon.GetBinContent(i)/lgluon.GetBinContent(i)
                            factor_quark = hquark.GetBinContent(i)/lquark.GetBinContent(i)
                            lquark.SetBinContent(i,lquark.GetBinContent(i)*factor_quark)
                            lgluon.SetBinContent(i,lgluon.GetBinContent(i)*factor_quark)
                            ldata.SetBinContent(i,ldata.GetBinContent(i)*factor_quark)

                            pass
                    pass
            pass
        
    if (doreweight=="Gluon"):
            for i in range(1,higher_quark.GetNbinsX()+1):
                    if (lower_quark.GetBinContent(i) > 0 and lower_gluon.GetBinContent(i) > 0):
                            #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                            factor_gluon = higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                            factor_quark = higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i)
                            lquark.SetBinContent(i,lower_quark.GetBinContent(i)*factor_gluon)
                            lgluon.SetBinContent(i,lower_gluon.GetBinContent(i)*factor_gluon)
                            ldata.SetBinContent(i,ldata.GetBinContent(i)*factor_gluon)
                            pass
                    pass
            pass
    
    if(hquark.Integral() != 0):
        hquark.Scale(1./hquark.Integral())
    if(lquark.Integral() != 0):
        lquark.Scale(1./lquark.Integral())
    if(hgluon.Integral() != 0):
        hgluon.Scale(1./hgluon.Integral())
    if(lgluon.Integral() != 0):
        lgluon.Scale(1./lgluon.Integral())
    if(hdata.Integral() != 0):
        hdata.Scale(1./hdata.Integral())
    if(ldata.Integral() != 0):
        ldata.Scale(1./ldata.Integral())
        
    for n in range(1,hq.GetNbinsX()+1):
        F_f = hdata.GetBinContent(n)
        C_f = ldata.GetBinContent(n)
        if((cg_f*fq_f-fg_f*cq_f) != 0 ):
            Q_f = -(C_f*fg_f-F_f*cg_f)/(cg_f*fq_f-fg_f*cq_f)
            G_f = (C_f*fq_f-F_f*cq_f)/(cg_f*fq_f-fg_f*cq_f)
            #print("Q",Q_f,hquark.GetBinContent(n))
            #print("G",G_f,hgluon.GetBinContent(n))
            Quark_f.SetBinContent(n,Q_f)
            Gluon_f.SetBinContent(n,G_f)

    return Quark_f,Gluon_f



def percentdifference(data_quark_1,data_gluon_1,mc_quark_1,mc_gluon_1,data_quark_2,data_gluon_2,mc_quark_2,mc_gluon_2): 
    sigma_q = np.zeros(data_quark_1.GetNbinsX())
    sigma_g = np.zeros(data_quark_1.GetNbinsX())

    for j in range(1,data_quark_1.GetNbinsX()+1):
                error = 0
                value = 0
                sq = 0
                if mc_quark_1.GetBinContent(j) != 0  and mc_quark_2.GetBinContent(j) != 0 :
                    error = data_quark_1.GetBinContent(j)/mc_quark_1.GetBinContent(j) - data_quark_2.GetBinContent(j)/mc_quark_2.GetBinContent(j)
                if mc_quark_1.GetBinContent(j) != 0:
                    value = data_quark_1.GetBinContent(j)/mc_quark_1.GetBinContent(j)
                if value != 0:
                    sq = abs(error/value)
                
                error = 0
                value = 0
                sg = 0
                if mc_gluon_1.GetBinContent(j) != 0  and mc_gluon_2.GetBinContent(j) != 0 :
                    error = data_gluon_1.GetBinContent(j)/mc_gluon_1.GetBinContent(j) - data_gluon_2.GetBinContent(j)/mc_gluon_2.GetBinContent(j)
                if mc_gluon_1.GetBinContent(j) != 0:
                    value = data_gluon_1.GetBinContent(j)/mc_gluon_1.GetBinContent(j)
                if value != 0:
                    sg = abs(error/value)

                sigma_q[j-1] = sq
                sigma_g[j-1] = sg

    return sigma_q,sigma_g

def bootstrap_result(data):
    w = 58.45/39.91 # prescale factor 
    n1 = data.Clone()
    n2 = data.Clone()
    result = data.Clone()
    for j in range(1,higher_quark.GetNbinsX()+1):
            n2.SetBinContent(j,int((data.GetBinContent(j)-data.GetBinError(j)**2)/(w-w**2)))#(events num for prescaled data)
            n1.SetBinContent(j,int((data.GetBinContent(j)-w*n2.GetBinContent(j))))
            result.SetBinContent(j,w*np.random.poisson(n2.GetBinContent(j))+np.random.poisson(n1.GetBinContent(j)))
    return(result)    

c1 = TCanvas("","",500,500)

bin = [0, 50, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 2000]


fpythia = TFile("../newroot/dijet_pythia_nominal.root")
ntrackall = TFile("../newroot/dijet_sherpa.root")
ntrackall3 = TFile("../newroot/dijet_data_prescale.root")
#ntrackall5 = TFile("../newroot/dijet-sherpa-pdf.root")
fherdipo = TFile("../newroot/dijet_herdipo.root")
fherang = TFile("../newroot/dijet_herang.root")
sherpa_lund = TFile("../newroot/dijet_sherpa_lund.root")
fpowpyt = TFile("../newroot/dijet_powheg.root")
#sv = TFile("../newroot/dijet-sv.root")

for i in range(7,13):   #for only dijet event, start frTom jet pT>500 GeV
#for i in range(13):    #for gamma+jet combined with dijet event, start from jet pT>0 GeV
#        if(bin[i] != 800):
            c = TCanvas("c","c",500,500)
            min = bin[i]
            max = bin[i+1]

            higher_quark = fpythia.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
            higher_quark2 = fpythia.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar)
            higher_gluon = fpythia.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
            higher_gluon2 = fpythia.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar)
            lower_quark = fpythia.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar)
            lower_quark2 = fpythia.Get(str(min)+"_SubJet_Central_Quark_"+inputvar)
            lower_gluon = fpythia.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
            lower_gluon2 = fpythia.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar)

            higher_data = ntrackall3.Get(str(min)+"_LeadingJet_Forward_Data_"+inputvar)
            higher_data2 = ntrackall3.Get(str(min)+"_SubJet_Forward_Data_"+inputvar)
            lower_data = ntrackall3.Get(str(min)+"_LeadingJet_Central_Data_"+inputvar)
            lower_data2 = ntrackall3.Get(str(min)+"_SubJet_Central_Data_"+inputvar)

            #add leading and subleading jet from only dijet event together,
            #note that for gammajet+dijet event, we need to add leading jet from gammajet and leading jet from dijet sample together
            higher_data.Add(higher_data2)
            lower_data.Add(lower_data2)
            higher_quark.Add(higher_quark2)
            higher_gluon.Add(higher_gluon2)
            lower_quark.Add(lower_quark2)
            lower_gluon.Add(lower_gluon2)
            
            higher_data = rebin(higher_data) 
            lower_data = rebin(lower_data) 
            higher_quark = rebin(higher_quark) 
            lower_quark = rebin(lower_quark) 
            higher_gluon = rebin(higher_gluon) 
            lower_gluon = rebin(lower_gluon) 
            
            
            nominal_extract_Q,nominal_extract_G = matrixmethod(higher_quark,lower_quark,higher_gluon,lower_gluon)
            extracted_data_nominal_Q,extracted_data_nominal_G = matrixmethod_data(higher_quark,lower_quark,higher_gluon,lower_gluon,higher_data,lower_data)


            higher_data_strap = higher_data.Clone("")     #Set aside for statistical uncertainty
            lower_data_strap = lower_data.Clone("")

            pdf_qvals = []
            pdf_gvals = []

            for j in range(1,higher_quark.GetNbinsX()+1):
                    pdf_qvals += [np.zeros(50)]
                    pdf_gvals += [np.zeros(50)]

            #uncertainty calculations
            #uncertainty lists, number-of-bins lists of 4 uncertainties.
            sigma_tot_q = []
            sigma_tot_g = []
            ebar_q = []
            ebar_g = []

            for j in range(0,higher_quark.GetNbinsX()):
                    sigma_tot_q += [np.zeros(8)]
                    sigma_tot_g += [np.zeros(8)]
                    ebar_q += [np.zeros(8)]
                    ebar_g += [np.zeros(8)]

            # do bootstrap
            #1. create lists to store bootstrapped values list of arrays of nstraps values
            
            
            
            nstraps = 5000
            SFQvals = []
            SFGvals = []
            
            for j in range(1,higher_quark.GetNbinsX()+1):
                    SFQvals += [np.zeros(nstraps)]
                    SFGvals += [np.zeros(nstraps)]

            ToT_Fq2 = 0.
            ToT_Fg2 = 0.

            ToT_Cq2 = 0.
            ToT_Cg2 = 0.
            
            higher_quark_bootstrap = higher_quark.Clone()
            higher_gluon_bootstrap = higher_gluon.Clone()
            lower_quark_bootstrap = lower_quark.Clone()
            lower_gluon_bootstrap = lower_gluon.Clone()
            higher_data_bootstrap = higher_data.Clone()
            lower_data_bootstrap = lower_data.Clone()
            
            for j in range(1,higher_quark.GetNbinsX()+1):
                    ToT_Fq2+=higher_quark_bootstrap.GetBinContent(j)
                    ToT_Cq2+=lower_quark_bootstrap.GetBinContent(j)
                    ToT_Fg2+=higher_gluon_bootstrap.GetBinContent(j)
                    ToT_Cg2+=lower_gluon_bootstrap.GetBinContent(j)

            fg=ToT_Fg2/(ToT_Fg2+ToT_Fq2)
            cg=ToT_Cg2/(ToT_Cq2+ToT_Cg2)
            fq=1.-fg
            cq=1.-cg
            

            #do bootsrapping
            for k in range(nstraps):
                    lower_data_bootstrap = lower_data.Clone()
                    higher_data_bootstrap = higher_data.Clone()
                    #print(higher_data_bootstrap.GetBinContent(3))
                    lower_data_bootstrap = bootstrap_result(lower_data_bootstrap)
                    higher_data_bootstrap = bootstrap_result(higher_data_bootstrap)

                                
                    if(lower_data_bootstrap.Integral() != 0):
                            lower_data_bootstrap.Scale(1./lower_data_bootstrap.Integral())
                    if(higher_data_bootstrap.Integral() != 0):
                            higher_data_bootstrap.Scale(1./higher_data_bootstrap.Integral())
                                       
                    if (doreweight=="Quark"):
                            for j in range(1,higher_quark_bootstrap.GetNbinsX()+1):
                                    if (lower_quark_bootstrap.GetBinContent(j) > 0 and lower_gluon_bootstrap.GetBinContent(j) > 0):
                                            #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                            factor_gluon = higher_gluon_bootstrap.GetBinContent(j)/lower_gluon_bootstrap.GetBinContent(j)
                                            factor_quark = higher_quark_bootstrap.GetBinContent(j)/lower_quark_bootstrap.GetBinContent(j)
                                            lower_data_bootstrap.SetBinContent(j,lower_data_bootstrap.GetBinContent(j)*factor_quark)
                                            pass
                                    pass
                            pass
                        
                    if (doreweight=="Gluon"):
                            for j in range(1,higher_quark_bootstrap.GetNbinsX()+1):
                                    if (lower_quark_bootstrap.GetBinContent(j) > 0 and lower_gluon_bootstrap.GetBinContent(j) > 0):
                                            #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                            factor_gluon = higher_gluon_bootstrap.GetBinContent(j)/lower_gluon_bootstrap.GetBinContent(j)
                                            factor_quark = higher_quark_bootstrap.GetBinContent(j)/lower_quark_bootstrap.GetBinContent(j)
                                            lower_data_bootstrap.SetBinContent(j,lower_data_bootstrap.GetBinContent(j)*factor_gluon)
                                            pass
                                    pass
                            pass
                        

                    if(lower_data_bootstrap.Integral() != 0):
                            lower_data_bootstrap.Scale(1./lower_data_bootstrap.Integral())
                    if(higher_data_bootstrap.Integral() != 0):
                            higher_data_bootstrap.Scale(1./higher_data_bootstrap.Integral())
                            
                    forward_data_strap = higher_data_bootstrap.Clone("f"+str(k))
                    central_data_strap = lower_data_bootstrap.Clone("c"+str(k))
                 
                    if(forward_data_strap.Integral() != 0):
                            forward_data_strap.Scale(1./forward_data_strap.Integral())
                    if(central_data_strap.Integral() != 0):
                            central_data_strap.Scale(1./central_data_strap.Integral())

                    for j in range(0,higher_quark.GetNbinsX()):
                            F = forward_data_strap.GetBinContent(j+1)
                            C = central_data_strap.GetBinContent(j+1)
                            Q = -(C*fg-F*cg)/(cg*fq-fg*cq)
                            G = (C*fq-F*cq)/(cg*fq-fg*cq)
                            if nominal_extract_Q.GetBinContent(j) != 0:
                                SFQvals[j][k] = Q/nominal_extract_Q.GetBinContent(j)
                            if nominal_extract_Q.GetBinContent(j) != 0:
                                SFGvals[j][k] = G/nominal_extract_G.GetBinContent(j)
                                
            quark_data = higher_quark.Clone()
            gluon_data = higher_gluon.Clone()

            #compute the uncertainty and plots
            quark_strap = higher_data.Clone()
            gluon_strap = higher_data.Clone()
            
            for j in range(0,quark_data.GetNbinsX()):
                    SFQvals[j].sort()
                    SFGvals[j].sort()
                    Q = np.median(SFQvals[j])
                    G = np.median(SFGvals[j])

                    sigmaQ = .5*(SFQvals[j][int(.84*len(SFQvals[j]))] - SFQvals[j][int(.16*len(SFQvals[j]))])
                    sigmaG = .5*(SFGvals[j][int(.84*len(SFGvals[j]))] - SFGvals[j][int(.16*len(SFGvals[j]))])


            #ebar_q[j][0] = sigmaQ/np.sum(sigmaQ[j])
            #ebar_g[j][0] = sigmaG/np.sum(sigmaG[j])

                    #print("statistical: q = "+str(sigmaQ)+" | g = "+str(sigmaG))

                    if(Q != 0):
                            sigmaQ = np.abs(sigmaQ/Q)
                    if(G != 0):
                            sigmaG = np.abs(sigmaG/G)

                    sigma_tot_q[j][0] = sigmaQ
                    sigma_tot_g[j][0] = sigmaG
                    #print("tsetflag:",j,sigma_tot_q[j][0])
                    quark_strap.SetBinContent(j+1,sigmaQ)
                    gluon_strap.SetBinContent(j+1,sigmaG)

            quark_negative = quark_strap.Clone()
            gluon_negative = gluon_strap.Clone()

            quark_negative = quark_negative * -1
            gluon_negative = gluon_negative * -1

            #mc uncertainty
            #uncertainty calculation percent difference
            higher_quark_copy = higher_quark.Clone()
            higher_gluon_copy = higher_gluon.Clone()
            
            
            higher_quark_copy.Scale(1./higher_quark_copy.Integral())
            higher_gluon_copy.Scale(1./higher_gluon_copy.Integral())
            quark,gluon = matrixmethod(higher_quark,lower_quark,higher_gluon,lower_gluon)
            qmc,gmc = percentdifference(extracted_data_nominal_Q,extracted_data_nominal_G,nominal_extract_Q,nominal_extract_G,extracted_data_nominal_Q,extracted_data_nominal_G,higher_quark_copy,higher_gluon_copy)

            quark_use = quark.Clone()
            gluon_use = gluon.Clone()

            for j in range(1,quark.GetNbinsX()+1):
                    quark_use.SetBinContent(j,qmc[j-1])
                    gluon_use.SetBinContent(j,gmc[j-1])
                    sigma_tot_q[j-1][1] = qmc[j-1]
                    sigma_tot_g[j-1][1] = gmc[j-1]
            #if i == 7:
            #    for j in range(1,quark.GetNbinsX()+1):
            #        print(j,gluon_use.GetBinContent(j),higher_gluon_copy.GetBinContent(j),gluon.GetBinContent(j))
            #    print("all",higher_gluon_copy.Integral(),gluon.Integral())
            quarkMC_negative = quark_use.Clone()
            gluonMC_negative = gluon_use.Clone()

            quarkMC_negative.Scale(-1)
            gluonMC_negative.Scale(-1)
            ## reweight factor uncertainty
            doreweight_origin = doreweight
            
            higher_quark_RF_Q = higher_quark.Clone()
            higher_gluon_RF_Q = higher_gluon.Clone()
            lower_quark_RF_Q = lower_quark.Clone()
            lower_gluon_RF_Q = lower_gluon.Clone()
            higher_data_RF_Q = higher_data.Clone()
            lower_data_RF_Q = lower_data.Clone()
            
            higher_quark_RF_G = higher_quark.Clone()
            higher_gluon_RF_G = higher_gluon.Clone()
            lower_quark_RF_G = lower_quark.Clone()
            lower_gluon_RF_G = lower_gluon.Clone()
            higher_data_RF_G = higher_data.Clone()
            lower_data_RF_G = lower_data.Clone()    
            
            
            doreweight = "Quark"
            quark_RF_Q_data,gluon_RF_Q_data = matrixmethod_data(higher_quark_RF_Q,lower_quark_RF_Q,higher_gluon_RF_Q,lower_gluon_RF_Q,higher_data_RF_Q,lower_data_RF_Q)
            quark_RF_Q_mc,gluon_RF_Q_mc = matrixmethod(higher_quark_RF_Q,lower_quark_RF_Q,higher_gluon_RF_Q,lower_gluon_RF_Q)

            doreweight = "Gluon"            
            quark_RF_G_data,gluon_RF_G_data = matrixmethod_data(higher_quark_RF_G,lower_quark_RF_G,higher_gluon_RF_G,lower_gluon_RF_G,higher_data_RF_G,lower_data_RF_G)
            quark_RF_G_mc,gluon_RF_G_mc = matrixmethod(higher_quark_RF_G,lower_quark_RF_G,higher_gluon_RF_G,lower_gluon_RF_G)
            
            doreweight = doreweight_origin
            qRF,gRF= percentdifference(quark_RF_Q_data,gluon_RF_Q_data,quark_RF_Q_mc,gluon_RF_Q_mc,quark_RF_G_data,gluon_RF_G_data,quark_RF_G_mc,gluon_RF_G_mc)
            quark_RF = quark.Clone()
            gluon_RF = gluon.Clone()

            for j in range(1,quark_RF.GetNbinsX()+1):
                    quark_RF.SetBinContent(j,qRF[j-1])
                    gluon_RF.SetBinContent(j,gRF[j-1])
                    #ebar_q[j-1][1] = np.abs(quark.GetBinContent(j) - higher_quark_copy.GetBinContent(j))
                    #ebar_g[j-1][1] = np.abs(gluon.GetBinContent(j) - higher_gluon_copy.GetBinContent(j))
                    #sigma_tot_q[j-1][1] = qmc[j-1]
                    #sigma_tot_g[j-1][1] = gmc[j-1]
            #if i == 7:
            #    for j in range(1,quark.GetNbinsX()+1):
            #        print(j,gluon_use.GetBinContent(j),higher_gluon_copy.GetBinContent(j),gluon.GetBinContent(j))
            #    print("all",higher_gluon_copy.Integral(),gluon.Integral())
            quarkRF_negative = quark_RF.Clone()
            gluonRF_negative = gluon_RF.Clone()

            quarkRF_negative.Scale(-1)
            gluonRF_negative.Scale(-1)            
            
            #pdf uncertainty. stdev of binvals
            #open the histograms for each pdf weight.
            """###
            index = -1
            for k in range(50):
                    index = int(index + 1)
                    
                    
                    higher_quark_pdf = ntrackall5.Get(str(min)+"_LeadingJet_Forward_Quark_"+str(k)+"_"+inputvar)
                    higher_quark1 = ntrackall5.Get(str(min)+"_SubJet_Forward_Quark"+str(k)+"_"+inputvar)
                    lower_quark_pdf = ntrackall5.Get(str(min)+"_LeadingJet_Central_Quark_"+str(k)+"_"+inputvar)
                    lower_quark1 = ntrackall5.Get(str(min)+"_SubJet_Central_Quark"+str(k)+"_"+inputvar)
                    higher_gluon_pdf = ntrackall5.Get(str(min)+"_LeadingJet_Forward_Gluon_"+str(k)+"_"+inputvar)
                    higher_gluon1 = ntrackall5.Get(str(min)+"_SubJet_Forward_Gluon"+str(k)+"_"+inputvar)
                    lower_gluon_pdf = ntrackall5.Get(str(min)+"_LeadingJet_Central_Gluon_"+str(k)+"_"+inputvar)
                    lower_gluon1 = ntrackall5.Get(str(min)+"_SubJet_Central_Gluon"+str(k)+"_"+inputvar)
                    higher_quark_pdf.SetName("test")
                    lower_quark_pdf.SetName("test")
                    higher_gluon_pdf.SetName("test")
                    lower_gluon_pdf.SetName("test")
                    higher_quark1.SetName("test")
                    lower_quark1.SetName("test")
                    higher_gluon1.SetName("test")
                    lower_gluon1.SetName("test")
                    higher_quark_pdf.Add(higher_quark1)
                    higher_gluon_pdf.Add(higher_gluon1)
                    lower_quark_pdf.Add(lower_quark1)
                    lower_gluon_pdf.Add(lower_gluon1)
                    higher_quark_pdf = rebin(higher_quark_pdf) 
                    lower_quark_pdf = rebin(lower_quark_pdf) 
                    higher_gluon_pdf = rebin(higher_gluon_pdf) 
                    lower_gluon_pdf = rebin(lower_gluon_pdf) 
                    
                    quark_pdf_mc,gluon_pdf_mc = matrixmethod(higher_quark_pdf,lower_quark_pdf,higher_gluon_pdf,lower_gluon_pdf)
                    quark_pdf_data,gluon_pdf_data = matrixmethod_data(higher_quark_pdf,lower_quark_pdf,higher_gluon_pdf,lower_gluon_pdf,higher_data,lower_data)
                    for j in range(1,quark_pdf_data.GetNbinsX()+1):
                        if quark_pdf_mc.GetBinContent(j) != 0: 
                            pdf_qvals[j-1][index] = quark_pdf_data.GetBinContent(j)/quark_pdf_mc.GetBinContent(j)
                        if gluon_pdf_mc.GetBinContent(j) != 0: 
                            pdf_gvals[j-1][index] = gluon_pdf_data.GetBinContent(j)/gluon_pdf_mc.GetBinContent(j)

            quark_pdf = quark.Clone("")
            gluon_pdf = quark.Clone("")

            for j in range(0,quark.GetNbinsX()):
                    pdf_qvals[j].sort()
                    pdf_gvals[j].sort()
                    Q = np.median(pdf_qvals[j])
                    G = np.median(pdf_gvals[j])

                    pdf_sigmaQ = .5*(pdf_qvals[j][int(.84*len(pdf_qvals[j]))] - pdf_qvals[j][int(.16*len(pdf_qvals[j]))])
                    pdf_sigmaG = .5*(pdf_gvals[j][int(.84*len(pdf_gvals[j]))] - pdf_gvals[j][int(.16*len(pdf_gvals[j]))])

                    #print("PDF: q = "+str(pdf_sigmaQ)+" | g = "+str(pdf_sigmaG))

                    if(Q != 0):
                            pdf_sigmaQ = np.abs(pdf_sigmaQ/Q)
                    if(G != 0):
                            pdf_sigmaG = np.abs(pdf_sigmaG/G)

                    sigma_tot_q[j][3] = pdf_sigmaQ
                    sigma_tot_g[j][3] = pdf_sigmaG

                    quark_pdf.SetBinContent(j+1,pdf_sigmaQ)
                    gluon_pdf.SetBinContent(j+1,pdf_sigmaG)

            quark_pdf_negative = quark_pdf.Clone("")
            gluon_pdf_negative = gluon_pdf.Clone("")

            quark_pdf_negative = quark_pdf_negative * -1
            gluon_pdf_negative = gluon_pdf_negative * -1
            """###


            #hadronization, difference in sherpa
            hqlund = sherpa_lund.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
            hqlund2 = sherpa_lund.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar)
            hglund = sherpa_lund.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
            hglund2 = sherpa_lund.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar)
            lqlund = sherpa_lund.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar)
            lqlund2 = sherpa_lund.Get(str(min)+"_SubJet_Central_Quark_"+inputvar)
            lglund = sherpa_lund.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
            lglund2 = sherpa_lund.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar)

            hqlund.Add(hqlund2)
            hglund.Add(hglund2)
            lqlund.Add(lqlund2)
            lglund.Add(lglund2)
            
            hqlund = rebin(hqlund) 
            hglund = rebin(hglund) 
            lqlund = rebin(lqlund) 
            lglund = rebin(lglund) 
            
            higher_quark_sherpa = ntrackall.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
            higher_quark2_sherpa = ntrackall.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar)
            higher_gluon_sherpa = ntrackall.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
            higher_gluon2_sherpa = ntrackall.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar)
            lower_quark_sherpa = ntrackall.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar)
            lower_quark2_sherpa = ntrackall.Get(str(min)+"_SubJet_Central_Quark_"+inputvar)
            lower_gluon_sherpa = ntrackall.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
            lower_gluon2_sherpa = ntrackall.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar)

            higher_quark_sherpa.Add(higher_quark2_sherpa)
            higher_gluon_sherpa.Add(higher_gluon2_sherpa)
            lower_quark_sherpa.Add(lower_quark2_sherpa)
            lower_gluon_sherpa.Add(lower_gluon2_sherpa)
            
            higher_quark_sherpa = rebin(higher_quark_sherpa) 
            lower_quark_sherpa = rebin(lower_quark_sherpa) 
            higher_gluon_sherpa = rebin(higher_gluon_sherpa) 
            lower_gluon_sherpa = rebin(lower_gluon_sherpa) 
            
            lower_quark_lund =  lqlund.Clone("")
            lower_gluon_lund =  lglund.Clone("")
            higher_quark_lund =  hqlund.Clone("")
            higher_gluon_lund =  hglund.Clone("")
            
            lund_extract_Q,lund_extract_G = matrixmethod(higher_quark_lund,lower_quark_lund,higher_gluon_lund,lower_gluon_lund)
            extracted_data_lund_Q,extracted_data_lund_G = matrixmethod_data(higher_quark_lund,lower_quark_lund,higher_gluon_lund,lower_gluon_lund,higher_data,lower_data)

            sherpa_extract_Q,sherpa_extract_G = matrixmethod(higher_quark_sherpa,lower_quark_sherpa,higher_gluon_sherpa,lower_gluon_sherpa)
            extracted_data_sherpa_Q,extracted_data_sherpa_G = matrixmethod_data(higher_quark_sherpa,lower_quark_sherpa,higher_gluon_sherpa,lower_gluon_sherpa,higher_data,lower_data)
            qhadunc = hqlund.Clone()
            ghadunc = hqlund.Clone()
            qhad,ghad = percentdifference(extracted_data_sherpa_Q,extracted_data_sherpa_G,sherpa_extract_Q,sherpa_extract_G,extracted_data_lund_Q,extracted_data_lund_G,lund_extract_Q,lund_extract_G)

            for j in range(1,hqlund.GetNbinsX()+1):
                    qhadunc.SetBinContent(j,qhad[j-1])
                    ghadunc.SetBinContent(j,ghad[j-1])
                    sigma_tot_q[j-1][5] = qhad[j-1]
                    sigma_tot_g[j-1][5] = ghad[j-1]

            qhadn = qhadunc.Clone()
            ghadn = ghadunc.Clone()

            qhadn.Scale(-1)
            ghadn.Scale(-1)

            for j in range(0,hqlund.GetNbinsX()):
                    sigma_tot_q[j][5] = qhadunc.GetBinContent(j+1)
                    sigma_tot_g[j][5] = ghadunc.GetBinContent(j+1)
            ### showering ###                     
            higher_quark_herdipo = fherdipo.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
            higher_quark2_herdipo = fherdipo.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar)
            higher_gluon_herdipo = fherdipo.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
            higher_gluon2_herdipo = fherdipo.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar)
            lower_quark_herdipo = fherdipo.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar)
            lower_quark2_herdipo = fherdipo.Get(str(min)+"_SubJet_Central_Quark_"+inputvar)
            lower_gluon_herdipo = fherdipo.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
            lower_gluon2_herdipo = fherdipo.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar)

            higher_quark_herang = fherang.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
            higher_quark2_herang = fherang.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar)
            higher_gluon_herang = fherang.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
            higher_gluon2_herang = fherang.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar)
            lower_quark_herang = fherang.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar)
            lower_quark2_herang = fherang.Get(str(min)+"_SubJet_Central_Quark_"+inputvar)
            lower_gluon_herang = fherang.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
            lower_gluon2_herang = fherang.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar)

            higher_quark_herdipo.Add(higher_quark2_herdipo)
            higher_gluon_herdipo.Add(higher_gluon2_herdipo)
            lower_quark_herdipo.Add(lower_quark2_herdipo)
            lower_gluon_herdipo.Add(lower_gluon2_herdipo)
            
            higher_quark_herdipo = rebin(higher_quark_herdipo) 
            lower_quark_herdipo = rebin(lower_quark_herdipo) 
            higher_gluon_herdipo = rebin(higher_gluon_herdipo) 
            lower_gluon_herdipo = rebin(lower_gluon_herdipo) 
            
            higher_quark_herang.Add(higher_quark2_herang)
            higher_gluon_herang.Add(higher_gluon2_herang)
            lower_quark_herang.Add(lower_quark2_herang)
            lower_gluon_herang.Add(lower_gluon2_herang)                    

            higher_quark_herang = rebin(higher_quark_herang) 
            lower_quark_herang = rebin(lower_quark_herang) 
            higher_gluon_herang = rebin(higher_gluon_herang) 
            lower_gluon_herang = rebin(lower_gluon_herang) 
    
            herdipo_extract_Q,herdipo_extract_G = matrixmethod(higher_quark_herdipo,lower_quark_herdipo,higher_gluon_herdipo,lower_gluon_herdipo)
            extracted_data_herdipo_Q,extracted_data_herdipo_G = matrixmethod_data(higher_quark_herdipo,lower_quark_herdipo,higher_gluon_herdipo,lower_gluon_herdipo,higher_data,lower_data)
            herang_extract_Q,herang_extract_G = matrixmethod(higher_quark_herang,lower_quark_herang,higher_gluon_herang,lower_gluon_herang)
            extracted_data_herang_Q,extracted_data_herang_G = matrixmethod_data(higher_quark_herang,lower_quark_herang,higher_gluon_herang,lower_gluon_herang,higher_data,lower_data)


            qshow,gshow = percentdifference(extracted_data_herdipo_Q,extracted_data_herdipo_G,herdipo_extract_Q,herdipo_extract_G,extracted_data_herang_Q,extracted_data_herang_G,herang_extract_Q,herang_extract_G)

            q_show_unc = herdipo_extract_Q.Clone()
            g_show_unc = herdipo_extract_G.Clone()

            for j in range(1,herdipo_extract_Q.GetNbinsX()+1):
                    #print(j,qshow[j-1])
                    q_show_unc.SetBinContent(j,qshow[j-1])
                    g_show_unc.SetBinContent(j,gshow[j-1])
                    sigma_tot_q[j-1][4] = qshow[j-1]
                    sigma_tot_g[j-1][4] = gshow[j-1]

            q_show_uncn = q_show_unc.Clone()
            g_show_uncn = g_show_unc.Clone()

            q_show_uncn.Scale(-1)
            g_show_uncn.Scale(-1)



            ### Matrix Element ###                     
            higher_quark_powpyt = fpowpyt.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
            higher_quark2_powpyt = fpowpyt.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar)
            higher_gluon_powpyt = fpowpyt.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
            higher_gluon2_powpyt = fpowpyt.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar)
            lower_quark_powpyt = fpowpyt.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar)
            lower_quark2_powpyt = fpowpyt.Get(str(min)+"_SubJet_Central_Quark_"+inputvar)
            lower_gluon_powpyt = fpowpyt.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
            lower_gluon2_powpyt = fpowpyt.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar)


            higher_quark_powpyt.Add(higher_quark2_powpyt)
            higher_gluon_powpyt.Add(higher_gluon2_powpyt)
            lower_quark_powpyt.Add(lower_quark2_powpyt)
            lower_gluon_powpyt.Add(lower_gluon2_powpyt)

            higher_quark_powpyt = rebin(higher_quark_powpyt) 
            lower_quark_powpyt = rebin(lower_quark_powpyt) 
            higher_gluon_powpyt = rebin(higher_gluon_powpyt) 
            lower_gluon_powpyt = rebin(lower_gluon_powpyt) 
            

    
            powpyt_extract_Q,powpyt_extract_G = matrixmethod(higher_quark_powpyt,lower_quark_powpyt,higher_gluon_powpyt,lower_gluon_powpyt)
            extracted_data_powpyt_Q,extracted_data_powpyt_G = matrixmethod_data(higher_quark_powpyt,lower_quark_powpyt,higher_gluon_powpyt,lower_gluon_powpyt,higher_data,lower_data)

            qme,gme = percentdifference(extracted_data_powpyt_Q,extracted_data_powpyt_G,powpyt_extract_Q,powpyt_extract_G,extracted_data_nominal_Q,extracted_data_nominal_G,nominal_extract_Q,nominal_extract_G)

            qmeunc = powpyt_extract_Q.Clone()
            gmeunc = powpyt_extract_G.Clone()

            for j in range(1,powpyt_extract_Q.GetNbinsX()+1):
                    #print(qme[j-1])
                    qmeunc.SetBinContent(j,qme[j-1])
                    gmeunc.SetBinContent(j,gme[j-1])
                    sigma_tot_q[j-1][6] = qme[j-1]
                    sigma_tot_g[j-1][6] = gme[j-1]

            qmen = qmeunc.Clone()
            gmen = gmeunc.Clone()

            qmen.Scale(-1)
            gmen.Scale(-1)  
            """###
            ### scale variation ### 
                
            for k in range(7):
                higher_quark_SV = sv.Get(str(min)+"_LeadingJet_Forward_Quark_"+str(k)+"_" + inputvar)
                higher_quark2_SV = sv.Get(str(min)+"_SubJet_Forward_Quark_"+str(k)+"_" + inputvar)
                higher_gluon_SV = sv.Get(str(min)+"_LeadingJet_Forward_Gluon_"+str(k)+"_" + inputvar)
                higher_gluon2_SV = sv.Get(str(min)+"_SubJet_Forward_Gluon_"+str(k)+"_" + inputvar)
                lower_quark_SV = sv.Get(str(min)+"_LeadingJet_Central_Quark_"+str(k)+"_" + inputvar)
                lower_quark2_SV = sv.Get(str(min)+"_SubJet_Central_Quark_"+str(k)+"_" + inputvar)
                lower_gluon_SV = sv.Get(str(min)+"_LeadingJet_Central_Gluon_"+str(k)+"_" + inputvar)
                lower_gluon2_SV = sv.Get(str(min)+"_SubJet_Central_Gluon_"+str(k)+"_" + inputvar)
                
                higher_quark_SV.Add(higher_quark2_SV)
                higher_gluon_SV.Add(higher_gluon2_SV)
                lower_quark_SV.Add(lower_quark2_SV)
                lower_gluon_SV.Add(lower_gluon2_SV)
                
                higher_quark_SV = rebin(higher_quark_SV) 
                lower_quark_SV = rebin(lower_quark_SV) 
                higher_gluon_SV = rebin(higher_gluon_SV) 
                lower_gluon_SV = rebin(lower_gluon_SV)
                
                SV_extract_Q,SV_extract_G = matrixmethod(higher_quark_SV,lower_quark_SV,higher_gluon_SV,lower_gluon_SV)
                extracted_data_SV_Q,extracted_data_SV_G = matrixmethod_data(higher_quark_SV,lower_quark_SV,higher_gluon_SV,lower_gluon_SV,higher_data,lower_data)
                if k == 0:
                    prevq,prevg = percentdifference(extracted_data_nominal_Q,extracted_data_nominal_G,nominal_extract_Q,nominal_extract_G,extracted_data_SV_Q,extracted_data_SV_G,SV_extract_Q,SV_extract_G)

                if k > 0:
                    diffq,diffg = percentdifference(extracted_data_nominal_Q,extracted_data_nominal_G,nominal_extract_Q,nominal_extract_G,extracted_data_SV_Q,extracted_data_SV_G,SV_extract_Q,SV_extract_G)
                    for l in range(0,SV_extract_Q.GetNbinsX()):
                        if(diffq[l]>prevq[l]):
                            prevq = diffq
                        if(diffg[l]>prevg[l]):
                            prevg = diffg
    

            svunq1,svung1 = prevq,prevg
            
            qsvunc = higher_quark_SV.Clone()
            gsvunc = higher_quark_SV.Clone()

            for j in range(1,higher_quark_SV.GetNbinsX()+1):
                    qsvunc.SetBinContent(j,svunq1[j-1])
                    gsvunc.SetBinContent(j,svung1[j-1])

                    sigma_tot_q[j-1][7] = svunq1[j-1]
                    sigma_tot_g[j-1][7] = svung1[j-1]
            qsvn = qsvunc.Clone()
            gsvn = gsvunc.Clone()

            qsvn.Scale(-1)
            gsvn.Scale(-1)
            """###
            #total uncertainty
            q_sigma_tot = higher_quark.Clone("")
            g_sigma_tot = higher_gluon.Clone("")

            for j in range(1,higher_quark_SV.GetNbinsX()+1):
                    qsvunc.SetBinContent(j,svunq1[j-1])
                    gsvunc.SetBinContent(j,svung1[j-1])

                    sigma_tot_q[j-1][7] = svunq1[j-1]
                    sigma_tot_g[j-1][7] = svung1[j-1]
            qsvn = qsvunc.Clone()
            gsvn = gsvunc.Clone()

            qsvn.Scale(-1)
            gsvn.Scale(-1)

            #total uncertainty
            q_sigma_tot = higher_quark.Clone("")
            g_sigma_tot = higher_gluon.Clone("")

            for j in range (0, quark.GetNbinsX()):
                    #print(sigma_tot_q[j][5],sigma_tot_g[j][5])
                    a = sigma_tot_q[j][0] #statistical
                    b = sigma_tot_q[j][1] # mc closure
                    c = sigma_tot_q[j][2] #reweight factor
                    d = sigma_tot_q[j][3] #PDF
                    e = sigma_tot_q[j][4] #Showering
                    f = sigma_tot_q[j][5] #Hadronization
                    g = sigma_tot_q[j][6] #ME
                    h = sigma_tot_q[j][7] #Scale variation
                    #print(a,b,c,d,e,f,g,h)
                    if var == "ntrk":
                        sigma_q_tot = np.sqrt((a**2)+(b**2)+(c**2)+(d**2)+(e**2)+(f**2)+(g**2)+(h**2))
                    if var == "bdt":
                        sigma_q_tot = np.sqrt((a**2)+(b**2)+(d**2)+(e**2)+(f**2)+(g**2)+(h**2))
                        
                    a = sigma_tot_g[j][0]
                    b = sigma_tot_g[j][1]
                    c = sigma_tot_g[j][2]
                    d = sigma_tot_g[j][3]
                    e = sigma_tot_g[j][4]
                    f = sigma_tot_g[j][5]
                    g = sigma_tot_g[j][6]
                    h = sigma_tot_g[j][7]
                    if var == "ntrk":
                        sigma_g_tot = np.sqrt((a**2)+(b**2)+(c**2)+(d**2)+(e**2)+(f**2)+(g**2)+(h**2))
                    if var == "bdt":
                        sigma_g_tot = np.sqrt((a**2)+(b**2)+(d**2)+(e**2)+(f**2)+(g**2)+(h**2))
                    q_sigma_tot.SetBinContent(j+1,sigma_q_tot)
                    g_sigma_tot.SetBinContent(j+1,sigma_g_tot)

                    #print("statistical: "+str(100*a)+" , MC Closure: "+str(100*b)+" , Showering: "+str(100*c)+" , PDF: "+str(100*d))

            q_sigma_tot.Scale(100)
            g_sigma_tot.Scale(100)

            q_sigma_tot_n = q_sigma_tot.Clone("")
            g_sigma_tot_n = g_sigma_tot.Clone("")
            q_sigma_tot_n.Scale(-1)
            g_sigma_tot_n.Scale(-1)
            
            quark_pdf.Scale(100)
            gluon_pdf.Scale(100)
            quark_pdf_negative.Scale(100)
            gluon_pdf_negative.Scale(100)
            
            quark_strap.Scale(100)
            gluon_strap.Scale(100)
            quark_negative.Scale(100)
            gluon_negative.Scale(100)

            quark_use.Scale(100)
            gluon_use.Scale(100)
            quarkMC_negative.Scale(100)
            gluonMC_negative.Scale(100)
 
            q_show_unc.Scale(100)
            g_show_unc.Scale(100)
            q_show_uncn.Scale(100)
            g_show_uncn.Scale(100)

            qhadunc.Scale(100)
            ghadunc.Scale(100)
            qhadn.Scale(100)
            ghadn.Scale(100)

            qmeunc.Scale(100)
            gmeunc.Scale(100)
            qmen.Scale(100)
            gmen.Scale(100)

            qsvunc.Scale(100)
            gsvunc.Scale(100)
            qsvn.Scale(100)
            gsvn.Scale(100)

            quark_RF.Scale(100)
            gluon_RF.Scale(100)
            quarkRF_negative.Scale(100)
            gluonRF_negative.Scale(100)
            
            ## below just do the ploting
            gPad.SetLeftMargin(0.15)
            gPad.SetTopMargin(0.05)
            gPad.SetBottomMargin(0.15)
            gPad.SetRightMargin(0.20)



            gStyle.SetOptStat(0)
            ######################## for ratio plo

            quark_strap.GetYaxis().SetRangeUser(-50,50)
            quark_strap.SetLineColor(2)
            quark_strap.SetLineStyle(2)
            #quark_strap.SetMarkerColor(8)
            #quark_strap.SetMarkerSize(0.8)
            quark_negative.SetLineColor(2)
            quark_negative.SetLineStyle(2)
            #quark_negative.SetMarkerSize(0.8)
            #quark_negative.SetMarkerColor(8)

            quark_use.SetLineColor(30)
            quark_use.SetLineStyle(2)
            #quark_use.SetMarkerColor(2)
            #quark_use.SetMarkerSize(0.8)
            quarkMC_negative.SetLineColor(30)
            quarkMC_negative.SetLineStyle(2)
            #quarkMC_negative.SetMarkerColor(2)
            #quarkMC_negative.SetMarkerSize(0.8)
            
            quark_pdf.SetLineColor(28)
            quark_pdf.SetLineStyle(2)
            quark_pdf_negative.SetLineColor(28)
            quark_pdf_negative.SetLineStyle(2)
            
            q_show_unc.SetLineColor(1)
            q_show_unc.SetLineStyle(2)
            q_show_uncn.SetLineColor(1)
            q_show_uncn.SetLineStyle(2)

            qhadunc.SetLineColor(9)
            qhadunc.SetLineStyle(2)
            qhadn.SetLineColor(9)
            qhadn.SetLineStyle(2)

            qmeunc.SetLineColor(7)
            qmeunc.SetLineStyle(2)
            qmen.SetLineColor(7)
            qmen.SetLineStyle(2)


            qsvunc.SetLineColor(6)
            qsvunc.SetLineStyle(2)
            qsvn.SetLineColor(6)
            qsvn.SetLineStyle(2)
            
            quark_RF.SetLineColor(3)
            quark_RF.SetLineStyle(2)
            quarkRF_negative.SetLineColor(3)
            quarkRF_negative.SetLineStyle(2)
            
            q_sigma_tot.SetLineColor(4)
            q_sigma_tot.SetLineStyle(1)
            q_sigma_tot.SetLineWidth(2)
            q_sigma_tot_n.SetLineColor(4)
            q_sigma_tot_n.SetLineStyle(1)
            q_sigma_tot_n.SetLineWidth(2)

            quark_strap.GetYaxis().SetTitle("Uncertainty (%)")

            quark_strap.Draw("HIST")
            quark_negative.Draw("HIST same")
            quark_use.Draw("HIST same")
            quarkMC_negative.Draw("HIST same")
            #quark_show_use.Draw("HIST same")
            #quark_show_negative.Draw("HIST same")
            quark_pdf.Draw("HIST same")
            quark_pdf_negative.Draw("HIST same")
            q_sigma_tot.Draw("HIST same")
            q_sigma_tot_n.Draw("HIST same")
            q_show_unc.Draw("HIST same")
            q_show_uncn.Draw("HIST same")
            qhadunc.Draw("hist same")
            qhadn.Draw("hist same")
            qmeunc.Draw("hist same")
            qmen.Draw("hist same")
            qsvunc.Draw("hist same")
            qsvn.Draw("hist same")
            if var == "ntrk":
                quark_RF.Draw("hist same")
                quarkRF_negative.Draw("hist same")
            leg = TLegend(0.80,0.5,0.999,0.9)##0.6,0.5,0.9,0.7
            leg.SetTextSize(0.022)
            
            leg.SetTextFont(42)
            leg.SetFillColor(0)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.SetNColumns(1)
            leg.AddEntry(quark_strap,"Statistical","l")
            leg.AddEntry(quark_use,"MC Closure","l")
            #leg.AddEntry(quark_show_use,"Showering","l")
            leg.AddEntry(quark_pdf,"PDF","l")
            leg.AddEntry(q_show_unc,"Showering","l")
            leg.AddEntry(qhadunc,"Hadronization","l")
            leg.AddEntry(qmeunc,"Matrix Element","l")
            leg.AddEntry(qsvunc,"Scale Variation","l")
            if var == "ntrk":
                leg.AddEntry(quark_RF,"Reweight Factor","l")
            leg.AddEntry(q_sigma_tot,"Total","l")

            myText(0.18,0.9,"#it{#bf{#scale[1.8]{#bf{ATLAS} Internal}}}")

            leg.Draw()

            myText(0.18,0.86,"#bf{#scale[1.5]{#sqrt{s} = 13 TeV}}")
            myText(0.18,0.82,"#bf{#scale[1.5]{pT range: "+str(min)+" - "+str(max)+" GeV}}")
            myText(0.18,0.78,"#bf{#scale[1.5]{Quark jet}}")

            if(inputvar == "ntrk"):
                line = TLine(0.,0,60,0)
                quark_strap.GetXaxis().SetTitle("n_{Track}")
            if(inputvar == "bdt"):
                line = TLine(-0.8,0,0.7,0)
                quark_strap.GetXaxis().SetTitle("BDT score")
#        line = TLine(0.,1,0.4,1)

#        quark_ratio.Draw()
            line.Draw("same")
            #c1.Print("./plots_bdt/quark_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"_fc.pdf")
            c1.Print("./plots_"+var+"/quark_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+".pdf")


            gluon_strap.GetYaxis().SetTitle("Uncertainty (%)")
            gluon_strap.GetYaxis().SetRangeUser(-50,50)


            gluon_strap.SetLineColor(2)
            gluon_strap.SetLineStyle(2)
            #gluon_strap.SetMarkerColor(2)
            #gluon_strap.SetMarkerSize(0.8)
            gluon_negative.SetLineColor(2)
            gluon_negative.SetLineStyle(2)
            #gluon_negative.SetMarkerColor(2)
            #gluon_negative.SetMarkerSize(0.8)


            gluon_use.SetLineColor(30)
            gluon_use.SetLineStyle(2)
            #gluon_use.SetMarkerColor(30)
            #gluon_use.SetMarkerSize(0.8)
            gluonMC_negative.SetLineColor(30)
            gluonMC_negative.SetLineStyle(2)
            #gluonMC_negative.SetMarkerColor(30)
            #gluonMC_negative.SetMarkerSize(0.8)

            #gluon_show_use.SetLineColor(6)
            #gluon_show_use.SetLineStyle(2)
            #gluon_show_negative.SetLineColor(6)
            #gluon_show_negative.SetLineStyle(2)

            gluon_pdf.SetLineColor(28)
            gluon_pdf.SetLineStyle(2)
            gluon_pdf_negative.SetLineColor(28)
            gluon_pdf_negative.SetLineStyle(2)

            g_show_unc.SetLineColor(1)
            g_show_unc.SetLineStyle(2)
            g_show_uncn.SetLineColor(1)
            g_show_uncn.SetLineStyle(2)

            ghadunc.SetLineColor(9)
            ghadunc.SetLineStyle(2)
            ghadn.SetLineColor(9)
            ghadn.SetLineStyle(2)

            gmeunc.SetLineColor(7)
            gmeunc.SetLineStyle(2)
            gmen.SetLineColor(7)
            gmen.SetLineStyle(2)

            gsvunc.SetLineColor(6)
            gsvunc.SetLineStyle(2)
            gsvn.SetLineColor(6)
            gsvn.SetLineStyle(2)

            gluon_RF.SetLineColor(3)
            gluon_RF.SetLineStyle(2)
            gluonRF_negative.SetLineColor(3)
            gluonRF_negative.SetLineStyle(2)
            
            g_sigma_tot.SetLineColor(4)
            g_sigma_tot.SetLineStyle(1)
            g_sigma_tot.SetLineWidth(2)
            g_sigma_tot_n.SetLineColor(4)
            g_sigma_tot_n.SetLineStyle(1)
            g_sigma_tot_n.SetLineWidth(2)

            gluon_strap.Draw("HIST")
            gluon_negative.Draw("HIST same")
            gluon_use.Draw("HIST same")
            gluonMC_negative.Draw("HIST same")
            gluon_pdf.Draw("HIST same")
            gluon_pdf_negative.Draw("HIST same")
            g_show_unc.Draw("HIST same")
            g_show_uncn.Draw("HIST same")
            g_sigma_tot.Draw("HIST same")
            g_sigma_tot_n.Draw("HIST same")
            ghadunc.Draw("HIST same")
            ghadn.Draw("hist same")
            gmeunc.Draw("hist same")
            gmen.Draw("hist same")
            gsvunc.Draw("hist same")
            gsvn.Draw("hist same")
            gluon_RF.Draw("hist same")
            if var == "ntrk":
                gluon_RF.Draw("hist same")
                gluonRF_negative.Draw("hist same")

            leg = TLegend(0.80,0.5,0.999,0.9)##0.6,0.5,0.9,0.7
            leg.SetTextSize(0.022)
            leg.SetTextFont(42)
            leg.SetFillColor(0)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.SetNColumns(1)
            leg.AddEntry(gluon_strap,"Statistical","l")
            leg.AddEntry(gluon_use,"MC Closure","l")
            leg.AddEntry(gluon_pdf,"PDF","l")
            leg.AddEntry(g_show_unc,"Showering","l")
            leg.AddEntry(ghadunc,"Hadronization","l")
            leg.AddEntry(gmeunc,"Matrix Element","l")
            leg.AddEntry(gsvunc,"Scale Variation","l")
            if var == "ntrk":
                leg.AddEntry(gluon_RF,"Reweight Factor","l")
            leg.AddEntry(g_sigma_tot,"Total","l")

            myText(0.18,0.9,"#it{#bf{#scale[1.8]{#bf{ATLAS} Internal}}}")

            leg.Draw()

            myText(0.18,0.86,"#bf{#scale[1.5]{#sqrt{s} = 13 TeV}}")
            myText(0.18,0.82,"#bf{#scale[1.5]{pT range: "+str(min)+" - "+str(max)+" GeV}}")
            myText(0.18,0.78,"#bf{#scale[1.5]{Gluon jet}}")

            if(inputvar == "ntrk"):
                line = TLine(0.,0,60,0)
                gluon_strap.GetXaxis().SetTitle("n_{Track}")
            if(inputvar == "bdt"):
                line = TLine(-0.8,0,0.7,0)
                gluon_strap.GetXaxis().SetTitle("BDT score")

#        bot.cd()
#        gluon_ratio.Draw()
            line.Draw("same")
            c1.Print("./plots_"+var+"/gluon_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+".pdf")

            #ufile = open("ufile"+str(min)+str(var)+".txt","w+")
            #for j in range(0,quark.GetNbinsX()):
            #        ufile.write("q,"+str(ebar_q[j][0])+","+str(ebar_q[j][1])+","+str(ebar_q[j][2])+","+str(ebar_q[j][3])+","+str(ebar_q[j][4])+","+str(ebar_q[j][5])+","+str(ebar_q[j][6])+","+str(ebar_q[j][7])+",")
            #        ufile.write("g,"+str(ebar_g[j][0])+","+str(ebar_g[j][1])+","+str(ebar_g[j][2])+","+str(ebar_g[j][3])+","+str(ebar_g[j][4])+","+str(ebar_g[j][5])+","+str(ebar_g[j][6])+","+str(ebar_g[j][7])+"\n")


            #This next part will plot each uncertainty separately in each pt bin.
            quark_strap.GetYaxis().SetRangeUser(-50,50)
            quark_use.GetYaxis().SetRangeUser(-50,50)
            q_show_unc.GetYaxis().SetRangeUser(-50,50)
            quark_pdf.GetYaxis().SetRangeUser(-50,50)
            qmeunc.GetYaxis().SetRangeUser(-50,50)
            qhadunc.GetYaxis().SetRangeUser(-50,50)
            qsvunc.GetYaxis().SetRangeUser(-50,50)
            quark_RF.GetYaxis().SetRangeUser(-50,50)
            
            gluon_strap.GetYaxis().SetRangeUser(-50,50)
            gluon_use.GetYaxis().SetRangeUser(-50,50)
            g_show_unc.GetYaxis().SetRangeUser(-50,50)
            gluon_pdf.GetYaxis().SetRangeUser(-50,50)
            gmeunc.GetYaxis().SetRangeUser(-50,50)
            ghadunc.GetYaxis().SetRangeUser(-50,50)
            gsvunc.GetYaxis().SetRangeUser(-50,50)
            gluon_RF.GetYaxis().SetRangeUser(-50,50)

        
            quark_strap.Draw("HIST")
            quark_negative.Draw("HIST same")
            quark_use.Draw("HIST same")
            quarkMC_negative.Draw("HIST same")
            quark_pdf.Draw("HIST same")
            quark_pdf_negative.Draw("HIST same")
            q_sigma_tot.Draw("HIST same")
            q_sigma_tot_n.Draw("HIST same")
            q_show_unc.Draw("HIST same")
            q_show_uncn.Draw("HIST same")
            qhadunc.Draw("hist same")
            qhadn.Draw("hist same")
            qmeunc.Draw("hist same")
            qmen.Draw("hist same")
            qsvunc.Draw("hist same")
            qsvn.Draw("hist same")
        


            quark_strap.Draw("HIST")
            quark_negative.Draw("HIST same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/quark_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"stat.pdf")

            quark_use.Draw("HIST")
            quarkMC_negative.Draw("HIST same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/quark_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"mc.pdf")

            q_show_unc.Draw("hist")
            q_show_uncn.Draw("hist same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/quark_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"shower.pdf")

            quark_pdf.Draw("HIST")
            quark_pdf_negative.Draw("HIST same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/quark_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"pdf.pdf")

            qhadunc.Draw("hist")
            qhadn.Draw("hist same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/quark_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"had.pdf")

            qmeunc.Draw("hist")
            qmen.Draw("hist same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/quark_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"me.pdf")

            qsvunc.Draw("hist")
            qsvn.Draw("hist same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/quark_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"sv.pdf")

            quark_RF.Draw("hist")
            quarkRF_negative.Draw("hist same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/quark_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"rf.pdf")
            
            gluon_strap.Draw("HIST")
            gluon_negative.Draw("HIST same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/gluon_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"stat.pdf")

            gluon_use.Draw("HIST")
            gluonMC_negative.Draw("HIST same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/gluon_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"mc.pdf")

            g_show_unc.Draw("hist")
            g_show_uncn.Draw("hist same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/gluon_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"shower.pdf")

            gluon_pdf.Draw("HIST")
            gluon_pdf_negative.Draw("HIST same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/gluon_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"pdf.pdf")

            ghadunc.Draw("hist")
            ghadn.Draw("hist same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/gluon_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"had.pdf")

            gmeunc.Draw("hist")
            gmen.Draw("hist same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/gluon_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"me.pdf")

            gsvunc.Draw("hist")
            gsvn.Draw("hist same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/guon_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"sv.pdf")

            
            gluon_RF.Draw("hist")
            gluonRF_negative.Draw("hist same")
            line.Draw("same")
            c1.Print("./plots_"+var+"/gluon_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"rf.pdf")



