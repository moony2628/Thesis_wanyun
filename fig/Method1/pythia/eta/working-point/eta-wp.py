from ROOT import *
import numpy as np
from prettytable import PrettyTable
import sys
import uncertainties as unc # propagate uncertainties
from uncertainties import unumpy as unp # array operations for type ufloat

# this script draw quark efficiency, gluon rejection , scale factor for given quark working point e.g. 0.6, 
# pdf is not included in the plot (can be added in line 684 and 686) because nan exist in old pdf sample 

wpoint_array = [float(sys.argv[1])]
var = sys.argv[2]  #change the var name according to the inputvar you want to read
inputvar = var  #by setting it as bdt (or ntrk,width,c1..), it will read the corresponding histogram, but remember to change the TLine range according to X-axis of different variable, one can check it by browsing the histograms in root file.

bins = np.array([500.,600.,800.,1000.,1200.,1500.,2000.])
bin = np.array([500,600,800,1000,1200,1500,2000])

#eta_bin = ["0-0.5_0.5-1","0.5-1_1-2.1"]
eta_bin = ["0-2.1_0-2.1","0-0.5_0.5-2.1"]



if var == "bdt":
    doreweight = 0

if var == "ntrk":
    doreweight = 0
def error_check(hist):
    for i in range(1,hist.GetNbinsX()+1):
        if hist.GetBinContent(i)>1:
            hist.SetBinContent(i,0)
    return(hist)


def unc_array_err(hist,err):
    value = np.zeros(hist.GetNbinsX())
    error = np.zeros(hist.GetNbinsX())
    for j in range(1,hist.GetNbinsX()+1):
        value[j-1] = hist.GetBinContent(j)
        error[j-1] = np.sqrt(err.GetBinContent(j))
    result = unp.uarray(value,error)
    return(result)

def myText(x,y,text, color = 1):
    l = TLatex()
    l.SetTextSize(0.025)
    l.SetNDC()
    l.SetTextColor(color)
    l.DrawLatex(x,y,text)
    pass
# input the mc sample(must be unnormalized)

#inclusive fraction calculation
def fraction(lower_quark,lower_gluon,higher_quark,higher_gluon):
    ToT_Fq2 = 0.
    ToT_Fg2 = 0.
    ToT_Cq2 =0.
    ToT_Cg2 = 0.

    for j in range(1,lower_quark.GetNbinsX()+1):
        ToT_Fq2+=higher_quark.GetBinContent(j)
        ToT_Cq2+=lower_quark.GetBinContent(j)
        ToT_Fg2+=higher_gluon.GetBinContent(j)
        ToT_Cg2+=lower_gluon.GetBinContent(j)

    # calculate the fraction of each sample 
    if ((ToT_Fg2+ToT_Fq2) != 0):
        fg=ToT_Fg2/(ToT_Fg2+ToT_Fq2)
        cg=ToT_Cg2/(ToT_Cq2+ToT_Cg2)

    fq=1.-fg
    cq=1.-cg
    return(fg,cg,fq,cq)


def mc_matrixmethod(lower_quark_input,lower_gluon_input,higher_quark_input,higher_gluon_input,fg,cg,fq,cq,higher_mc_input,lower_mc_input):
        lower_quark = lower_quark_input.Clone()
        lower_gluon = lower_gluon_input.Clone()
        higher_quark = higher_quark_input.Clone()
        higher_gluon = higher_gluon_input.Clone()
        lower_mc = lower_mc_input.Clone()
        higher_mc = higher_mc_input.Clone()
        if (doreweight=="Quark"):
                for i in range(1,higher_quark.GetNbinsX()+1):
                        if (lower_quark.GetBinContent(i) > 0 and lower_gluon.GetBinContent(i) > 0):
                                #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                factor_gluon = higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                factor_quark = higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i)
                                lower_quark.SetBinContent(i,lower_quark.GetBinContent(i)*factor_quark)
                                lower_gluon.SetBinContent(i,lower_gluon.GetBinContent(i)*factor_quark)
                                lower_mc.SetBinContent(i,lower_mc.GetBinContent(i)*factor_quark)


            
        if (doreweight=="Gluon"):
                for i in range(1,higher_quark.GetNbinsX()+1):
                        if (lower_quark.GetBinContent(i) > 0 and lower_gluon.GetBinContent(i) > 0):
                                #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                factor_gluon = higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                factor_quark = higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i)
                                lower_quark.SetBinContent(i,lower_quark.GetBinContent(i)*factor_gluon)
                                lower_gluon.SetBinContent(i,lower_gluon.GetBinContent(i)*factor_gluon)
                                lower_mc.SetBinContent(i,lower_mc.GetBinContent(i)*factor_gluon)

        
        if (lower_quark.Integral() != 0):
            lower_quark.Scale(1./lower_quark.Integral())
        if(lower_gluon.Integral() != 0):
            lower_gluon.Scale(1./lower_gluon.Integral())
        if(higher_quark.Integral() != 0):
            higher_quark.Scale(1./higher_quark.Integral())
        if(higher_gluon.Integral() != 0):
            higher_gluon.Scale(1./higher_gluon.Integral()) 
        if(lower_mc.Integral() != 0):
                lower_mc.Scale(1./lower_mc.Integral())
        if(higher_mc.Integral() != 0):
                higher_mc.Scale(1./higher_mc.Integral())

        higher = higher_mc.Clone()
        lower = lower_mc.Clone()          


        #Now, let's solve.
        quark_extracted = higher_quark.Clone()
        gluon_extracted = higher_quark.Clone()

        #Matrix method here
        for i in range(1,higher.GetNbinsX()+1):
                F = higher.GetBinContent(i)
                C = lower.GetBinContent(i)
                if((cg*fq-fg*cq) != 0 ):
                        Q = -(C*fg-F*cg)/(cg*fq-fg*cq)
                        G = (C*fq-F*cq)/(cg*fq-fg*cq)
                        quark_extracted.SetBinContent(i,Q)
                        gluon_extracted.SetBinContent(i,G)
                        #print "   ",i,G,higher_gluon.GetBinContent(i),lower_gluon.GetBinContent(i)
                pass
            
        return(quark_extracted,gluon_extracted)




def data_matrixmethod(lower_quark_input,lower_gluon_input,higher_quark_input,higher_gluon_input,higher_input,lower_input,fg,cg,fq,cq):
        lower_quark = lower_quark_input.Clone()
        lower_gluon = lower_gluon_input.Clone()
        higher_quark = higher_quark_input.Clone()
        higher_gluon = higher_gluon_input.Clone()
        lower = lower_input.Clone()
        higher = higher_input.Clone()
        if (doreweight=="Quark"):
                for i in range(1,higher_quark.GetNbinsX()+1):
                        if (lower_quark.GetBinContent(i) > 0 and lower_gluon.GetBinContent(i) > 0):
                                #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                factor_quark = higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i)
                                lower_quark.SetBinContent(i,lower_quark.GetBinContent(i)*factor_quark)
                                lower_gluon.SetBinContent(i,lower_gluon.GetBinContent(i)*factor_quark)
                                #print("testflag4:",i,lower.GetBinContent(i),higher_quark.GetBinContent(i),lower_quark.GetBinContent(i),factor_quark)
                                lower.SetBinContent(i,lower.GetBinContent(i)*factor_quark)    
        if (doreweight=="Gluon"):
                for i in range(1,higher_quark.GetNbinsX()+1):
                        if (lower_quark.GetBinContent(i) > 0 and lower_gluon.GetBinContent(i) > 0):
                                #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                factor_gluon = higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                factor_quark = higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i)
                                lower_quark.SetBinContent(i,lower_quark.GetBinContent(i)*factor_gluon)
                                lower_gluon.SetBinContent(i,lower_gluon.GetBinContent(i)*factor_gluon)
                                lower.SetBinContent(i,lower.GetBinContent(i)*factor_gluon)                                
        if (lower_quark.Integral() != 0):
            lower_quark.Scale(1./lower_quark.Integral())
        if(lower_gluon.Integral() != 0):
            lower_gluon.Scale(1./lower_gluon.Integral())
        if(higher_quark.Integral() != 0):
            higher_quark.Scale(1./higher_quark.Integral())
        if(higher_gluon.Integral() != 0):
            higher_gluon.Scale(1./higher_gluon.Integral()) 
        if(lower.Integral() != 0):
            lower.Scale(1./lower.Integral()) 
        if(higher.Integral() != 0):
            higher.Scale(1./higher.Integral())         
        #Now, let's solve.
        quark_extracted = higher_quark.Clone()
        gluon_extracted = higher_quark.Clone()
        #if min == 500:
        #    for j in range(1,61):
        #        print("testflag3:",j,higher.GetBinContent(j),lower.GetBinContent(j))

        #Matrix method here
        for i in range(1,higher.GetNbinsX()+1):
                F = higher.GetBinContent(i)
                C = lower.GetBinContent(i)
                if((cg*fq-fg*cq) != 0 ):
                        Q = -(C*fg-F*cg)/(cg*fq-fg*cq)
                        G = (C*fq-F*cq)/(cg*fq-fg*cq)
                        quark_extracted.SetBinContent(i,Q)
                        gluon_extracted.SetBinContent(i,G)
                        #print "   ",i,G,higher_gluon.GetBinContent(i),lower_gluon.GetBinContent(i)
                pass
            
        return(quark_extracted,gluon_extracted)
    
    
def wp_bin(wpoint,quark_mc,gluon_mc,quark_data,gluon_data):
    quark_mc_cumsum = np.cumsum(quark_mc)
    gluon_mc_cumsum = np.cumsum(gluon_mc)
    quark_data_cumsum = np.cumsum(quark_data)
    gluon_data_cumsum = np.cumsum(gluon_data)
    #print(quark_mc_cumsum[-1],quark_data_cumsum[-1])
    
    mc_bin = np.abs(quark_mc_cumsum -  wpoint*quark_mc.Integral(1,quark_mc.GetNbinsX()+1)).argmin()
    data_bin = np.abs(quark_data_cumsum -  wpoint*quark_data.Integral(1,quark_data.GetNbinsX()+1)).argmin()
    sf_q = quark_data_cumsum[mc_bin]/quark_mc_cumsum[mc_bin]
    sf_g = (gluon_data_cumsum[-1] - gluon_data_cumsum[mc_bin])/(gluon_mc_cumsum[-1] - gluon_mc_cumsum[mc_bin])
    q_eff_mc = quark_mc_cumsum[mc_bin]/quark_mc_cumsum[-1]
    g_rej_mc = (gluon_mc_cumsum[-1]-gluon_mc_cumsum[mc_bin])/gluon_mc_cumsum[-1]
    q_eff_data = quark_data_cumsum[mc_bin]/quark_data_cumsum[-1]
    g_rej_data = (gluon_data_cumsum[-1]-gluon_data_cumsum[mc_bin])/gluon_data_cumsum[-1]
    #print(gluon_data_cumsum,gluon_mc_cumsum)
    #print(gluon_data_cumsum[data_bin],gluon_mc_cumsum[mc_bin])
    return(mc_bin,data_bin,sf_q,sf_g,q_eff_mc,g_rej_mc,q_eff_data,g_rej_data)


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


gStyle.SetOptStat(0)
c = TCanvas("","",500,500)
#gPad.SetTickx()
#gPad.SetTicky()


ntrackall = TFile("../roots/pythia_etas_all.root") #mc sample
ntrackall3  = TFile("../roots/data_etas_prescale.root") #data sample


histq = TH1F("histq","",6,bins)
histg = TH1F("histg","",6,bins)
histqsf = TH1F("histqsf","",6,bins)
histgsf = TH1F("histgsf","",6,bins)




for wpoint in wpoint_array:
    for e in range(len(eta_bin)):
        qeff = []
        grej = []
        qsf_array = []
        gsf_array = []
        unc_array_q = []
        unc_array_g = []
        unc_form_q = []  
        unc_form_g = []
        for i in range(0,6):
                unc_form_q_pt = []
                unc_form_g_pt = []
    #for i in range(13):	#for gamma+jet combined with dijet event, start from jet pT>0 GeV
    #        if(bin[i] != 800):
                min = bin[i]
                max = bin[i+1]
                print(min)
    
                higher_quark2 = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
                higher_gluon2 = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
                higher_data2 = ntrackall3.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Forward_Data_"+inputvar)
                lower_quark2 = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Central_Quark_"+inputvar)
                lower_gluon2 = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
                lower_data2 = ntrackall3.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Central_Data_"+inputvar)

                higher_quark = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Forward_Quark_"+inputvar)
                higher_gluon = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Forward_Gluon_"+inputvar)

                lower_quark = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Central_Quark_"+inputvar)
                lower_gluon = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Central_Gluon_"+inputvar)

                higher_data = ntrackall3.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Forward_Data_"+inputvar)
                lower_data = ntrackall3.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Central_Data_"+inputvar)
           
                """
                higher_quark_pythia = ntrackall4.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
                higher_quark2_pythia = ntrackall4.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar)
                higher_gluon_pythia = ntrackall4.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
                higher_gluon2_pythia = ntrackall4.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar)
                lower_quark_pythia = ntrackall4.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar)
                lower_quark2_pythia = ntrackall4.Get(str(min)+"_SubJet_Central_Quark_"+inputvar)
                lower_gluon_pythia = ntrackall4.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
                lower_gluon2_pythia = ntrackall4.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar)
                """
                #add leading and subleading jet from only dijet event together,
                #note that for gammajet+dijet event, we need to add leading jet from gammajet and leading jet from dijet sample together
                higher_data.Add(higher_data2)
                lower_data.Add(lower_data2)
                quark_data = higher_data.Clone()
                gluon_data = higher_data.Clone()
                higher_quark.Add(higher_quark2)
                higher_gluon.Add(higher_gluon2)
                lower_quark.Add(lower_quark2)
                lower_gluon.Add(lower_gluon2)
                higher_mc = higher_quark.Clone()
                higher_mc.Add(higher_gluon)
                lower_mc = lower_quark.Clone()
                lower_mc.Add(lower_gluon)


                #uncertainty propagation to get mc statistics uncertainty
                
                higher_quark2_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Forward_Quark_"+inputvar + "_err")
                higher_gluon2_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Forward_Gluon_"+inputvar + "_err")
                lower_quark2_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Central_Quark_"+inputvar + "_err")
                lower_gluon2_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Central_Gluon_"+inputvar + "_err")

                higher_quark_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Forward_Quark_"+inputvar + "_err")
                higher_gluon_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Forward_Gluon_"+inputvar + "_err")

                lower_quark_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Central_Quark_"+inputvar + "_err")
                lower_gluon_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Central_Gluon_"+inputvar + "_err")

                higher_quark_err.Add(higher_quark2_err)
                higher_gluon_err.Add(higher_gluon2_err)
                lower_quark_err.Add(lower_quark2_err)
                lower_gluon_err.Add(lower_gluon2_err)
        
                
                
                higher_mc_err = higher_quark_err.Clone()
                higher_mc_err.Add(higher_gluon_err)
                lower_mc_err = lower_quark_err.Clone()
                lower_mc_err.Add(lower_gluon_err)         
                
                higher_quark_unc = unc_array_err(higher_quark,higher_quark_err)
                higher_gluon_unc = unc_array_err(higher_gluon,higher_gluon_err)
                lower_quark_unc = unc_array_err(lower_quark,lower_quark_err)
                lower_gluon_unc = unc_array_err(lower_gluon,lower_gluon_err)
                higher_mc_unc = unc_array_err(higher_mc,higher_mc_err)
                lower_mc_unc = unc_array_err(lower_mc,lower_mc_err)
                
                ToT_Fq2_unc = higher_quark_unc.sum()
                ToT_Fg2_unc = higher_gluon_unc.sum()
        
                ToT_Cq2_unc = lower_quark_unc.sum()
                ToT_Cg2_unc = lower_gluon_unc.sum()
        
                # calculate the fraction of forward(higher) / central(lower) quark or gluon jet
                fg_unc=ToT_Fg2_unc/(ToT_Fg2_unc+ToT_Fq2_unc)
                cg_unc=ToT_Cg2_unc/(ToT_Cq2_unc+ToT_Cg2_unc)
                fq_unc=1.-fg_unc
                cq_unc=1.-cg_unc
        
                factor_quark_unc = lower_quark_unc
                factor_gluon_unc = lower_gluon_unc
        
                higher_quark_unc = higher_quark_unc/higher_quark_unc.sum()
                higher_gluon_unc = higher_gluon_unc/higher_gluon_unc.sum()
                lower_quark_unc = lower_quark_unc/lower_quark_unc.sum()
                lower_gluon_unc = lower_gluon_unc/lower_gluon_unc.sum()
    
                higher_mc_unc = higher_mc_unc/higher_mc_unc.sum()
                lower_mc_unc = lower_mc_unc/lower_mc_unc.sum()        
            
                if (doreweight=="Quark"):
                    for j in range(1,higher_quark.GetNbinsX()+1):
                            if (lower_quark.GetBinContent(j) > 0 and lower_gluon.GetBinContent(j) > 0):
                                    factor_gluon_unc[j-1] = higher_gluon_unc[j-1]/lower_gluon_unc[j-1]
                                    factor_quark_unc[j-1] = higher_quark_unc[j-1]/lower_quark_unc[j-1]
                            else:
                                    factor_gluon_unc[j-1] = unc.ufloat(1, 0)
                                    factor_quark_unc[j-1] = unc.ufloat(1, 0)
                    lower_quark_unc=lower_quark_unc*factor_quark_unc
                    lower_gluon_unc=lower_gluon_unc*factor_quark_unc
                    lower_mc_unc=lower_mc_unc*factor_quark_unc
    
                if (doreweight=="Gluon"):
                    for j in range(1,higher_quark.GetNbinsX()+1):
                            if (lower_quark.GetBinContent(j) > 0 and lower_gluon.GetBinContent(j) > 0):
                                    #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                    factor_gluon_unc[j-1] = higher_gluon_unc[j-1]/lower_gluon_unc[j-1]
                                    factor_quark_unc[j-1] = higher_quark_unc[j-1]/lower_quark_unc[j-1]
                            else:
                                    factor_gluon_unc[j-1] = unc.ufloat(1, 0)
                                    factor_quark_unc[j-1] = unc.ufloat(1, 0)
                    lower_quark_unc=lower_quark_unc*factor_gluon_unc
                    lower_gluon_unc=lower_gluon_unc*factor_gluon_unc
                    lower_mc_unc=lower_mc_unc*factor_quark_unc
                
                higher_quark_unc = higher_quark_unc/higher_quark_unc.sum()
                higher_gluon_unc = higher_gluon_unc/higher_gluon_unc.sum()
                lower_quark_unc = lower_quark_unc/lower_quark_unc.sum()
                lower_gluon_unc = lower_gluon_unc/lower_gluon_unc.sum()
    
                higher_mc_unc = higher_mc_unc/higher_mc_unc.sum()
                lower_mc_unc = lower_mc_unc/lower_mc_unc.sum()   
                
                higher = higher_mc.Clone()
                lower = lower_mc.Clone()
        
                higher_unc = higher_mc_unc
                lower_unc  = lower_mc_unc
        
                
                # matrix method here
                #Now, let's solve.
                F_unc = higher_unc
                C_unc = lower_unc
                Q_unc = -(C_unc*fg_unc-F_unc*cg_unc)/(cg_unc*fq_unc-fg_unc*cq_unc)
                G_unc = (C_unc*fq_unc-F_unc*cq_unc)/(cg_unc*fq_unc-fg_unc*cq_unc)
    

                """
                higher_quark_pythia.Add(higher_quark2_pythia)
                higher_gluon_pythia.Add(higher_gluon2_pythia)
                lower_quark_pythia.Add(lower_quark2_pythia)
                lower_gluon_pythia.Add(lower_gluon2_pythia)
                
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
    
                higher_quark_herang.Add(higher_quark2_herang)
                higher_gluon_herang.Add(higher_gluon2_herang)
                lower_quark_herang.Add(lower_quark2_herang)
                lower_gluon_herang.Add(lower_gluon2_herang)
                
                
                hqlund = sherpa_lund.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
                hqlund2 = sherpa_lund.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar)
                hglund = sherpa_lund.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
                hglund2 = sherpa_lund.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar)
                lqlund = sherpa_lund.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar)
                lqlund2 = sherpa_lund.Get(str(min)+"_SubJet_Central_Quark_"+inputvar)
                lglund = sherpa_lund.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
                lglund2 = sherpa_lund.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar)
    
                #print(lqlund.GetBinContent(58),lglund.GetBinContent(58),hqlund.GetBinContent(58),hglund.GetBinContent(58))    
                #print(lqlund.GetBinContent(59),lglund.GetBinContent(59),hqlund.GetBinContent(59),hglund.GetBinContent(59))    
                
            
        #Matrix element uncertainty: pythia - powheg+pythia
                
                higher_quark_pow = powpyt.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
                higher_quark2_pow = powpyt.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar)
                higher_gluon_pow = powpyt.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
                higher_gluon2_pow = powpyt.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar)
                lower_quark_pow = powpyt.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar)
                lower_quark2_pow = powpyt.Get(str(min)+"_SubJet_Central_Quark_"+inputvar)
                lower_gluon_pow = powpyt.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
                lower_gluon2_pow = powpyt.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar)
    
                higher_quark_pythia = ntrackall4.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
                higher_quark2_pythia = ntrackall4.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar)
                higher_gluon_pythia = ntrackall4.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
                higher_gluon2_pythia = ntrackall4.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar)
                lower_quark_pythia = ntrackall4.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar)
                lower_quark2_pythia = ntrackall4.Get(str(min)+"_SubJet_Central_Quark_"+inputvar)
                lower_gluon_pythia = ntrackall4.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
                lower_gluon2_pythia = ntrackall4.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar)
    
                higher_quark_pow.Add(higher_quark2_pow)
                higher_gluon_pow.Add(higher_gluon2_pow)
                lower_quark_pow.Add(lower_quark2_pow)
                lower_gluon_pow.Add(lower_gluon2_pow)
    
                higher_quark_pythia.Add(higher_quark2_pythia)
                higher_gluon_pythia.Add(higher_gluon2_pythia)
                lower_quark_pythia.Add(lower_quark2_pythia)
                lower_gluon_pythia.Add(lower_gluon2_pythia)
                """
                #uncertainty calculations
                #uncertainty lists, number-of-bins lists of 4 uncertainties.
          
                # do matrix method to extract the distribution of sherpa first
                fg,cg,fq,cq = fraction(lower_quark,lower_gluon,higher_quark,higher_gluon)
    
                # normalize the sherpa mc
                if (lower_quark.Integral() != 0):
                    lower_quark.Scale(1./lower_quark.Integral())
                if(lower_gluon.Integral() != 0):
                    lower_gluon.Scale(1./lower_gluon.Integral())
                if(higher_quark.Integral() != 0):
                    higher_quark.Scale(1./higher_quark.Integral())
                if(higher_gluon.Integral() != 0):
                    higher_gluon.Scale(1./higher_gluon.Integral())            
                if(lower_mc.Integral() != 0):
                    lower_mc.Scale(1./lower_mc.Integral())
                if(higher_mc.Integral() != 0):
                    higher_mc.Scale(1./higher_mc.Integral())                               
                      
                sherpa_extract_Q,sherpa_extract_G = mc_matrixmethod(lower_quark,lower_gluon,higher_quark,higher_gluon,fg,cg,fq,cq,higher_mc,lower_mc)
                quark_data = higher_data.Clone()
                gluon_data = higher_data.Clone()
                
                doreweight_normal = doreweight
                doreweight = "Quark"
                Quark_Factor_extract_Q,Quark_Factor_extract_G = mc_matrixmethod(lower_quark,lower_gluon,higher_quark,higher_gluon,fg,cg,fq,cq,higher_mc,lower_mc)
                doreweight = "Gluon"
                Gluon_Factor_extract_Q,Gluon_Factor_extract_G = mc_matrixmethod(lower_quark,lower_gluon,higher_quark,higher_gluon,fg,cg,fq,cq,higher_mc,lower_mc)
                doreweight = doreweight_normal
    
                higher_data_strap = higher_data.Clone()     #Set aside for statistical uncertainty
                lower_data_strap = lower_data.Clone()
                #statistical
                # do bootstrap(not normalized yet)
                # create lists to store bootstrapped values list of arrays of nstraps values
                nstraps = 5000
                SF_Qvals = []
                SF_Gvals = []
                Qvals = np.zeros(5000)
                Gvals = np.zeros(5000)
    
                higher_quark_bootstrap = higher_quark.Clone()
                higher_gluon_bootstrap = higher_gluon.Clone()
                lower_quark_bootstrap = lower_quark.Clone()
                lower_gluon_bootstrap = lower_gluon.Clone()

                
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
                                
                        if forward_data_strap.Integral() != 0:
                            forward_data_strap.Scale(1/forward_data_strap.Integral())
                        if central_data_strap.Integral() != 0:
                            central_data_strap.Scale(1/central_data_strap.Integral())                            
                                
                        bootstrap_extracted_Q = central_data_strap.Clone()
                        bootstrap_extracted_G = central_data_strap.Clone()
    
                        # get extracted data Q/G with sherpa sample
                        for j in range(1,higher_quark.GetNbinsX()+1):
                                F_data = forward_data_strap.GetBinContent(j)
                                C_data = central_data_strap.GetBinContent(j)
                                Q_data = -(C_data*fg-F_data*cg)/(cg*fq-fg*cq)
                                G_data = (C_data*fq-F_data*cq)/(cg*fq-fg*cq)
                                bootstrap_extracted_Q.SetBinContent(j,Q_data)
                                bootstrap_extracted_G.SetBinContent(j,G_data)
                        
                        mc_bin,data_bin,sf_q,sf_g,q_eff_mc,g_rej_mc,q_eff_data,g_rej_data = wp_bin(wpoint,sherpa_extract_Q,sherpa_extract_G,bootstrap_extracted_Q,bootstrap_extracted_G)
                        
                        SF_Qvals.append(q_eff_data/q_eff_mc)
                        SF_Gvals.append(g_rej_data/g_rej_mc)
                        
                #compute the uncertainty and plots
                quark_strap = quark_data.Clone()
                gluon_strap = gluon_data.Clone()
    
                
                #do matrix method on data (sherpa fraction)
                # first normalize it
                if (higher_data.Integral() != 0):
                    higher_data.Scale(1/higher_data.Integral())
                if (lower_data.Integral() != 0):
                    lower_data.Scale(1/lower_data.Integral())
                  
                extracted_data_sherpa_Q,extracted_data_sherpa_G = data_matrixmethod(lower_quark,lower_gluon,higher_quark,higher_gluon,higher_data,lower_data,fg,cg,fq,cq)            
                mc_bin,data_bin,sf_q,sf_g,q_eff_mc,g_rej_mc,q_eff_data,g_rej_data = wp_bin(wpoint,sherpa_extract_Q,sherpa_extract_G,extracted_data_sherpa_Q,extracted_data_sherpa_G)
                doreweight_normal = doreweight
                doreweight = "Quark"
                extracted_data_Quark_Factor_Q,extracted_data_Quark_Factor_G = data_matrixmethod(lower_quark,lower_gluon,higher_quark,higher_gluon,higher_data,lower_data,fg,cg,fq,cq)            
                doreweight = "Gluon"
                extracted_data_Gluon_Factor_Q,extracted_data_Gluon_Factor_G = data_matrixmethod(lower_quark,lower_gluon,higher_quark,higher_gluon,higher_data,lower_data,fg,cg,fq,cq)            
                doreweight = doreweight_normal
                
                SF_Qvals.sort()
                SF_Gvals.sort()
                Q = np.median(Qvals)
                G = np.median(Gvals)
                
                sigmaQ = .5*(SF_Qvals[int(.84*len(SF_Qvals))] - SF_Qvals[int(.16*len(SF_Qvals))])
                sigmaG = .5*(SF_Gvals[int(.84*len(SF_Gvals))] - SF_Gvals[int(.16*len(SF_Gvals))])
    
                Q_eff_unc = unp.std_devs(Q_unc[0:mc_bin].sum())
                G_rej_unc = unp.std_devs(G_unc[mc_bin:-1].sum())
                #print("statistical: q = "+str(sigmaQ)+" | g = "+str(sigmaG))                  
                
                sigmaQ = np.abs(sigmaQ)
                sigmaG = np.abs(sigmaG)
                                
                ##Reweight Factor Uncertainty
                mc_bin_QF,data_bin_QF,sf_q_QF,sf_g_QF,q_eff_mc_QF,g_rej_mc_QF,q_eff_data_QF,g_rej_data_QF = wp_bin(wpoint,Quark_Factor_extract_Q,Quark_Factor_extract_G,extracted_data_Quark_Factor_Q,extracted_data_Quark_Factor_G)
                mc_bin_GF,data_bin_GF,sf_q_GF,sf_g_GF,q_eff_mc_GF,g_rej_mc_GF,q_eff_data_GF,g_rej_data_GF = wp_bin(wpoint,Gluon_Factor_extract_Q,Gluon_Factor_extract_G,extracted_data_Gluon_Factor_Q,extracted_data_Gluon_Factor_G)
                RF_unc_Q = abs(sf_q_QF-sf_q_GF)
                RF_unc_G = abs(sf_g_QF-sf_q_GF)
    
                #print(extracted_data_sherpa_Q.Integral(),sherpa_extract_Q.Integral())
                #Reweight Factor unc
             
                """
                #Matrix element uncertainty: pythia - powheg+pythia
                fg_pythia,cg_pythia,fq_pythia,cq_pythia = fraction(lower_quark_pythia,lower_gluon_pythia,higher_quark_pythia,higher_gluon_pythia)
                if (higher_quark_pythia.Integral() != 0):
                    higher_quark_pythia.Scale(1./higher_quark_pythia.Integral())
                if(higher_gluon_pythia.Integral() != 0):
                    higher_gluon_pythia.Scale(1./higher_gluon_pythia.Integral())
                if(lower_quark_pythia.Integral() != 0):
                    lower_quark_pythia.Scale(1./lower_quark_pythia.Integral())
                if(lower_gluon_pythia.Integral() != 0):
                    lower_gluon_pythia.Scale(1./lower_gluon_pythia.Integral())
                    
                
                pythia_extract_Q,pythia_extract_G = mc_matrixmethod(lower_quark_pythia,lower_gluon_pythia,higher_quark_pythia,higher_gluon_pythia,fg_pythia,cg_pythia,fq_pythia,cq_pythia)
                extracted_data_pythia_Q,extracted_data_pythia_G = data_matrixmethod(lower_quark_pythia,lower_gluon_pythia,higher_quark_pythia,higher_gluon_pythia,higher_data,lower_data,fg_pythia,cg_pythia,fq_pythia,cq_pythia)
                
                
                fg_pow,cg_pow,fq_pow,cq_pow = fraction(lower_quark_pow,lower_gluon_pow,higher_quark_pow,higher_gluon_pow)
                if (higher_quark_pow.Integral() != 0):
                    higher_quark_pow.Scale(1./higher_quark_pow.Integral())
                if(higher_gluon_pow.Integral() != 0):
                    higher_gluon_pow.Scale(1./higher_gluon_pow.Integral())
                if(lower_quark_pow.Integral() != 0):
                    lower_quark_pow.Scale(1./lower_quark_pow.Integral())
                if(lower_gluon_pow.Integral() != 0):
                    lower_gluon_pow.Scale(1./lower_gluon_pow.Integral())
                    
                    
                pow_extract_Q,pow_extract_G = mc_matrixmethod(lower_quark_pow,lower_gluon_pow,higher_quark_pow,higher_gluon_pow,fg_pow,cg_pow,fq_pow,cq_pow)
                me_qerr,me_gerr = mc_error(pow_extract_Q,pow_extract_G,pythia_extract_Q,pythia_extract_G)  
                
                #hadronization, difference in sherpa
    
                hqlund.Add(hqlund2)
                hglund.Add(hglund2)
                lqlund.Add(lqlund2)
                lglund.Add(lglund2)
                higher_lund = hqlund.Clone()
                higher_lund.Add(hglund)
                lower_lund = lqlund.Clone()
                lower_lund.Add(lglund)
                
                fg_lund,cg_lund,fq_lund,cq_lund = fraction(lqlund,lglund,hqlund,hglund)
                
                if (hqlund.Integral() != 0):
                    hqlund.Scale(1./hqlund.Integral())
                if(hglund.Integral() != 0):
                    hglund.Scale(1./hglund.Integral())
                if(lqlund.Integral() != 0):
                    lqlund.Scale(1./lqlund.Integral())
                if(lglund.Integral() != 0):
                    lglund.Scale(1./lglund.Integral())
                if(higher_lund.Integral() != 0):
                    higher_lund.Scale(1./higher_lund.Integral())            
                if(lower_lund.Integral() != 0):
                    lower_lund.Scale(1./lower_lund.Integral())    
                    
                lower_quark_lund =  lqlund.Clone()
                lower_gluon_lund =  lglund.Clone()
                higher_quark_lund =  hqlund.Clone()
                higher_gluon_lund =  hglund.Clone()
                    
    
                lund_extract_Q,lund_extract_G = mc_matrixmethod(lower_quark_lund,lower_gluon_lund,higher_quark_lund,higher_gluon_lund,fg_lund,cg_lund,fq_lund,cq_lund,higher_lund,lower_lund)            
                extracted_data_lund_Q,extracted_data_lund_G = data_matrixmethod(lower_quark_lund,lower_gluon_lund,higher_quark_lund,higher_gluon_lund,higher_data,lower_data,fg_lund,cg_lund,fq_lund,cq_lund)            
    
                lund_extract_G = error_check(lund_extract_G)
                lund_extract_Q = error_check(lund_extract_Q)
                   
                
                #print("testflag:",lund_extract_G.GetBinContent(59),sherpa_extract_G.GetBinContent(59))
                mc_bin_lund,data_bin_lund,sf_q_lund,sf_g_lund,q_eff_lund,g_rej_lund,q_eff_data_lund,g_rej_data_lund = wp_bin(wpoint,lund_extract_Q,lund_extract_G,extracted_data_lund_Q,extracted_data_lund_G)
                had_unc_Q = abs(sf_q_lund-sf_q)
                had_unc_G = abs(sf_g_lund-sf_g)
    
                # Showering
                higher_herang = higher_quark_herang.Clone()
                higher_herang.Add(higher_gluon_herang)
                lower_herang = lower_quark_herang.Clone()
                lower_herang.Add(lower_gluon_herang)
                
                fg_herang,cg_herang,fq_herang,cq_herang = fraction(lower_quark_herang,lower_gluon_herang,higher_quark_herang,higher_gluon_herang)
                
                if (higher_quark_herang.Integral() != 0):
                    higher_quark_herang.Scale(1./higher_quark_herang.Integral())
                if(higher_gluon_herang.Integral() != 0):
                    higher_gluon_herang.Scale(1./higher_gluon_herang.Integral())
                if(lower_quark_herang.Integral() != 0):
                    lower_quark_herang.Scale(1./lower_quark_herang.Integral())
                if(lower_gluon_herang.Integral() != 0):
                    lower_gluon_herang.Scale(1./lower_gluon_herang.Integral())
                if(higher_herang.Integral() != 0):
                    higher_herang.Scale(1./higher_herang.Integral())            
                if(lower_herang.Integral() != 0):
                    lower_herang.Scale(1./lower_herang.Integral())    
                    
    
                herang_extract_Q,herang_extract_G = mc_matrixmethod(lower_quark_herang,lower_gluon_herang,higher_quark_herang,higher_gluon_herang,fg_herang,cg_herang,fq_herang,cq_herang,higher_herang,lower_herang)            
                extracted_data_herang_Q,extracted_data_herang_G = data_matrixmethod(lower_quark_herang,lower_gluon_herang,higher_quark_herang,higher_gluon_herang,higher_data,lower_data,fg_herang,cg_herang,fq_herang,cq_herang)            
    
                herang_extract_G = error_check(herang_extract_G)
                herang_extract_Q = error_check(herang_extract_Q)
                   
                
                mc_bin_herang,data_bin_herang,sf_q_herang,sf_g_herang,q_eff_herang,g_rej_herang,q_eff_data_herang,g_rej_data_herang = wp_bin(wpoint,herang_extract_Q,herang_extract_G,extracted_data_herang_Q,extracted_data_herang_G)
                
                higher_herdipo = higher_quark_herdipo.Clone()
                higher_herdipo.Add(higher_gluon_herdipo)
                lower_herdipo = lower_quark_herdipo.Clone()
                lower_herdipo.Add(lower_gluon_herdipo)
                
                fg_herdipo,cg_herdipo,fq_herdipo,cq_herdipo = fraction(lower_quark_herdipo,lower_gluon_herdipo,higher_quark_herdipo,higher_gluon_herdipo)
                
                if (higher_quark_herdipo.Integral() != 0):
                    higher_quark_herdipo.Scale(1./higher_quark_herdipo.Integral())
                if(higher_gluon_herdipo.Integral() != 0):
                    higher_gluon_herdipo.Scale(1./higher_gluon_herdipo.Integral())
                if(lower_quark_herdipo.Integral() != 0):
                    lower_quark_herdipo.Scale(1./lower_quark_herdipo.Integral())
                if(lower_gluon_herdipo.Integral() != 0):
                    lower_gluon_herdipo.Scale(1./lower_gluon_herdipo.Integral())
                if(higher_herdipo.Integral() != 0):
                    higher_herdipo.Scale(1./higher_herdipo.Integral())            
                if(lower_herdipo.Integral() != 0):
                    lower_herdipo.Scale(1./lower_herdipo.Integral())    
                    
    
                herdipo_extract_Q,herdipo_extract_G = mc_matrixmethod(lower_quark_herdipo,lower_gluon_herdipo,higher_quark_herdipo,higher_gluon_herdipo,fg_herdipo,cg_herdipo,fq_herdipo,cq_herdipo,higher_herdipo,lower_herdipo)            
                extracted_data_herdipo_Q,extracted_data_herdipo_G = data_matrixmethod(lower_quark_herdipo,lower_gluon_herdipo,higher_quark_herdipo,higher_gluon_herdipo,higher_data,lower_data,fg_herdipo,cg_herdipo,fq_herdipo,cq_herdipo)            
    
                herdipo_extract_G = error_check(herdipo_extract_G)
                herdipo_extract_Q = error_check(herdipo_extract_Q)
                
                
                mc_bin_herdipo,data_bin_herdipo,sf_q_herdipo,sf_g_herdipo,q_eff_herdipo,g_rej_herdipo,q_eff_data_herdipo,g_rej_data_herdipo = wp_bin(wpoint,herdipo_extract_Q,herdipo_extract_G,extracted_data_herdipo_Q,extracted_data_herdipo_G)
                
                show_unc_Q = abs(sf_q_herang-sf_q_herdipo)
                show_unc_G = abs(sf_g_herang-sf_g_herdipo)
                
                #pdf
                pdf_qerr = np.zeros(60)
                pdf_gerr = np.zeros(60)
                pdf_qvals = []
                pdf_gvals = []
    
                sf_pdf_qvals = np.zeros(53)
                sf_pdf_gvals = np.zeros(53)           
                for k in range(53):
                        higher_quark_pdf = ntrackall5.Get(str(min)+"_LeadingJet_Forward_Quark_"+str(k)+"_"+inputvar)
                        higher_quark1_pdf = ntrackall5.Get(str(min)+"_SubJet_Forward_Quark"+str(k)+"_"+inputvar)
                        lower_quark_pdf = ntrackall5.Get(str(min)+"_LeadingJet_Central_Quark_"+str(k)+"_"+inputvar)
                        lower_quark1_pdf = ntrackall5.Get(str(min)+"_SubJet_Central_Quark"+str(k)+"_"+inputvar)
                        higher_gluon_pdf = ntrackall5.Get(str(min)+"_LeadingJet_Forward_Gluon_"+str(k)+"_"+inputvar)
                        higher_gluon1_pdf = ntrackall5.Get(str(min)+"_SubJet_Forward_Gluon"+str(k)+"_"+inputvar)
                        lower_gluon_pdf = ntrackall5.Get(str(min)+"_LeadingJet_Central_Gluon_"+str(k)+"_"+inputvar)
                        lower_gluon1_pdf = ntrackall5.Get(str(min)+"_SubJet_Central_Gluon"+str(k)+"_"+inputvar)
                        higher_quark_pdf.SetName("test")
                        lower_quark_pdf.SetName("test")
                        higher_gluon_pdf.SetName("test")
                        lower_gluon_pdf.SetName("test")
                        higher_quark1_pdf.SetName("test")
                        lower_quark1_pdf.SetName("test")
                        higher_gluon1_pdf.SetName("test")
                        lower_gluon1_pdf.SetName("test")
    
                        higher_quark_pdf.Add(higher_quark1_pdf)
                        higher_gluon_pdf.Add(higher_gluon1_pdf)
                        lower_quark_pdf.Add(lower_quark1_pdf)
                        lower_gluon_pdf.Add(lower_gluon1_pdf)
                        higher_pdf = higher_quark_pdf.Clone()
                        higher_pdf.Add(higher_gluon_pdf)
                        lower_pdf = lower_quark_pdf.Clone()
                        lower_pdf.Add(lower_gluon_pdf)
                        
                        ToT_Fq2 = 0.
                        ToT_Fg2 = 0.
            
                        ToT_Cq2 = 0.
                        ToT_Cg2 = 0.
            
                        for j in range(1,lower_quark_pdf.GetNbinsX()+1):
                                #print(j,higher_gluon_pdf.GetBinContent(j))
    
                                ToT_Fq2+=higher_quark_pdf.GetBinContent(j)
                                ToT_Cq2+=lower_quark_pdf.GetBinContent(j)
                                ToT_Fg2+=higher_gluon_pdf.GetBinContent(j)
                                ToT_Cg2+=lower_gluon_pdf.GetBinContent(j)
            
                        fg_pdf=ToT_Fg2/(ToT_Fg2+ToT_Fq2)
                        cg_pdf=ToT_Cg2/(ToT_Cq2+ToT_Cg2)
                        fq_pdf=1.-fg
                        cq_pdf=1.-cg
                        
                        if(lower_quark_pdf.Integral() != 0):
                                lower_quark_pdf.Scale(1./lower_quark_pdf.Integral())
                        if(lower_gluon_pdf.Integral() != 0):
                                lower_gluon_pdf.Scale(1./lower_gluon_pdf.Integral())
                        if(higher_quark_pdf.Integral() != 0):
                                higher_quark_pdf.Scale(1./higher_quark_pdf.Integral())
                        if(higher_gluon_pdf.Integral() != 0):
                                higher_gluon_pdf.Scale(1./higher_gluon_pdf.Integral())
                        if(higher_pdf.Integral() != 0):
                                higher_pdf.Scale(1./higher_pdf.Integral())                            
                        if(lower_gluon_pdf.Integral() != 0):
                                lower_pdf.Scale(1./lower_pdf.Integral())
                                
                        #for i in range(1,higher_quark.GetNbinsX()+1):
                                #print(higher_pdf.GetBinContent(i))
                        pdf_extract_Q,pdf_extract_G = mc_matrixmethod(lower_quark_pdf,lower_gluon_pdf,higher_quark_pdf,higher_gluon_pdf,fg_pdf,cg_pdf,fq_pdf,cq_pdf,higher_pdf,lower_pdf)
                        extracted_data_pdf_Q,extracted_data_pdf_G = data_matrixmethod(lower_quark_pdf,lower_gluon_pdf,higher_quark_pdf,higher_gluon_pdf,higher_data,lower_data,fg_pdf,cg_pdf,fq_pdf,cq_pdf)
                        mc_bin_pdf,data_bin_pdf,sf_q_pdf,sf_g_pdf,q_eff_pdf,g_rej_pdf,q_eff_data_pdf,g_rej_data_pdf = wp_bin(wpoint,pdf_extract_Q,pdf_extract_G,extracted_data_pdf_Q,extracted_data_pdf_G)
                        sf_pdf_qvals[k] = sf_q_pdf
                        sf_pdf_gvals[k] = sf_g_pdf
    
                sf_pdf_qvals.sort()
                sf_pdf_gvals.sort()
                Q = np.median(sf_pdf_qvals)
                G = np.median(sf_pdf_gvals)
                 
                pdf_sigmaQ = abs(np.std(sf_pdf_qvals))
                pdf_sigmaG = abs(np.std(sf_pdf_gvals))
                pdf_unc_Q = pdf_sigmaQ
                pdf_unc_G = pdf_sigmaG
                
                ### scale variation ### 
                
                for k in range(7):
                    higher_quark_1 = sv.Get(str(min)+"_LeadingJet_Forward_Quark_"+str(k)+"_" + inputvar)
                    higher_quark2_1 = sv.Get(str(min)+"_SubJet_Forward_Quark_"+str(k)+"_" + inputvar)
                    higher_gluon_1 = sv.Get(str(min)+"_LeadingJet_Forward_Gluon_"+str(k)+"_" + inputvar)
                    higher_gluon2_1 = sv.Get(str(min)+"_SubJet_Forward_Gluon_"+str(k)+"_" + inputvar)
                    lower_quark_1 = sv.Get(str(min)+"_LeadingJet_Central_Quark_"+str(k)+"_" + inputvar)
                    lower_quark2_1 = sv.Get(str(min)+"_SubJet_Central_Quark_"+str(k)+"_" + inputvar)
                    lower_gluon_1 = sv.Get(str(min)+"_LeadingJet_Central_Gluon_"+str(k)+"_" + inputvar)
                    lower_gluon2_1 = sv.Get(str(min)+"_SubJet_Central_Gluon_"+str(k)+"_" + inputvar)
                    
                    higher_quark_1.Add(higher_quark2_1)
                    higher_gluon_1.Add(higher_gluon2_1)
                    lower_quark_1.Add(lower_quark2_1)
                    lower_gluon_1.Add(lower_gluon2_1)
                    
                    higher_1 = higher_quark_1.Clone()
                    higher_1.Add(higher_gluon_1)
                    lower_1 = lower_quark_1.Clone()
                    lower_1.Add(lower_gluon_1)
                    fg_1,cg_1,fq_1,cq_1 = fraction(lower_quark_1,lower_gluon_1,higher_quark_1,higher_gluon_1)
                    
                    if (higher_quark_1.Integral() != 0):
                        higher_quark_1.Scale(1./higher_quark_1.Integral())
                    if(higher_gluon_1.Integral() != 0):
                        higher_gluon_1.Scale(1./higher_gluon_1.Integral())
                    if(lower_quark_1.Integral() != 0):
                        lower_quark_1.Scale(1./lower_quark_1.Integral())
                    if(lower_gluon_1.Integral() != 0):
                        lower_gluon_1.Scale(1./lower_gluon_1.Integral())
                        
                    SV_extract_Q,SV_extract_G = mc_matrixmethod(lower_quark_1,lower_gluon_1,higher_quark_1,higher_gluon_1,fg_1,cg_1,fq_1,cq_1,higher_1,lower_1)
                    extracted_data_SV_Q,extracted_data_SV_G = data_matrixmethod(lower_quark_1,lower_gluon_1,higher_quark_1,higher_gluon_1,higher_data,lower_data,fg_1,cg_1,fq_1,cq_1)
                    mc_bin_SV,data_bin_SV,sf_q_SV,sf_g_SV,q_eff_SV,g_rej_SV,q_eff_data_SV,g_rej_data_SV = wp_bin(wpoint,SV_extract_Q,SV_extract_G,extracted_data_SV_Q,extracted_data_SV_G)
    
                    if k == 1:
                        SV_unc_Q = abs(sf_q_SV-sf_q)
                        SV_unc_G = abs(sf_g_SV-sf_g)
                        print(SV_unc_G)
                    if k > 1:
                        diffq = abs(sf_q_SV-sf_q)
                        diffg = abs(sf_g_SV-sf_g)
                        if(diffq>SV_unc_Q):
                            SV_unc_Q = diffq
                        if(diffg>SV_unc_G):
                            SV_unc_G = diffg
                """                        
                #mc closure
                mc_bin_closure,data_bin_closure,sf_q_closure,sf_g_closure,q_eff_closure,g_rej_closure,q_eff_data_closure,g_rej_data_closure = wp_bin(wpoint,higher_quark,higher_gluon,extracted_data_sherpa_Q,extracted_data_sherpa_G)
                closure_unc_Q = abs(sf_q_closure-sf_q)
                closure_unc_G = abs(sf_g_closure-sf_g)            
                
                #calculate the SF uncertainty
                qeff.append(q_eff_mc)
                grej.append(g_rej_mc)
                qsf_array.append(sf_q)
                gsf_array.append(sf_g)
                if var == "ntrk":
                    unc_q = np.sqrt(sigmaQ**2)
                    unc_g = np.sqrt(sigmaG**2)
                if var == "bdt":
                    unc_q = np.sqrt(sigmaQ**2)
                    unc_g = np.sqrt(sigmaG**2)
                histq.SetBinContent(i+1,q_eff_mc)
                histg.SetBinContent(i+1,g_rej_mc)
                histq.SetBinError(i+1,Q_eff_unc)
                histg.SetBinError(i+1,G_rej_unc) 
                histqsf.SetBinContent(i+1,sf_q)
                histgsf.SetBinContent(i+1,sf_g)
                histqsf.SetBinError(i+1,unc_q)
                histgsf.SetBinError(i+1,unc_g) 
                """
                systematic_qerr_hist = extracted_data_sherpa_Q.Clone()
                for j in range(1,61):
                    if systematic_qerr_hist.GetBinContent(j)!=0:
                        systematic_qerr_hist.SetBinContent(j,systematic_qerr[j-1]/extracted_data_sherpa_Q.GetBinContent(j))
                    else:
                        systematic_qerr_hist.SetBinContent(j,0)
                systematic_qerr_hist.SetMaximum(0.5)
                c1 = TCanvas("","",500,500)
                systematic_qerr_hist.Draw("hist")
                c1.Print("dijet-wp-unc" + str(bin[i]) + ".pdf")
            
                """
                pt_range = str(str(min)+"-"+str(max))
                if var == "ntrk":
                    unc_form_q_pt = [pt_range,round(sf_q,5),round(unc_q,5),round(sigmaQ,5),round(closure_unc_Q,5)]
                    unc_form_g_pt = [pt_range,round(sf_g,5),round(unc_g,5),round(sigmaG,5),round(closure_unc_G,5)]
                if var == "bdt":
                    unc_form_q_pt = [pt_range,round(sf_q,5),round(unc_q,5),round(sigmaQ,5),round(closure_unc_Q,5)]
                    unc_form_g_pt = [pt_range,round(sf_g,5),round(unc_g,5),round(sigmaG,5),round(closure_unc_G,5)]
                unc_form_q.append(unc_form_q_pt)
                unc_form_g.append(unc_form_g_pt)
    
                
        rmax= 4/3
        histq.SetMaximum(1.2)
        histq.SetMinimum(0.3)
        histq.GetXaxis().SetRangeUser(500.,2000.)
        histq.GetXaxis().SetTitle("Jet p_{T} (GeV)")
        histq.GetYaxis().SetTitle("Efficiency")
        histq.SetLabelSize(0.04,"Y")
        histq.SetTitleOffset(1.1,"Y")
        histq.SetTitleSize(0.04,"Y")
        
        histq.SetMarkerStyle(24)
        histq.SetLineColor(4)
        histq.SetMarkerColor(1)
        histq.SetMarkerSize(1)
        histq.SetFillColor(4)
        histq.SetFillStyle(3005)
        
        
        histg.SetLineColor(2)
        histg.SetMarkerStyle(32)
        histg.SetMarkerSize(1)
        histg.SetMarkerColor(1)
        histg.SetFillColor(2)
        histg.SetFillStyle(3005)
        
        histqsf.SetLineColor(4)
        histqsf.SetMarkerStyle(25)
        histqsf.SetMarkerSize(1)
        histqsf.SetMarkerColor(1)
        histqsf.SetLineStyle(7)
        
        histgsf.SetLineColor(2)
        histgsf.SetMarkerStyle(26)
        histgsf.SetMarkerSize(1)
        histgsf.SetMarkerColor(1)
        histgsf.SetLineStyle(7)
        
        scale = 1/rmax
        histqsf.Scale(scale)
        histgsf.Scale(scale)
        from ROOT import *
        leg = TLegend(0.6,0.7,0.84,0.9)

        leg.SetTextFont(42)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.AddEntry(histq,"Quark Efficiency","lp")
        leg.AddEntry(histqsf,"Quark Scale Factor","lp")
        leg.AddEntry(histg,"Gluon Rejection","lp")
        leg.AddEntry(histgsf,"Gluon Scale Factor","lp")
        
        histq.Draw("L P0 E2")
        histg.Draw("L P0 E2 same")
        histqsf.Draw("hist E0 L P0 same")
        histgsf.Draw("hist E0 L P0 same")
        leg.Draw("same")
        myText(0.16,0.84,"#it{#bf{#scale[1.8]{#bf{ATLAS} Internal}}}")
        myText(0.16,0.80,"#bf{#scale[1.5]{#sqrt{s} = 13 TeV ,139 fb^{-1}}}")
        #if var == "ntrk":
        #    myText(0.18,0.76,"#bf{#scale[1.5]{N_{track} Quark Efficiency" + str(100 * wpoint) +"%  }}")
        #if var == "bdt":
        #    myText(0.18,0.76,"#bf{#scale[1.5]{BDT Quark Efficiency" + str(100 * wpoint) +"%  }}")
        axis = TGaxis(2000.,0.3,2000.,1.2,0.4,1.6,510,"+L")

        axis.SetTitle("Scale factor")
        axis.Draw()
        
        c.Print(eta_bin[e]+"-dijet-wp-unc" + str(wpoint) +"-"+ str(var)+ "-rej.pdf")
        tbq = PrettyTable()
        tbq.field_names = ["Pt(GeV)","Quark Scale Factor","Uncertainty","Statistical","MC closure"]
        for i in range(6):
            tbq.add_row(unc_form_q[i])
        print(tbq)
        tbg = PrettyTable()
        tbg.field_names = ["Pt(GeV)","Gluon Scale Factor","Uncertainty","Statistical","MC closure"]
        for i in range(6):
            tbg.add_row(unc_form_g[i])
        print(tbg)
        foutput = open(str(eta_bin[e])+"_"+str(wpoint) +"-"+ str(var)+".txt", "w")
        foutput.write(str(eta_bin[e])+" "+str(wpoint) + " " + str(var) +" "+ str(doreweight) +  "\n")
        foutput.write(tbq.get_string())
        foutput.write("\n")
        foutput.write(tbg.get_string())
        foutput.close()
