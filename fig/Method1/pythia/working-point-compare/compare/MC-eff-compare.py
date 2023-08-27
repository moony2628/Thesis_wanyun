from ROOT import *
import numpy as np
from prettytable import PrettyTable
import sys
import uncertainties as unc # propagate uncertainties
from uncertainties import unumpy as unp # array operations for type ufloat

# this script draw quark efficiency, gluon rejection , scale factor for given quark working point e.g. 0.6, 
# pdf is not included in the plot (can be added in line 684 and 686) because nan exist in old pdf sample 
doreweight = "Quark"   #decide if we want to do the reweighting process

wpoint_array = [float(sys.argv[1])]

var = sys.argv[2]  #change the var name according to the inputvar you want to read

inputvar = var  #by setting it as bdt (or ntrk,width,c1..), it will read the corresponding histogram, but remember to change the TLine range according to X-axis of different variable, one can check it by browsing the histograms in root file.

bins = np.array([500.,600.,800.,1000.,1200.,1500.,2000.])
bin = np.array([500,600,800,1000,1200,1500,2000])

if var == "bdt":
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

def unc_data_err(hist):
    value = np.zeros(hist.GetNbinsX())
    error = np.zeros(hist.GetNbinsX())
    for j in range(1,hist.GetNbinsX()+1):
        value[j-1] = hist.GetBinContent(j)
        error[j-1] = hist.GetBinError(j)
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

#return the abs systematic uncertainty
def mc_error(mc_quark_1,mc_gluon_1,mc_quark_2,mc_gluon_2): 
    error_q = np.zeros(mc_quark_1.GetNbinsX())
    error_g = np.zeros(mc_gluon_1.GetNbinsX())

    for j in range(1,mc_gluon_1.GetNbinsX()+1):
                q = abs(mc_quark_1.GetBinContent(j) - mc_quark_2.GetBinContent(j))
                g = abs(mc_gluon_1.GetBinContent(j) - mc_gluon_2.GetBinContent(j))
                error_q[j-1] = q
                error_g[j-1] = g

    return error_q,error_g


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
                                factor_quark = higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i)
                                lower_quark.SetBinContent(i,lower_quark.GetBinContent(i)*factor_quark)
                                lower_gluon.SetBinContent(i,lower_gluon.GetBinContent(i)*factor_quark)
                                lower.SetBinContent(i,lower.GetBinContent(i)*factor_quark)    
        if (doreweight=="Gluon"):
                for i in range(1,higher_quark.GetNbinsX()+1):
                        if (lower_quark.GetBinContent(i) > 0 and lower_gluon.GetBinContent(i) > 0):
                                factor_gluon = higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
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


qeff = []
grej = []
qsf_array = []
gsf_array = []
unc_array_q = []
unc_array_g = []


mc_sample_1 = "pythia"
mc_sample_2 = "herang"

mc_sample_list = [TFile("../../newroot/dijet_pythia_nominal.root"),TFile("../../newroot/dijet_herang.root")] 
ntrackall3 = TFile("../../newroot/dijet_data_prescale.root")

histq_array = []  
histg_array = []

for wpoint in wpoint_array:
    for ntrackall4 in mc_sample_list:
        histq = TH1F(str(ntrackall4)+"q","",6,bins)
        histg = TH1F(str(ntrackall4)+"g","",6,bins)
        print(histq.GetName())
        for i in range(0,6):
            
            min = bin[i]
            max = bin[i+1]
            print(min)

            higher_quark = ntrackall4.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
            higher_quark2 = ntrackall4.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar)
            higher_gluon = ntrackall4.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
            higher_gluon2 = ntrackall4.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar)
            lower_quark = ntrackall4.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar)
            lower_quark2 = ntrackall4.Get(str(min)+"_SubJet_Central_Quark_"+inputvar)
            lower_gluon = ntrackall4.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
            lower_gluon2 = ntrackall4.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar)

            higher_data = ntrackall3.Get(str(min)+"_LeadingJet_Forward_Data_"+inputvar)
            higher_data2 = ntrackall3.Get(str(min)+"_SubJet_Forward_Data_"+inputvar)
            lower_data = ntrackall3.Get(str(min)+"_LeadingJet_Central_Data_"+inputvar)
            lower_data2 = ntrackall3.Get(str(min)+"_SubJet_Central_Data_"+inputvar)
            higher_data.SetName(str(ntrackall4)+"hdata")
            lower_data.SetName(str(ntrackall4)+"ldata")
            higher_data2.SetName(str(ntrackall4)+"hdata2")
            lower_data2.SetName(str(ntrackall4)+"ldata2")
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
            
            higher_mc.SetName(str(ntrackall4)+"hmc")
            lower_mc.SetName(str(ntrackall4)+"lmc")

            #uncertainty propagation to get mc statistics uncertainty
            
            higher_quark2_err = ntrackall4.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar + "_err")
            higher_gluon2_err = ntrackall4.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar + "_err")
            lower_quark2_err = ntrackall4.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar + "_err")
            lower_gluon2_err = ntrackall4.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar + "_err")
    
            higher_quark_err = ntrackall4.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar + "_err")
            higher_gluon_err = ntrackall4.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar + "_err")
    
            lower_quark_err = ntrackall4.Get(str(min)+"_SubJet_Central_Quark_"+inputvar + "_err")
            lower_gluon_err = ntrackall4.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar + "_err")
            
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
            higher_mc_unc = unc_array_err(higher_mc,higher_mc_err)
            lower_mc_unc = unc_array_err(lower_mc,lower_mc_err)
            higher_data_unc = unc_data_err(higher_data)
            lower_data_unc = unc_data_err(lower_data)


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
            higher_data_unc = higher_data_unc/higher_data_unc.sum()
            lower_data_unc = lower_data_unc/lower_data_unc.sum()   
       
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
                lower_data_unc=lower_data_unc*factor_quark_unc

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
                lower_mc_unc=lower_mc_unc*factor_gluon_unc
                lower_data_unc=lower_data_unc*factor_gluon_unc

            higher_quark_unc = higher_quark_unc/higher_quark_unc.sum()
            higher_gluon_unc = higher_gluon_unc/higher_gluon_unc.sum()
            lower_quark_unc = lower_quark_unc/lower_quark_unc.sum()
            lower_gluon_unc = lower_gluon_unc/lower_gluon_unc.sum()

            higher_mc_unc = higher_mc_unc/higher_mc_unc.sum()
            lower_mc_unc = lower_mc_unc/lower_mc_unc.sum()   

            higher_data_unc = higher_data_unc/higher_data_unc.sum()
            lower_data_unc = lower_data_unc/lower_data_unc.sum()  
            
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

            F_data_unc = higher_data_unc
            C_data_unc = lower_data_unc
            Q_data_unc = -(C_data_unc*fg_unc-F_data_unc*cg_unc)/(cg_unc*fq_unc-fg_unc*cq_unc)
            G_data_unc = (C_data_unc*fq_unc-F_data_unc*cq_unc)/(cg_unc*fq_unc-fg_unc*cq_unc)
 
            
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
                
            nominal_extract_Q,nominal_extract_G = mc_matrixmethod(lower_quark,lower_gluon,higher_quark,higher_gluon,fg,cg,fq,cq,higher_mc,lower_mc)
            quark_data = higher_data.Clone()
            gluon_data = higher_data.Clone()
            
            
            #do matrix method on data (nominal fraction)
            # first normalize it
            if (higher_data.Integral() != 0):
                higher_data.Scale(1/higher_data.Integral())
            if (lower_data.Integral() != 0):
                lower_data.Scale(1/lower_data.Integral())
              
            extracted_data_nominal_Q,extracted_data_nominal_G = data_matrixmethod(lower_quark,lower_gluon,higher_quark,higher_gluon,higher_data,lower_data,fg,cg,fq,cq)            
            mc_bin,data_bin,sf_q,sf_g,q_eff_mc,g_rej_mc,q_eff_data,g_rej_data = wp_bin(wpoint,nominal_extract_Q,nominal_extract_G,extracted_data_nominal_Q,extracted_data_nominal_G)
            

            Q_eff_data_unc = unp.std_devs(Q_data_unc[0:mc_bin].sum())
            G_rej_data_unc = unp.std_devs(G_data_unc[mc_bin:-1].sum())
            
            Q_eff_unc = unp.std_devs(Q_unc[0:mc_bin].sum())
            G_rej_unc = unp.std_devs(G_unc[mc_bin:-1].sum())



            #calculate the SF uncertainty
            qeff.append(q_eff_mc)
            grej.append(g_rej_mc)

            print(q_eff_data)
            histq.SetBinContent(i+1,q_eff_mc)
            histg.SetBinContent(i+1,g_rej_mc)
            histq.SetBinError(i+1,Q_eff_unc)
            histg.SetBinError(i+1,G_rej_unc) 
            
        histq_array.append(histq)
        histg_array.append(histg)


            

    histq_array[0].SetMaximum(1.2)
    histq_array[0].SetMinimum(0.3)
    histq_array[0].GetXaxis().SetRangeUser(500.,2000.)
    histq_array[0].GetXaxis().SetTitle("Jet p_{T} (GeV)")
    histq_array[0].GetYaxis().SetTitle("Efficiency")
    histq_array[0].SetLabelSize(0.04,"Y")
    histq_array[0].SetTitleOffset(1.1,"Y")
    histq_array[0].SetTitleSize(0.04,"Y")
    histq_array[0].SetMarkerStyle(24)
    histq_array[0].SetLineColor(4)
    histq_array[0].SetMarkerColor(1)
    histq_array[0].SetMarkerSize(1)
    histq_array[0].SetFillColor(4)
    histq_array[0].SetFillStyle(3005)
    
    
    histg_array[0].SetLineColor(2)
    histg_array[0].SetMarkerStyle(32)
    histg_array[0].SetMarkerSize(1)
    histg_array[0].SetMarkerColor(1)
    histg_array[0].SetFillColor(2)
    histg_array[0].SetFillStyle(3005)
    
    histq_array[1].SetMaximum(1.2)
    histq_array[1].SetMinimum(0.3)
    histq_array[1].SetMarkerStyle(24)
    histq_array[1].SetLineColor(5)
    histq_array[1].SetMarkerColor(1)
    histq_array[1].SetMarkerSize(1)
    histq_array[1].SetFillColor(5)
    histq_array[1].SetFillStyle(3005)
    
    
    histg_array[1].SetLineColor(3)
    histg_array[1].SetMarkerStyle(32)
    histg_array[1].SetMarkerSize(1)
    histg_array[1].SetMarkerColor(1)
    histg_array[1].SetFillColor(3)
    histg_array[1].SetFillStyle(3005)
    
    from ROOT import *
    c = TCanvas("","",500,500)
    gStyle.SetOptStat(0)
    leg = TLegend(0.6,0.7,0.84,0.9)
    leg.SetTextFont(42)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    
    leg.AddEntry(histq_array[0],"Quark Efficiency("+mc_sample_1+")","lp")
    leg.AddEntry(histg_array[0],"Gluon Rejection("+mc_sample_1+")","lp")
    leg.AddEntry(histq_array[1],"Quark Efficiency("+mc_sample_2+")","lp")
    leg.AddEntry(histg_array[1],"Gluon Rejection("+mc_sample_2+")","lp")
    
    histq_array[0].Draw("L P0 E2")
    histg_array[0].Draw("L P0 E2 same")
    histq_array[0].Draw("L  same")
    histg_array[0].Draw("L  same")
    histq_array[1].Draw("L P0 E2 same")
    histg_array[1].Draw("L P0 E2 same")
    histq_array[1].Draw("L  same")
    histg_array[1].Draw("L  same")    

    leg.Draw("same")
    myText(0.16,0.84,"#it{#bf{#scale[1.8]{#bf{ATLAS} Internal}}}")
    myText(0.16,0.80,"#bf{#scale[1.5]{#sqrt{s} = 13 TeV ,139 fb^{-1}}}")
    #if var == "ntrk":
    #    myText(0.18,0.76,"#bf{#scale[1.5]{N_{track} Quark Efficiency" + str(100 * wpoint) +"%  }}")
    #if var == "bdt":
    #    myText(0.18,0.76,"#bf{#scale[1.5]{BDT Quark Efficiency" + str(100 * wpoint) +"%  }}")
    
    c.Print(mc_sample_1 +"-"+ mc_sample_2 +"-"+ str(wpoint) +"-"+ str(var)+ "-rej.pdf")
