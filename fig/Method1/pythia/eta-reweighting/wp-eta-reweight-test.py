from ROOT import *
import numpy as np
from prettytable import PrettyTable
import sys
import uncertainties as unc # propagate uncertainties
from uncertainties import unumpy as unp # array operations for type ufloat

    
var = str(sys.argv[2])  #change the var name according to the inputvar you want to read
mc = "pythia"   #by setting it as "SF" or "MC", it will automatically making scale factor plots or MC closure plots

if var == "ntrk":
    doreweight = "Quark"  
if var == "bdt":
    doreweight = "Quark"
    
wpoint_array = [float(sys.argv[1])] 
eta_reweight_array = [0,1,2,3]

    

    
#eta_bin = [0,0.5,1,2.1]


#fdijetMC = TFile("/eos/user/r/rqian/dijet-mono-result/pythia-eta-test/dijet_pythia_etatest.root")
#fdijetData = TFile("/eos/user/r/rqian/dijet-mono-result/data-eta-test/dijet_data_eta_test.root")


    
eta_bin = [0,0.4,0.8,1.4]

bin_width = [0.8,0.8,1.2]

fdijetMC = TFile("/eos/user/r/rqian/dijet-mono-result/pythia-eta-test/dijet_pythia_etatest_1-4.root")

fdijetData = TFile("/eos/user/r/rqian/dijet-mono-result/data-eta-test/dijet_data_1-4.root")

bins = np.array([500.,600.,800.,1000.,1200.,1500.,2000.])
bin = np.array([500,600,800,1000,1200,1500,2000])

def myText(x,y,text, color = 1):
    l = TLatex()
    l.SetTextSize(0.025)
    l.SetNDC()
    l.SetTextColor(color)
    l.DrawLatex(x,y,text)
    pass

def error_check(hist):
    for i in range(1,hist.GetNbinsX()+1):
        if hist.GetBinContent(i)>1:
            hist.SetBinContent(i,0)
    return(hist)

# convert histogram into uncertainty array
def unc_array(hist):
    value = np.zeros(hist.GetNbinsX())
    error = np.zeros(hist.GetNbinsX())
    for j in range(1,hist.GetNbinsX()+1):
        value[j-1] = hist.GetBinContent(j)
        error[j-1] = hist.GetBinError(j)
    result = unp.uarray(value,error)
    return(result)

# convert uncertainty array into histogram
def set_hist(hist,unc):
    for i in range(1,hist.GetNbinsX()+1):
        hist.SetBinError(i,unp.std_devs(unc[i-1]))
        hist.SetBinContent(i,unp.nominal_values(unc[i-1]))


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

# do matrix on data for given normalized distribution
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
    
    
# Get the bin(cut) of given working point for the nominal result and the Quark efficiency , Gluon rejection , scale factor.
def wp_sf_nominal(wpoint,quark_mc,gluon_mc,quark_data,gluon_data):
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

# Get the Quark efficiency , Gluon rejection , scale factor for systematics.
def wp_sf(mc_bin,quark_mc,gluon_mc,quark_data,gluon_data):
    quark_mc_cumsum = np.cumsum(quark_mc)
    gluon_mc_cumsum = np.cumsum(gluon_mc)
    quark_data_cumsum = np.cumsum(quark_data)
    gluon_data_cumsum = np.cumsum(gluon_data)
    #print(quark_mc_cumsum[-1],quark_data_cumsum[-1])

    sf_q = quark_data_cumsum[mc_bin]/quark_mc_cumsum[mc_bin]
    sf_g = (gluon_data_cumsum[-1] - gluon_data_cumsum[mc_bin])/(gluon_mc_cumsum[-1] - gluon_mc_cumsum[mc_bin])
    q_eff_mc = quark_mc_cumsum[mc_bin]/quark_mc_cumsum[-1]
    g_rej_mc = (gluon_mc_cumsum[-1]-gluon_mc_cumsum[mc_bin])/gluon_mc_cumsum[-1]
    q_eff_data = quark_data_cumsum[mc_bin]/quark_data_cumsum[-1]
    g_rej_data = (gluon_data_cumsum[-1]-gluon_data_cumsum[mc_bin])/gluon_data_cumsum[-1]
    #print(gluon_data_cumsum,gluon_mc_cumsum)
    #print(gluon_data_cumsum[data_bin],gluon_mc_cumsum[mc_bin])
    return(sf_q,sf_g,q_eff_mc,g_rej_mc,q_eff_data,g_rej_data)


    
def ReadingDijetDataEta(file0,pt,var,name,index):
    histo_dic = {}
    higher_data = file0.Get(pt+"_LeadingJet_Forward_Data_"+str(index)+"_"+var)
    lower_data = file0.Get(pt+"_LeadingJet_Central_Data_"+str(index)+"_"+var)

    higher_data.SetName(pt+"_LeadingJet_Forward_Data_"+str(index)+var+name)
    lower_data.SetName(pt+"_LeadingJet_Central_Data_"+str(index)+var+name)
    
    higher_data_2 = file0.Get(pt+"_SubJet_Forward_Data_"+str(index)+"_"+var)
    lower_data_2 = file0.Get(pt+"_SubJet_Central_Data_"+str(index)+"_"+var)
    
    higher_data_2.SetName(pt+"_SubJet_Forward_Data_"+str(index)+var+name)
    lower_data_2.SetName(pt+"_SubJet_Central_Data_"+str(index)+var+name)
    
    higher_data.Add(higher_data_2)
    lower_data.Add(lower_data_2)
    

    higher_data_err = file0.Get(pt+"_LeadingJet_Forward_Data_"+str(index)+"_"+var+"_err")
    lower_data_err = file0.Get(pt+"_LeadingJet_Central_Data_"+str(index)+"_"+var+"_err")

    higher_data_err.SetName(pt+"_LeadingJet_Forward_Data_"+var+"_err"+str(index)+var+name)
    lower_data_err.SetName(pt+"_LeadingJet_Central_Data_"+var+"_err"+str(index)+var+name)

    higher_data_2_err = file0.Get(pt+"_SubJet_Forward_Data_"+str(index)+"_"+var+"_err")
    lower_data_2_err = file0.Get(pt+"_SubJet_Central_Data_"+str(index)+"_"+var+"_err")

    higher_data_2_err.SetName(pt+"_SubJet_Forward_Data_"+var+"_err"+str(index)+var+name)
    lower_data_2_err.SetName(pt+"_SubJet_Central_Data_"+var+"_err"+str(index)+var+name)

    higher_data_err.Add(higher_data_2_err)
    lower_data_err.Add(lower_data_2_err)
    
    for i in range(1,higher_data.GetNbinsX()+1):
        higher_data.SetBinError(i,np.sqrt(higher_data_err.GetBinContent(i)))
        lower_data.SetBinError(i,np.sqrt(lower_data_err.GetBinContent(i)))
        
    histo_dic["hd"] =  higher_data
    histo_dic["ld"] =  lower_data

    return(histo_dic)        


def ReadingDijetMCEta(file0,pt,var,name,index):
    histo_dic = {}
    higher_quark = file0.Get(pt+"_LeadingJet_Forward_Quark_"+str(index)+"_"+var)
    lower_quark = file0.Get(pt+"_LeadingJet_Central_Quark_"+str(index)+"_"+var)
    higher_gluon = file0.Get(pt+"_LeadingJet_Forward_Gluon_"+str(index)+"_"+var)
    lower_gluon = file0.Get(pt+"_LeadingJet_Central_Gluon_"+str(index)+"_"+var)
    higher_b_quark = file0.Get(pt+"_LeadingJet_Forward_B_Quark_"+str(index)+"_"+var)
    lower_b_quark = file0.Get(pt+"_LeadingJet_Central_B_Quark_"+str(index)+"_"+var)
    higher_c_quark = file0.Get(pt+"_LeadingJet_Forward_C_Quark_"+str(index)+"_"+var)
    lower_c_quark = file0.Get(pt+"_LeadingJet_Central_C_Quark_"+str(index)+"_"+var)


    higher_quark.SetName(pt+"_LeadingJet_Forward_Quark_"+str(index)+var+name)
    lower_quark.SetName(pt+"_LeadingJet_Central_Quark_"+str(index)+var+name)
    higher_gluon.SetName(pt+"_LeadingJet_Forward_Gluon_"+str(index)+var+name)
    lower_gluon.SetName(pt+"_LeadingJet_Central_Gluon_"+str(index)+var+name)
    higher_b_quark.SetName(pt+"_LeadingJet_Forward_B_Quark_"+str(index)+var+name)
    lower_b_quark.SetName(pt+"_LeadingJet_Central_B_Quark_"+str(index)+var+name)
    higher_c_quark.SetName(pt+"_LeadingJet_Forward_C_Quark_"+str(index)+var+name)
    lower_c_quark.SetName(pt+"_LeadingJet_Central_C_Quark_"+str(index)+var+name)
    

    higher_quark_2 = file0.Get(pt+"_SubJet_Forward_Quark_"+str(index)+"_"+var)
    lower_quark_2 = file0.Get(pt+"_SubJet_Central_Quark_"+str(index)+"_"+var)
    higher_gluon_2 = file0.Get(pt+"_SubJet_Forward_Gluon_"+str(index)+"_"+var)
    lower_gluon_2 = file0.Get(pt+"_SubJet_Central_Gluon_"+str(index)+"_"+var)
    higher_b_quark_2 = file0.Get(pt+"_SubJet_Forward_B_Quark_"+str(index)+"_"+var)
    lower_b_quark_2 = file0.Get(pt+"_SubJet_Central_B_Quark_"+str(index)+"_"+var)
    higher_c_quark_2 = file0.Get(pt+"_SubJet_Forward_C_Quark_"+str(index)+"_"+var)
    lower_c_quark_2 = file0.Get(pt+"_SubJet_Central_C_Quark_"+str(index)+"_"+var)

    higher_quark_2.SetName(pt+"_SubJet_Forward_Quark_"+str(index)+var+name)
    lower_quark_2.SetName(pt+"_SubJet_Central_Quark_"+str(index)+var+name)
    higher_gluon_2.SetName(pt+"_SubJet_Forward_Gluon_"+str(index)+var+name)
    lower_gluon_2.SetName(pt+"_SubJet_Central_Gluon_"+str(index)+var+name)
    higher_b_quark_2.SetName(pt+"_SubJet_Forward_B_Quark_"+str(index)+var+name)
    lower_b_quark_2.SetName(pt+"_SubJet_Central_B_Quark_"+str(index)+var+name)
    higher_c_quark_2.SetName(pt+"_SubJet_Forward_C_Quark_"+str(index)+var+name)
    lower_c_quark_2.SetName(pt+"_SubJet_Central_C_Quark_"+str(index)+var+name)

    higher_quark.Add(higher_quark_2)
    higher_gluon.Add(higher_gluon_2)
    higher_b_quark.Add(higher_b_quark_2)
    higher_c_quark.Add(higher_c_quark_2)

    lower_quark.Add(lower_quark_2)
    lower_gluon.Add(lower_gluon_2)
    lower_b_quark.Add(lower_b_quark_2)
    lower_c_quark.Add(lower_c_quark_2)
 
            
    higher_b_quark.Add(higher_c_quark)
    lower_b_quark.Add(lower_c_quark)

    higher_quark_err = file0.Get(pt+"_LeadingJet_Forward_Quark_"+str(index)+"_"+var+"_err")
    lower_quark_err = file0.Get(pt+"_LeadingJet_Central_Quark_"+str(index)+"_"+var+"_err")
    higher_gluon_err = file0.Get(pt+"_LeadingJet_Forward_Gluon_"+str(index)+"_"+var+"_err")
    lower_gluon_err = file0.Get(pt+"_LeadingJet_Central_Gluon_"+str(index)+"_"+var+"_err")
    higher_b_quark_err = file0.Get(pt+"_LeadingJet_Forward_B_Quark_"+str(index)+"_"+var+"_err")
    lower_b_quark_err = file0.Get(pt+"_LeadingJet_Central_B_Quark_"+str(index)+"_"+var+"_err")
    higher_c_quark_err = file0.Get(pt+"_LeadingJet_Forward_C_Quark_"+str(index)+"_"+var+"_err")
    lower_c_quark_err = file0.Get(pt+"_LeadingJet_Central_C_Quark_"+str(index)+"_"+var+"_err")



    
    higher_quark_err.SetName(pt+"_LeadingJet_Forward_Quark_"+str(index)+"_"+var+"_err"+name)
    lower_quark_err.SetName(pt+"_LeadingJet_Central_Quark_"+str(index)+"_"+var+"_err"+name)
    higher_gluon_err.SetName(pt+"_LeadingJet_Forward_Gluon_"+str(index)+"_"+var+"_err"+name)
    lower_gluon_err.SetName(pt+"_LeadingJet_Central_Gluon_"+str(index)+"_"+var+"_err"+name)
    higher_b_quark_err.SetName(pt+"_LeadingJet_Forward_B_Quark_"+str(index)+"_"+var+"_err"+name)
    lower_b_quark_err.SetName(pt+"_LeadingJet_Central_B_Quark_"+str(index)+"_"+var+"_err"+name)
    higher_c_quark_err.SetName(pt+"_LeadingJet_Forward_C_Quark_"+str(index)+"_"+var+"_err"+name)
    lower_c_quark_err.SetName(pt+"_LeadingJet_Central_C_Quark_"+str(index)+"_"+var+"_err"+name)

    higher_quark_2_err = file0.Get(pt+"_SubJet_Forward_Quark_"+str(index)+"_"+var+"_err")
    lower_quark_2_err = file0.Get(pt+"_SubJet_Central_Quark_"+str(index)+"_"+var+"_err")
    higher_gluon_2_err = file0.Get(pt+"_SubJet_Forward_Gluon_"+str(index)+"_"+var+"_err")
    lower_gluon_2_err = file0.Get(pt+"_SubJet_Central_Gluon_"+str(index)+"_"+var+"_err")
    higher_b_quark_2_err = file0.Get(pt+"_SubJet_Forward_B_Quark_"+str(index)+"_"+var+"_err")
    lower_b_quark_2_err = file0.Get(pt+"_SubJet_Central_B_Quark_"+str(index)+"_"+var+"_err")
    higher_c_quark_2_err = file0.Get(pt+"_SubJet_Forward_C_Quark_"+str(index)+"_"+var+"_err")
    lower_c_quark_2_err = file0.Get(pt+"_SubJet_Central_C_Quark_"+str(index)+"_"+var+"_err")
    
    higher_quark_2_err.SetName(pt+"_SubJet_Forward_Quark_"+str(index)+"_"+var+"_err"+name)
    lower_quark_2_err.SetName(pt+"_SubJet_Central_Quark_"+str(index)+"_"+var+"_err"+name)
    higher_gluon_2_err.SetName(pt+"_SubJet_Forward_Gluon_"+str(index)+"_"+var+"_err"+name)
    lower_gluon_2_err.SetName(pt+"_SubJet_Central_Gluon_"+str(index)+"_"+var+"_err"+name)
    higher_b_quark_2_err.SetName(pt+"_SubJet_Forward_B_Quark_"+str(index)+"_"+var+"_err"+name)
    lower_b_quark_2_err.SetName(pt+"_SubJet_Central_B_Quark_"+str(index)+"_"+var+"_err"+name)
    higher_c_quark_2_err.SetName(pt+"_SubJet_Forward_C_Quark_"+str(index)+"_"+var+"_err"+name)
    lower_c_quark_2_err.SetName(pt+"_SubJet_Central_C_Quark_"+str(index)+"_"+var+"_err"+name)

    higher_quark_err.Add(higher_quark_2_err)
    higher_gluon_err.Add(higher_gluon_2_err)
    higher_b_quark_err.Add(higher_b_quark_2_err)
    higher_c_quark_err.Add(higher_c_quark_2_err)

    lower_quark_err.Add(lower_quark_2_err)
    lower_gluon_err.Add(lower_gluon_2_err)
    lower_b_quark_err.Add(lower_b_quark_2_err)
    lower_c_quark_err.Add(lower_c_quark_2_err)
 
            
    higher_b_quark_err.Add(higher_c_quark_err)
    lower_b_quark_err.Add(lower_c_quark_err)
    
    for i in range(1,higher_quark.GetNbinsX()+1):
        higher_quark.SetBinError(i,np.sqrt(higher_quark_err.GetBinContent(i)))
        higher_gluon.SetBinError(i,np.sqrt(higher_gluon_err.GetBinContent(i)))
        higher_b_quark.SetBinError(i,np.sqrt(higher_b_quark_err.GetBinContent(i)))
        lower_quark.SetBinError(i,np.sqrt(lower_quark_err.GetBinContent(i)))
        lower_gluon.SetBinError(i,np.sqrt(lower_gluon_err.GetBinContent(i)))
        lower_b_quark.SetBinError(i,np.sqrt(lower_b_quark_err.GetBinContent(i)))

    
    histo_dic["hq"] =  higher_quark
    histo_dic["hg"] =  higher_gluon
    histo_dic["ho"] =  higher_b_quark
    
    histo_dic["hmc"] =  higher_quark.Clone()
    histo_dic["hmc"].SetName(pt+"higher_all_jets"+var+"_"+name)
    histo_dic["hmc"].Add(higher_gluon)
    histo_dic["hmc"].Add(higher_b_quark)
    
    histo_dic["lq"] =  lower_quark
    histo_dic["lg"] =  lower_gluon
    histo_dic["lo"] =  lower_b_quark
    
    histo_dic["lmc"] =  lower_quark.Clone()
    histo_dic["lmc"].SetName(pt+"lower_all_jets"+var+"_"+name)
    histo_dic["lmc"].Add(lower_gluon)
    histo_dic["lmc"].Add(lower_b_quark)
    
    return(histo_dic)

# Get the bin(cut) of given working point for the nominal result and the Quark efficiency , Gluon rejection , scale factor.
def wp_sf_nominal(wpoint,quark_mc,gluon_mc,quark_data,gluon_data):
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

# Get the Quark efficiency , Gluon rejection , scale factor for systematics.
def wp_sf(mc_bin,quark_mc,gluon_mc,quark_data,gluon_data):
    quark_mc_cumsum = np.cumsum(quark_mc)
    gluon_mc_cumsum = np.cumsum(gluon_mc)
    quark_data_cumsum = np.cumsum(quark_data)
    gluon_data_cumsum = np.cumsum(gluon_data)
    #print(quark_mc_cumsum[-1],quark_data_cumsum[-1])

    sf_q = quark_data_cumsum[mc_bin]/quark_mc_cumsum[mc_bin]
    sf_g = (gluon_data_cumsum[-1] - gluon_data_cumsum[mc_bin])/(gluon_mc_cumsum[-1] - gluon_mc_cumsum[mc_bin])
    q_eff_mc = quark_mc_cumsum[mc_bin]/quark_mc_cumsum[-1]
    g_rej_mc = (gluon_mc_cumsum[-1]-gluon_mc_cumsum[mc_bin])/gluon_mc_cumsum[-1]
    q_eff_data = quark_data_cumsum[mc_bin]/quark_data_cumsum[-1]
    g_rej_data = (gluon_data_cumsum[-1]-gluon_data_cumsum[mc_bin])/gluon_data_cumsum[-1]
    #print(gluon_data_cumsum,gluon_mc_cumsum)
    #print(gluon_data_cumsum[data_bin],gluon_mc_cumsum[mc_bin])
    return(sf_q,sf_g,q_eff_mc,g_rej_mc,q_eff_data,g_rej_data)


def bootstrap_result(data):
    n1 = data.Clone()
    result = data.Clone()
    for j in range(1,higher_quark.GetNbinsX()+1):
            result.SetBinContent(j,np.random.poisson(n1.GetBinContent(j)))
    return(result)  


qeff = []
grej = []
qsf_array = []
gsf_array = []
unc_array_q = []
unc_array_g = []

for wpoint in wpoint_array:

    mc_bin_array = []
    for eta_reweight in eta_reweight_array:
        unc_form_q = []  
        unc_form_g = []
        from ROOT import *
        histq = TH1F("histq"+str(eta_reweight),"",6,bins)
        histg = TH1F("histg"+str(eta_reweight),"",6,bins)
        histqsf = TH1F("histqsf"+str(eta_reweight),"",6,bins)
        histgsf = TH1F("histgsf"+str(eta_reweight),"",6,bins)
        
        if eta_reweight == 0:
            nevents_eta = []
            for i in range(6):
                nevents_eta.append(np.ones(3))
    

        if eta_reweight == 1:
            nevents_eta = np.load("./mc_eta_events.npy")
            
        if eta_reweight == 2:
            nevents_eta = []
            for i in range(6):
                nevents_eta.append([1,2,3]) 
            
        if eta_reweight == 3:
            nevents_eta = []
            for i in range(6):
                nevents_eta.append([1,10,20]) 
        
        
        for i in range(0,6):  #for only dijet event, start frTom jet pT>500 GeV
                    unc_form_q_pt = []
                    unc_form_g_pt = []
        
                    c = TCanvas("c","c",500,500)
                    min = bin[i]
                    max = bin[i+1]
                    print(min)
                    
                    nominal_jets =  ReadingDijetMCEta(fdijetMC,str(min),var,"0"+str(eta_reweight),"0")
                    
                    data_jets = ReadingDijetDataEta(fdijetData,str(min),var,"0"+str(eta_reweight),"0")
        
                            
                    higher_quark = nominal_jets["hq"]
                    lower_quark = nominal_jets["lq"]
                    higher_gluon = nominal_jets["hg"]
                    lower_gluon = nominal_jets["lg"]
                    higher_other = nominal_jets["ho"]
                    lower_other = nominal_jets["lo"]        
                    higher_mc = nominal_jets["hmc"]
                    lower_mc = nominal_jets["lmc"]        
        
                    
                    higher_data_TH1I = data_jets["hd"]
                    higher_data = higher_quark.Clone(higher_data_TH1I.GetName()+"TH1F"+str(eta_reweight))
                    for j in range(1,higher_data.GetNbinsX()+1):
                            higher_data.SetBinContent(j,float(higher_data_TH1I.GetBinContent(j)))
                            higher_data.SetBinError(j,float(higher_data_TH1I.GetBinError(j)))
                            
                            
                    lower_data_TH1I = data_jets["ld"]
                    lower_data = lower_quark.Clone(lower_data_TH1I.GetName()+"TH1F"+str(eta_reweight))
                    for j in range(1,higher_data.GetNbinsX()+1):
                            lower_data.SetBinContent(j,float(lower_data_TH1I.GetBinContent(j)))
                            lower_data.SetBinError(j,float(lower_data_TH1I.GetBinError(j))) 
                            
                    higher_quark.Scale(nevents_eta[i][0])
                    higher_gluon.Scale(nevents_eta[i][0])
                    higher_other.Scale(nevents_eta[i][0])
                    higher_mc.Scale(nevents_eta[i][0])
                    lower_quark.Scale(nevents_eta[i][0])
                    lower_gluon.Scale(nevents_eta[i][0])
                    lower_other.Scale(nevents_eta[i][0])
                    lower_mc.Scale(nevents_eta[i][0])
                    higher_data.Scale(nevents_eta[i][0])
                    lower_data.Scale(nevents_eta[i][0])
                    
                    for e in range(1,len(eta_bin)-1):
        
                            print(var,str(eta_bin[e]))
                            
                            nominal_jets_add =  ReadingDijetMCEta(fdijetMC,str(min),var,str(eta_bin[e])+str(eta_reweight),str(eta_bin[e]))
                            
        
                            data_jets_add = ReadingDijetDataEta(fdijetData,str(min),var,str(eta_bin[e])+str(eta_reweight),str(eta_bin[e]))
                            
                            higher_data_TH1I_add = data_jets_add["hd"]
                            higher_data_add = higher_quark.Clone(higher_data_TH1I_add.GetName()+"TH1F"+str(eta_reweight))
                            for j in range(1,higher_data.GetNbinsX()+1):
                                    higher_data_add.SetBinContent(j,float(higher_data_TH1I_add.GetBinContent(j)))
                                    higher_data_add.SetBinError(j,float(higher_data_TH1I_add.GetBinError(j)))
                                    
                                    
                            lower_data_TH1I_add = data_jets_add["ld"]
                            lower_data_add = lower_quark.Clone(lower_data_TH1I_add.GetName()+"TH1F"+str(eta_reweight))
                            for j in range(1,higher_data.GetNbinsX()+1):
                                    lower_data_add.SetBinContent(j,float(lower_data_TH1I_add.GetBinContent(j)))
                                    lower_data_add.SetBinError(j,float(lower_data_TH1I_add.GetBinError(j)))   
                            
                            # flatten the eta distribution
                            nominal_jets_add["hq"].Scale(nevents_eta[i][e])
                            nominal_jets_add["hg"].Scale(nevents_eta[i][e])
                            nominal_jets_add["ho"].Scale(nevents_eta[i][e])
                            nominal_jets_add["hmc"].Scale(nevents_eta[i][e])
                            nominal_jets_add["lq"].Scale(nevents_eta[i][e])
                            nominal_jets_add["lg"].Scale(nevents_eta[i][e])
                            nominal_jets_add["lo"].Scale(nevents_eta[i][e])
                            nominal_jets_add["lmc"].Scale(nevents_eta[i][e])
                            higher_data_add.Scale(nevents_eta[i][e])
                            lower_data_add.Scale(nevents_eta[i][e])
                            
                            higher_quark.Add(nominal_jets_add["hq"])
                            lower_quark.Add(nominal_jets_add["lq"])
                            higher_gluon.Add(nominal_jets_add["hg"])
                            lower_gluon.Add(nominal_jets_add["lg"])
                            higher_other.Add(nominal_jets_add["ho"])
                            lower_other.Add(nominal_jets_add["lo"])        
                            higher_mc.Add(nominal_jets_add["hmc"])
                            lower_mc.Add(nominal_jets_add["lmc"])  
                            
                            higher_data.Add(higher_data_add)
                            lower_data.Add(lower_data_add) 
                    
                    higher_quark_unc = unc_array(higher_quark)
                    higher_gluon_unc = unc_array(higher_gluon)
                    lower_quark_unc = unc_array(lower_quark)
                    lower_gluon_unc = unc_array(lower_gluon)
                    higher_mc_unc = unc_array(higher_mc)
                    lower_mc_unc = unc_array(lower_mc)
                    
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
                        lower_mc_unc=lower_mc_unc*factor_gluon_unc
                    
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
                    F_unc = higher_unc
                    C_unc = lower_unc
                    Q_unc = -(C_unc*fg_unc-F_unc*cg_unc)/(cg_unc*fq_unc-fg_unc*cq_unc)
                    G_unc = (C_unc*fq_unc-F_unc*cq_unc)/(cg_unc*fq_unc-fg_unc*cq_unc)
                    
        
                    
                    #compute the uncertainty and plots
        
        
                    higher_data_nominal = higher_data.Clone()
                    lower_data_nominal = lower_data.Clone()
                    
                    
                    #do matrix method on data (nominal fraction)
                    # first normalize it
                    if (higher_data_nominal.Integral() != 0):
                        higher_data_nominal.Scale(1/higher_data_nominal.Integral())
                    if (lower_data_nominal.Integral() != 0):
                        lower_data_nominal.Scale(1/lower_data_nominal.Integral())
                    # do matrix method to extract the distribution of sherpa first
                    fg,cg,fq,cq = fraction(lower_quark,lower_gluon,higher_quark,higher_gluon)
                    
                    
                    # normalize the nominal mc
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
                    
                    
                    #Get nominal distribution , cut bin , efficiency/rejection and scale factors
                    nominal_extract_Q,nominal_extract_G = mc_matrixmethod(lower_quark,lower_gluon,higher_quark,higher_gluon,fg,cg,fq,cq,higher_mc,lower_mc)
                    extracted_data_nominal_Q,extracted_data_nominal_G = data_matrixmethod(lower_quark,lower_gluon,higher_quark,higher_gluon,higher_data_nominal,lower_data_nominal,fg,cg,fq,cq)
                    
                    if eta_reweight == 0:
                        mc_bin,data_bin,sf_q,sf_g,q_eff_mc,g_rej_mc,q_eff_data,g_rej_data=wp_sf_nominal(wpoint,nominal_extract_Q,nominal_extract_G,extracted_data_nominal_Q,extracted_data_nominal_G)
                        mc_bin_array.append(mc_bin)
                    else :
                        mc_bin = mc_bin_array[i]
                        sf_q,sf_g,q_eff_mc,g_rej_mc,q_eff_data,g_rej_data = wp_sf(mc_bin,nominal_extract_Q,nominal_extract_G,extracted_data_nominal_Q,extracted_data_nominal_G)
                    
                    # reweighting systematics calculation
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
        
        
                    higher_quark_bootstrap = higher_quark.Clone()
                    higher_gluon_bootstrap = higher_gluon.Clone()
                    lower_quark_bootstrap = lower_quark.Clone()
                    lower_gluon_bootstrap = lower_gluon.Clone()
                    higher_data_bootstrap = higher_data.Clone()
                    lower_data_bootstrap = lower_data.Clone() 
                    
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
        
                            for j in range(1,higher_quark.GetNbinsX()+1):
                                    F_data = forward_data_strap.GetBinContent(j)
                                    C_data = central_data_strap.GetBinContent(j)
                                    Q_data = -(C_data*fg-F_data*cg)/(cg*fq-fg*cq)
                                    G_data = (C_data*fq-F_data*cq)/(cg*fq-fg*cq)
                                    bootstrap_extracted_Q.SetBinContent(j,Q_data)
                                    bootstrap_extracted_G.SetBinContent(j,G_data)
                            
                            sf_q_bootstrap,sf_g_bootstrap,q_eff_mc_bootstrap,g_rej_mc_bootstrap,q_eff_data_bootstrap,g_rej_data_bootstrap = wp_sf(mc_bin,nominal_extract_Q,nominal_extract_G,bootstrap_extracted_Q,bootstrap_extracted_G)
        
                            
                            SF_Qvals.append(sf_q_bootstrap)
                            SF_Gvals.append(sf_g_bootstrap)
        
                    #compute the uncertainty and plots
        
                    
                    SF_Qvals.sort()
                    SF_Gvals.sort()
        
                    
                    sigmaQ = .5*(SF_Qvals[int(.84*len(SF_Qvals))] - SF_Qvals[int(.16*len(SF_Qvals))])
                    sigmaG = .5*(SF_Gvals[int(.84*len(SF_Gvals))] - SF_Gvals[int(.16*len(SF_Gvals))])
        
        
                    #print("statistical: q = "+str(sigmaQ)+" | g = "+str(sigmaG))                  
                    
                    sigmaQ = np.abs(sigmaQ)
                    sigmaG = np.abs(sigmaG)
                    
                    
                    ##Reweight Factor Uncertainty
                    # first normalize it
                    if (higher_data.Integral() != 0):
                        higher_data.Scale(1/higher_data.Integral())
                    if (lower_data.Integral() != 0):
                        lower_data.Scale(1/lower_data.Integral())
                      
                    doreweight_normal = doreweight
                    doreweight = "Quark"
                    extracted_data_Quark_Factor_Q,extracted_data_Quark_Factor_G = data_matrixmethod(lower_quark,lower_gluon,higher_quark,higher_gluon,higher_data,lower_data,fg,cg,fq,cq)            
                    doreweight = "Gluon"
                    extracted_data_Gluon_Factor_Q,extracted_data_Gluon_Factor_G = data_matrixmethod(lower_quark,lower_gluon,higher_quark,higher_gluon,higher_data,lower_data,fg,cg,fq,cq)            
                    doreweight = doreweight_normal
                    
        
                    
                    Q_eff_unc = unp.std_devs(Q_unc[0:mc_bin].sum())
                    G_rej_unc = unp.std_devs(G_unc[mc_bin:-1].sum())
        
                    ##Reweight Factor Uncertainty
                    sf_q_QF,sf_g_QF,q_eff_mc_QF,g_rej_mc_QF,q_eff_data_QF,g_rej_data_QF = wp_sf(mc_bin,Quark_Factor_extract_Q,Quark_Factor_extract_G,extracted_data_Quark_Factor_Q,extracted_data_Quark_Factor_G)
                    sf_q_GF,sf_g_GF,q_eff_mc_GF,g_rej_mc_GF,q_eff_data_GF,g_rej_data_GF = wp_sf(mc_bin,Gluon_Factor_extract_Q,Gluon_Factor_extract_G,extracted_data_Gluon_Factor_Q,extracted_data_Gluon_Factor_G)
                    RF_unc_Q = abs(sf_q_QF-sf_q_GF)
                    RF_unc_G = abs(sf_g_QF-sf_g_GF)
                                
                    #mc closure
                    sf_q_closure,sf_g_closure,q_eff_closure,g_rej_closure,q_eff_data_closure,g_rej_data_closure = wp_sf(mc_bin,higher_quark,higher_gluon,extracted_data_nominal_Q,extracted_data_nominal_G)
                    closure_unc_Q = abs(sf_q_closure-sf_q)
                    closure_unc_G = abs(sf_g_closure-sf_g)            
                    
                    #calculate the SF uncertainty
                    qeff.append(q_eff_mc)
                    grej.append(g_rej_mc)
                    qsf_array.append(sf_q)
                    gsf_array.append(sf_g)
        
        
                    unc_q = np.sqrt(closure_unc_Q**2+sigmaQ**2+RF_unc_Q**2)
                    unc_g = np.sqrt(closure_unc_G**2+sigmaG**2+RF_unc_G**2)
        
        
                    pt_range = str(str(min)+"-"+str(max))
            
                  
            
                    unc_form_q_pt = [pt_range,round(sf_q,5),round(unc_q,5),round(sigmaQ,5),round(closure_unc_Q,5),round(RF_unc_Q,5)]
                    unc_form_g_pt = [pt_range,round(sf_g,5),round(unc_g,5),round(sigmaG,5),round(closure_unc_G,5),round(RF_unc_G,5)]

                    unc_form_q.append(unc_form_q_pt)
                    unc_form_g.append(unc_form_g_pt)

                    histq.SetBinContent(i+1,q_eff_mc)
                    histg.SetBinContent(i+1,g_rej_mc)
                    histq.SetBinError(i+1,Q_eff_unc)
                    histg.SetBinError(i+1,G_rej_unc) 
                    histqsf.SetBinContent(i+1,sf_q)
                    histgsf.SetBinContent(i+1,sf_g)
                    histqsf.SetBinError(i+1,unc_q)
                    histgsf.SetBinError(i+1,unc_g)    
    



     
    
        histq.SetMaximum(1.2)
        histq.SetMinimum(0.3)
        histq.GetXaxis().SetRangeUser(500.,2000.)
        histq.GetXaxis().SetTitle("Jet p_{T} (GeV)")
        histq.GetYaxis().SetTitle("Quark Efficiency,Gluon Rejection")
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
                
        histqsf.SetMaximum(1.6)
        histqsf.SetMinimum(0.4)
        histqsf.GetXaxis().SetRangeUser(500.,2000.)
        histqsf.GetXaxis().SetTitle("Jet p_{T} (GeV)")
        histqsf.GetYaxis().SetTitle("Scale Factor")
        histqsf.SetLabelSize(0.04,"Y")
        histqsf.SetTitleOffset(1.1,"Y")
        histqsf.SetTitleSize(0.04,"Y")
                
        histqsf.SetLineColor(4)
        histqsf.SetMarkerStyle(25)
        histqsf.SetMarkerSize(1)
        histqsf.SetMarkerColor(1)
        histqsf.SetLineStyle(7)
        histqsf.SetFillColor(4)
        histqsf.SetFillStyle(3005)
                
        histgsf.SetLineColor(2)
        histgsf.SetMarkerStyle(26)
        histgsf.SetMarkerSize(1)
        histgsf.SetMarkerColor(1)
        histgsf.SetLineStyle(7)
        histgsf.SetFillColor(2)
        histgsf.SetFillStyle(3005)
                
        from ROOT import *
        c = TCanvas("","",500,500)
        gStyle.SetOptStat(0)
        leg = TLegend(0.6,0.7,0.84,0.9)
        leg.SetTextFont(42)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.AddEntry(histq,"Quark Efficiency","lp")
        leg.AddEntry(histg,"Gluon Rejection","lp")
        histq.Draw("L P0 E2")
        histg.Draw("L P0 E2 same")
        histq.Draw("L  same")
        histg.Draw("L  same")
        leg.Draw("same")
        myText(0.16,0.84,"#it{#bf{#scale[1.8]{#bf{ATLAS} Internal}}}")
        myText(0.16,0.80,"#bf{#scale[1.5]{#sqrt{s} = 13 TeV ,139 fb^{-1}}}")
        #if var == "ntrk":
        #    myText(0.18,0.76,"#bf{#scale[1.5]{N_{track} Quark Efficiency" + str(100 * wpoint) +"%  }}")
        #if var == "bdt":
        #    myText(0.18,0.76,"#bf{#scale[1.5]{BDT Quark Efficiency" + str(100 * wpoint) +"%  }}")
        
        c.Print("dijet-wp-" + str(wpoint) +"-"+ str(var)+ "_"+str(eta_reweight)+".pdf")
                
                
        c1 = TCanvas("","",500,500)
        
        leg2 = TLegend(0.6,0.7,0.84,0.9)
        leg2.SetTextFont(42)
        leg2.SetFillColor(0)
        leg2.SetBorderSize(0)
        leg2.SetFillStyle(0)
        leg2.AddEntry(histqsf,"Quark Scale Factor","lp")
        leg2.AddEntry(histgsf,"Gluon Scale Factor","lp")
        histqsf.Draw("L P0 E2")
        histgsf.Draw("L P0 E2 same")
        histqsf.Draw("L  same")
        histgsf.Draw("L  same")
        leg2.Draw("same")
        myText(0.16,0.84,"#it{#bf{#scale[1.8]{#bf{ATLAS} Internal}}}")
        myText(0.16,0.80,"#bf{#scale[1.5]{#sqrt{s} = 13 TeV ,139 fb^{-1}}}")
        #if var == "ntrk":
        #    myText(0.18,0.76,"#bf{#scale[1.5]{N_{track} Quark Efficiency" + str(100 * wpoint) +"%  }}")
        #if var == "bdt":
        #    myText(0.18,0.76,"#bf{#scale[1.5]{BDT Quark Efficiency" + str(100 * wpoint) +"%  }}")
        
        c1.Print("dijet-sf-" + str(wpoint) +"-"+ str(var)+ "_"+str(eta_reweight)+".pdf")
        tbq = PrettyTable()
        #tbq.title=(str(wpoint*100)+" quark efficiency,"+"eta reweighting "+str(eta_reweight)+" " +var)
        print(str(wpoint*100)+" quark efficiency,"+"eta reweighting "+str(eta_reweight)+" " +var)

        tbq.field_names = ["Pt(GeV)","Quark Scale Factor","Uncertainty","Statistical","MC Closure","Reweight Factor"]
        
        for i in range(6):
            tbq.add_row(unc_form_q[i])
        print(tbq)
        
        tbg = PrettyTable()
        #tbq.title=(str(wpoint*100)+" quark efficiency,"+"eta reweighting "+str(eta_reweight)+" " +var)
        print(str(wpoint*100)+" quark efficiency,"+"eta reweighting "+str(eta_reweight)+" " +var)
        tbg.field_names = ["Pt(GeV)","Gluon Scale Factor","Uncertainty","Statistical","MC Closure","Reweight Factor"]
        
        for i in range(6):
            tbg.add_row(unc_form_g[i])
        print(tbg)