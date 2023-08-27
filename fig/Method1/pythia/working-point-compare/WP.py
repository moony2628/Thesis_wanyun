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
        
def ReadingDijetMCindex(file0,pt,var,name,index):
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


def ReadingDijetMC(file0,pt,var,name):
    histo_dic = {}
    higher_quark = file0.Get(pt+"_LeadingJet_Forward_Quark_"+var)
    lower_quark = file0.Get(pt+"_LeadingJet_Central_Quark_"+var)
    higher_gluon = file0.Get(pt+"_LeadingJet_Forward_Gluon_"+var)
    lower_gluon = file0.Get(pt+"_LeadingJet_Central_Gluon_"+var)
    higher_b_quark = file0.Get(pt+"_LeadingJet_Forward_B_Quark_"+var)
    lower_b_quark = file0.Get(pt+"_LeadingJet_Central_B_Quark_"+var)
    higher_c_quark = file0.Get(pt+"_LeadingJet_Forward_C_Quark_"+var)
    lower_c_quark = file0.Get(pt+"_LeadingJet_Central_C_Quark_"+var)

    
    higher_quark.SetName(pt+"_LeadingJet_Forward_Quark_"+var+name)
    lower_quark.SetName(pt+"_LeadingJet_Central_Quark_"+var+name)
    higher_gluon.SetName(pt+"_LeadingJet_Forward_Gluon_"+var+name)
    lower_gluon.SetName(pt+"_LeadingJet_Central_Gluon_"+var+name)
    higher_b_quark.SetName(pt+"_LeadingJet_Forward_B_Quark_"+var+name)
    lower_b_quark.SetName(pt+"_LeadingJet_Central_B_Quark_"+var+name)
    higher_c_quark.SetName(pt+"_LeadingJet_Forward_C_Quark_"+var+name)
    lower_c_quark.SetName(pt+"_LeadingJet_Central_C_Quark_"+var+name)
    

    higher_quark_2 = file0.Get(pt+"_SubJet_Forward_Quark_"+var)
    lower_quark_2 = file0.Get(pt+"_SubJet_Central_Quark_"+var)
    higher_gluon_2 = file0.Get(pt+"_SubJet_Forward_Gluon_"+var)
    lower_gluon_2 = file0.Get(pt+"_SubJet_Central_Gluon_"+var)
    higher_b_quark_2 = file0.Get(pt+"_SubJet_Forward_B_Quark_"+var)
    lower_b_quark_2 = file0.Get(pt+"_SubJet_Central_B_Quark_"+var)
    higher_c_quark_2 = file0.Get(pt+"_SubJet_Forward_C_Quark_"+var)
    lower_c_quark_2 = file0.Get(pt+"_SubJet_Central_C_Quark_"+var)


    higher_quark_2.SetName(pt+"_SubJet_Forward_Quark_"+var+name)
    lower_quark_2.SetName(pt+"_SubJet_Central_Quark_"+var+name)
    higher_gluon_2.SetName(pt+"_SubJet_Forward_Gluon_"+var+name)
    lower_gluon_2.SetName(pt+"_SubJet_Central_Gluon_"+var+name)
    higher_b_quark_2.SetName(pt+"_SubJet_Forward_B_Quark_"+var+name)
    lower_b_quark_2.SetName(pt+"_SubJet_Central_B_Quark_"+var+name)
    higher_c_quark_2.SetName(pt+"_SubJet_Forward_C_Quark_"+var+name)
    lower_c_quark_2.SetName(pt+"_SubJet_Central_C_Quark_"+var+name)

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

    higher_quark_err = file0.Get(pt+"_LeadingJet_Forward_Quark_"+var+"_err")
    lower_quark_err = file0.Get(pt+"_LeadingJet_Central_Quark_"+var+"_err")
    higher_gluon_err = file0.Get(pt+"_LeadingJet_Forward_Gluon_"+var+"_err")
    lower_gluon_err = file0.Get(pt+"_LeadingJet_Central_Gluon_"+var+"_err")
    higher_b_quark_err = file0.Get(pt+"_LeadingJet_Forward_B_Quark_"+var+"_err")
    lower_b_quark_err = file0.Get(pt+"_LeadingJet_Central_B_Quark_"+var+"_err")
    higher_c_quark_err = file0.Get(pt+"_LeadingJet_Forward_C_Quark_"+var+"_err")
    lower_c_quark_err = file0.Get(pt+"_LeadingJet_Central_C_Quark_"+var+"_err")

    
    higher_quark_err.SetName(pt+"_LeadingJet_Forward_Quark_"+var+"_err"+name)
    lower_quark_err.SetName(pt+"_LeadingJet_Central_Quark_"+var+"_err"+name)
    higher_gluon_err.SetName(pt+"_LeadingJet_Forward_Gluon_"+var+"_err"+name)
    lower_gluon_err.SetName(pt+"_LeadingJet_Central_Gluon_"+var+"_err"+name)
    higher_b_quark_err.SetName(pt+"_LeadingJet_Forward_B_Quark_"+var+"_err"+name)
    lower_b_quark_err.SetName(pt+"_LeadingJet_Central_B_Quark_"+var+"_err"+name)
    higher_c_quark_err.SetName(pt+"_LeadingJet_Forward_C_Quark_"+var+"_err"+name)
    lower_c_quark_err.SetName(pt+"_LeadingJet_Central_C_Quark_"+var+"_err"+name)
    
    higher_quark_2_err = file0.Get(pt+"_SubJet_Forward_Quark_"+var+"_err")
    lower_quark_2_err = file0.Get(pt+"_SubJet_Central_Quark_"+var+"_err")
    higher_gluon_2_err = file0.Get(pt+"_SubJet_Forward_Gluon_"+var+"_err")
    lower_gluon_2_err = file0.Get(pt+"_SubJet_Central_Gluon_"+var+"_err")
    higher_b_quark_2_err = file0.Get(pt+"_SubJet_Forward_B_Quark_"+var+"_err")
    lower_b_quark_2_err = file0.Get(pt+"_SubJet_Central_B_Quark_"+var+"_err")
    higher_c_quark_2_err = file0.Get(pt+"_SubJet_Forward_C_Quark_"+var+"_err")
    lower_c_quark_2_err = file0.Get(pt+"_SubJet_Central_C_Quark_"+var+"_err")


    higher_quark_2_err.SetName(pt+"_SubJet_Forward_Quark_"+var+"_err"+name)
    lower_quark_2_err.SetName(pt+"_SubJet_Central_Quark_"+var+"_err"+name)
    higher_gluon_2_err.SetName(pt+"_SubJet_Forward_Gluon_"+var+"_err"+name)
    lower_gluon_2_err.SetName(pt+"_SubJet_Central_Gluon_"+var+"_err"+name)
    higher_b_quark_2_err.SetName(pt+"_SubJet_Forward_B_Quark_"+var+"_err"+name)
    lower_b_quark_2_err.SetName(pt+"_SubJet_Central_B_Quark_"+var+"_err"+name)
    higher_c_quark_2_err.SetName(pt+"_SubJet_Forward_C_Quark_"+var+"_err"+name)
    lower_c_quark_2_err.SetName(pt+"_SubJet_Central_C_Quark_"+var+"_err"+name)

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

def ReadingDijetData(file0,pt,var,name):
    histo_dic = {}
    higher_data = file0.Get(pt+"_LeadingJet_Forward_Data_"+var)
    lower_data = file0.Get(pt+"_LeadingJet_Central_Data_"+var)

    higher_data.SetName(pt+"_LeadingJet_Forward_Data_"+var+name)
    lower_data.SetName(pt+"_LeadingJet_Central_Data_"+var+name)
    
    higher_data_2 = file0.Get(pt+"_SubJet_Forward_Data_"+var)
    lower_data_2 = file0.Get(pt+"_SubJet_Central_Data_"+var)
    
    higher_data_2.SetName(pt+"_SubJet_Forward_Data_"+var+name)
    lower_data_2.SetName(pt+"_SubJet_Central_Data_"+var+name)
    
    higher_data.Add(higher_data_2)
    lower_data.Add(lower_data_2)
    

    higher_data_err = file0.Get(pt+"_LeadingJet_Forward_Data_"+var+"_err")
    lower_data_err = file0.Get(pt+"_LeadingJet_Central_Data_"+var+"_err")

    higher_data_err.SetName(pt+"_LeadingJet_Forward_Data_"+var+"_err"+name)
    lower_data_err.SetName(pt+"_LeadingJet_Central_Data_"+var+"_err"+name)

    higher_data_2_err = file0.Get(pt+"_SubJet_Forward_Data_"+var+"_err")
    lower_data_2_err = file0.Get(pt+"_SubJet_Central_Data_"+var+"_err")

    higher_data_2_err.SetName(pt+"_SubJet_Forward_Data_"+var+"_err"+name)
    lower_data_2_err.SetName(pt+"_SubJet_Central_Data_"+var+"_err"+name)

    higher_data_err.Add(higher_data_2_err)
    lower_data_err.Add(lower_data_2_err)
    
    for i in range(1,higher_data.GetNbinsX()+1):
        higher_data.SetBinError(i,np.sqrt(higher_data_err.GetBinContent(i)))
        lower_data.SetBinError(i,np.sqrt(lower_data_err.GetBinContent(i)))
        
    histo_dic["hd"] =  higher_data
    histo_dic["ld"] =  lower_data

    return(histo_dic)        

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

## do matrix on mc for given normalized distribution
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

#def bootstrap_result(data):
#    w = 58.45/39.91 # prescale factor 
#    n1 = data.Clone()
#    n2 = data.Clone()
#    result = data.Clone()
#    for j in range(1,higher_quark.GetNbinsX()+1):
#            n2.SetBinContent(j,int((data.GetBinContent(j)-data.GetBinError(j)**2)/(w-w**2)))#(events num for prescaled data)
#            n1.SetBinContent(j,int((data.GetBinContent(j)-w*n2.GetBinContent(j))))
#            result.SetBinContent(j,w*np.random.poisson(n2.GetBinContent(j))+np.random.poisson(n1.GetBinContent(j)))
#    return(result)  

def bootstrap_result(data):
    n1 = data.Clone()
    result = data.Clone()
    for j in range(1,higher_quark.GetNbinsX()+1):
            result.SetBinContent(j,np.random.poisson(n1.GetBinContent(j)))
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

                
ntrackall = TFile("/eos/user/r/rqian/dijet-mono-result/sherpa/dijet_sherpa.root")
ntrackall3 = TFile("/eos/user/r/rqian/dijet-mono-result/data/dijet_data.root")
ntrackall4 = TFile("/eos/user/r/rqian/dijet-mono-result/pythia/dijet_pythia_mono.root")
ntrackall5 = TFile("/eos/user/r/rqian/dijet-mono-result/pdf/dijet_pythia_pdf.root")
fherdipo = TFile("/eos/user/r/rqian/dijet-mono-result/herdipo/dijet_herdipo.root")
fherang = TFile("/eos/user/r/rqian/dijet-mono-result/herang/dijet_herang.root")
sherpa_lund = TFile("/eos/user/r/rqian/dijet-mono-result/sherpa_lund/dijet_sherpa_lund.root")
sv = TFile("/eos/user/r/rqian/dijet-mono-result/sv/dijet_pythia_sv.root")
powpyt = TFile("/eos/user/r/rqian/dijet-mono-result/powpyt/dijet_powpyt.root")

histq = TH1F("histq","",6,bins)
histg = TH1F("histg","",6,bins)
histqsf = TH1F("histqsf","",6,bins)
histgsf = TH1F("histgsf","",6,bins)


unc_form_q = []  
unc_form_g = []

for wpoint in wpoint_array:
    for i in range(0,6):
            unc_form_q_pt = []
            unc_form_g_pt = []

            min = bin[i]
            max = bin[i+1]
            print(min)

            nominal_jets =  ReadingDijetMC(ntrackall4,str(min),var,"mc")
            
            higher_quark = nominal_jets["hq"]
            lower_quark = nominal_jets["lq"]
            higher_gluon = nominal_jets["hg"]
            lower_gluon = nominal_jets["lg"]
            higher_other = nominal_jets["ho"]
            lower_other = nominal_jets["lo"]        
            higher_mc = nominal_jets["hmc"]
            lower_mc = nominal_jets["lmc"]    

            data_jets = ReadingDijetData(ntrackall3,str(min),var,"data")

            higher_data_TH1I = data_jets["hd"]
            higher_data = higher_quark.Clone(higher_data_TH1I.GetName()+"TH1F")
            for j in range(1,higher_data.GetNbinsX()+1):
                    higher_data.SetBinContent(j,float(higher_data_TH1I.GetBinContent(j)))
                    higher_data.SetBinError(j,float(higher_data_TH1I.GetBinError(j)))             
            
            lower_data_TH1I = data_jets["ld"]
            lower_data = lower_quark.Clone(lower_data_TH1I.GetName()+"TH1F")
            for j in range(1,higher_data.GetNbinsX()+1):
                    lower_data.SetBinContent(j,float(lower_data_TH1I.GetBinContent(j)))
                    lower_data.SetBinError(j,float(lower_data_TH1I.GetBinError(j))) 
                    
            
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
            
            
            #Reading from other MC sample
            sherpa_jets =  ReadingDijetMC(ntrackall,str(min),var,"sherpa")
            
            higher_quark_sherpa = sherpa_jets["hq"]
            lower_quark_sherpa = sherpa_jets["lq"]
            higher_gluon_sherpa = sherpa_jets["hg"]
            lower_gluon_sherpa = sherpa_jets["lg"]
            higher_other_sherpa = sherpa_jets["ho"]
            lower_other_sherpa = sherpa_jets["lo"]        
            higher_sherpa = sherpa_jets["hmc"]
            lower_sherpa = sherpa_jets["lmc"]
            

            herang_jets =  ReadingDijetMC(fherang,str(min),var,"herang")
            
            higher_quark_herang = herang_jets["hq"]
            lower_quark_herang = herang_jets["lq"]
            higher_gluon_herang = herang_jets["hg"]
            lower_gluon_herang = herang_jets["lg"]
            higher_other_herang = herang_jets["ho"]
            lower_other_herang = herang_jets["lo"]        
            higher_herang = herang_jets["hmc"]
            lower_herang = herang_jets["lmc"]  

            herdipo_jets =  ReadingDijetMC(fherdipo,str(min),var,"herdipo")
            
            higher_quark_herdipo = herdipo_jets["hq"]
            lower_quark_herdipo = herdipo_jets["lq"]
            higher_gluon_herdipo = herdipo_jets["hg"]
            lower_gluon_herdipo = herdipo_jets["lg"]
            higher_other_herdipo = herdipo_jets["ho"]
            lower_other_herdipo = herdipo_jets["lo"]        
            higher_herdipo = herdipo_jets["hmc"]
            lower_herdipo = herdipo_jets["lmc"]  

            powpyt_jets =  ReadingDijetMC(powpyt,str(min),var,"powpyt")
            
            higher_quark_powpyt = powpyt_jets["hq"]
            lower_quark_powpyt = powpyt_jets["lq"]
            higher_gluon_powpyt = powpyt_jets["hg"]
            lower_gluon_powpyt = powpyt_jets["lg"]
            higher_other_powpyt = powpyt_jets["ho"]
            lower_other_powpyt = powpyt_jets["lo"]        
            higher_powpyt = powpyt_jets["hmc"]
            lower_powpyt = powpyt_jets["lmc"]  


            lund_jets =  ReadingDijetMC(sherpa_lund,str(min),var,"lund")
            
            higher_quark_lund = lund_jets["hq"]
            lower_quark_lund = lund_jets["lq"]
            higher_gluon_lund = lund_jets["hg"]
            lower_gluon_lund = lund_jets["lg"]
            higher_other_lund = lund_jets["ho"]
            lower_other_lund = lund_jets["lo"]        
            higher_lund = lund_jets["hmc"]
            lower_lund = lund_jets["lmc"]             

            
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
            mc_bin,data_bin,sf_q,sf_g,q_eff_mc,g_rej_mc,q_eff_data,g_rej_data=wp_sf_nominal(wpoint,nominal_extract_Q,nominal_extract_G,extracted_data_nominal_Q,extracted_data_nominal_G)
            
            
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

                    # get extracted data Q/G with sherpa sample
                    for j in range(1,higher_quark.GetNbinsX()+1):
                            F_data = forward_data_strap.GetBinContent(j)
                            C_data = central_data_strap.GetBinContent(j)
                            Q_data = -(C_data*fg-F_data*cg)/(cg*fq-fg*cq)
                            G_data = (C_data*fq-F_data*cq)/(cg*fq-fg*cq)
                            bootstrap_extracted_Q.SetBinContent(j,Q_data)
                            bootstrap_extracted_G.SetBinContent(j,G_data)
                    
                    sf_q_bootstrap,sf_g_bootstrap,q_eff_mc,g_rej_mc,q_eff_data,g_rej_data = wp_sf(mc_bin,nominal_extract_Q,nominal_extract_G,bootstrap_extracted_Q,bootstrap_extracted_G)

                    
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
            
            
            ##Matrix element uncertainty: pythia - powheg+pythia

            fg_powpyt,cg_powpyt,fq_powpyt,cq_powpyt = fraction(lower_quark_powpyt,lower_gluon_powpyt,higher_quark_powpyt,higher_gluon_powpyt)
            
            if (higher_quark_powpyt.Integral() != 0):
                higher_quark_powpyt.Scale(1./higher_quark_powpyt.Integral())
            if(higher_gluon_powpyt.Integral() != 0):
                higher_gluon_powpyt.Scale(1./higher_gluon_powpyt.Integral())
            if(lower_quark_powpyt.Integral() != 0):
                lower_quark_powpyt.Scale(1./lower_quark_powpyt.Integral())
            if(lower_gluon_powpyt.Integral() != 0):
                lower_gluon_powpyt.Scale(1./lower_gluon_powpyt.Integral())
            if(higher_powpyt.Integral() != 0):
                higher_powpyt.Scale(1./higher_powpyt.Integral())            
            if(lower_powpyt.Integral() != 0):
                lower_powpyt.Scale(1./lower_powpyt.Integral())    
                

            powpyt_extract_Q,powpyt_extract_G = mc_matrixmethod(lower_quark_powpyt,lower_gluon_powpyt,higher_quark_powpyt,higher_gluon_powpyt,fg_powpyt,cg_powpyt,fq_powpyt,cq_powpyt,higher_powpyt,lower_powpyt)            
            extracted_data_powpyt_Q,extracted_data_powpyt_G = data_matrixmethod(lower_quark_powpyt,lower_gluon_powpyt,higher_quark_powpyt,higher_gluon_powpyt,higher_data,lower_data,fg_powpyt,cg_powpyt,fq_powpyt,cq_powpyt)            

            powpyt_extract_G = error_check(powpyt_extract_G)
            powpyt_extract_Q = error_check(powpyt_extract_Q)
               
            
            sf_q_powpyt,sf_g_powpyt,q_eff_powpyt,g_rej_powpyt,q_eff_data_powpyt,g_rej_data_powpyt = wp_sf(mc_bin,powpyt_extract_Q,powpyt_extract_G,extracted_data_powpyt_Q,extracted_data_powpyt_G)
           
            
            me_unc_Q = abs(sf_q_powpyt-sf_q)
            me_unc_G = abs(sf_g_powpyt-sf_g)
            
            ##hadronization, difference in sherpa
            
            fg_lund,cg_lund,fq_lund,cq_lund = fraction(higher_quark_lund,higher_gluon_lund,higher_quark_lund,higher_gluon_lund)
            
            if (higher_quark_lund.Integral() != 0):
                higher_quark_lund.Scale(1./higher_quark_lund.Integral())
            if(higher_gluon_lund.Integral() != 0):
                higher_gluon_lund.Scale(1./higher_gluon_lund.Integral())
            if(higher_quark_lund.Integral() != 0):
                higher_quark_lund.Scale(1./higher_quark_lund.Integral())
            if(higher_gluon_lund.Integral() != 0):
                higher_gluon_lund.Scale(1./higher_gluon_lund.Integral())
            if(higher_lund.Integral() != 0):
                higher_lund.Scale(1./higher_lund.Integral())            
            if(lower_lund.Integral() != 0):
                lower_lund.Scale(1./lower_lund.Integral())    

            lund_extract_Q,lund_extract_G = mc_matrixmethod(lower_quark_lund,lower_gluon_lund,higher_quark_lund,higher_gluon_lund,fg_lund,cg_lund,fq_lund,cq_lund,higher_lund,lower_lund)            
            extracted_data_lund_Q,extracted_data_lund_G = data_matrixmethod(lower_quark_lund,lower_gluon_lund,higher_quark_lund,higher_gluon_lund,higher_data,lower_data,fg_lund,cg_lund,fq_lund,cq_lund)            

            lund_extract_G = error_check(lund_extract_G)
            lund_extract_Q = error_check(lund_extract_Q)
               

            fg_sherpa,cg_sherpa,fq_sherpa,cq_sherpa = fraction(lower_quark_sherpa,lower_gluon_sherpa,higher_quark_sherpa,higher_gluon_sherpa)
            
            if (higher_quark_sherpa.Integral() != 0):
                higher_quark_sherpa.Scale(1./higher_quark_sherpa.Integral())
            if(higher_gluon_sherpa.Integral() != 0):
                higher_gluon_sherpa.Scale(1./higher_gluon_sherpa.Integral())
            if(lower_quark_sherpa.Integral() != 0):
                lower_quark_sherpa.Scale(1./lower_quark_sherpa.Integral())
            if(lower_gluon_sherpa.Integral() != 0):
                lower_gluon_sherpa.Scale(1./lower_gluon_sherpa.Integral())
            if(higher_sherpa.Integral() != 0):
                higher_sherpa.Scale(1./higher_sherpa.Integral())            
            if(lower_sherpa.Integral() != 0):
                lower_sherpa.Scale(1./lower_sherpa.Integral())    
                

            sherpa_extract_Q,sherpa_extract_G = mc_matrixmethod(lower_quark_sherpa,lower_gluon_sherpa,higher_quark_sherpa,higher_gluon_sherpa,fg_sherpa,cg_sherpa,fq_sherpa,cq_sherpa,higher_sherpa,lower_sherpa)            
            extracted_data_sherpa_Q,extracted_data_sherpa_G = data_matrixmethod(lower_quark_sherpa,lower_gluon_sherpa,higher_quark_sherpa,higher_gluon_sherpa,higher_data,lower_data,fg_sherpa,cg_sherpa,fq_sherpa,cq_sherpa)            

            sherpa_extract_G = error_check(sherpa_extract_G)
            sherpa_extract_Q = error_check(sherpa_extract_Q)
               
            
            sf_q_sherpa,sf_g_sherpa,q_eff_sherpa,g_rej_sherpa,q_eff_data_sherpa,g_rej_data_sherpa = wp_sf(mc_bin,sherpa_extract_Q,sherpa_extract_G,extracted_data_sherpa_Q,extracted_data_sherpa_G)
            sf_q_lund,sf_g_lund,q_eff_lund,g_rej_lund,q_eff_data_lund,g_rej_data_lund = wp_sf(mc_bin,lund_extract_Q,lund_extract_G,extracted_data_lund_Q,extracted_data_lund_G)
            
            # Showering

            had_unc_Q = abs(sf_q_lund-sf_q_sherpa)
            had_unc_G = abs(sf_g_lund-sf_g_sherpa)


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
               
            
            sf_q_herang,sf_g_herang,q_eff_herang,g_rej_herang,q_eff_data_herang,g_rej_data_herang =  wp_sf(mc_bin,herang_extract_Q,herang_extract_G,extracted_data_herang_Q,extracted_data_herang_G)
            

            
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
            
            
            sf_q_herdipo,sf_g_herdipo,q_eff_herdipo,g_rej_herdipo,q_eff_data_herdipo,g_rej_data_herdipo = wp_sf(mc_bin,herdipo_extract_Q,herdipo_extract_G,extracted_data_herdipo_Q,extracted_data_herdipo_G)
            
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
                    pdf_jets =  ReadingDijetMCindex(ntrackall5,str(min),var,"pdf"+str(k),k)
            
                    higher_quark_pdf = pdf_jets["hq"]
                    lower_quark_pdf = pdf_jets["lq"]
                    higher_gluon_pdf = pdf_jets["hg"]
                    lower_gluon_pdf = pdf_jets["lg"]
                    higher_other_pdf = pdf_jets["ho"]
                    lower_other_pdf = pdf_jets["lo"]        
                    higher_pdf = pdf_jets["hmc"]
                    lower_pdf = pdf_jets["lmc"]  
                    
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
                            
                    #for j in range(1,higher_quark.GetNbinsX()+1):
                            #print(higher_pdf.GetBinContent(j))
                    pdf_extract_Q,pdf_extract_G = mc_matrixmethod(lower_quark_pdf,lower_gluon_pdf,higher_quark_pdf,higher_gluon_pdf,fg_pdf,cg_pdf,fq_pdf,cq_pdf,higher_pdf,lower_pdf)
                    extracted_data_pdf_Q,extracted_data_pdf_G = data_matrixmethod(lower_quark_pdf,lower_gluon_pdf,higher_quark_pdf,higher_gluon_pdf,higher_data,lower_data,fg_pdf,cg_pdf,fq_pdf,cq_pdf)
                    sf_q_pdf,sf_g_pdf,q_eff_pdf,g_rej_pdf,q_eff_data_pdf,g_rej_data_pdf = wp_sf(mc_bin,pdf_extract_Q,pdf_extract_G,extracted_data_pdf_Q,extracted_data_pdf_G)
                    sf_pdf_qvals[k] = sf_q_pdf
                    sf_pdf_gvals[k] = sf_g_pdf

            sf_pdf_qvals.sort()
            sf_pdf_gvals.sort()
            Q = np.median(sf_pdf_qvals)
            G = np.median(sf_pdf_gvals)
            pdf_sigmaQ = .5*(sf_pdf_qvals[int(.84*len(sf_pdf_qvals))] - sf_pdf_qvals[int(.16*len(sf_pdf_qvals))])
            pdf_sigmaG = .5*(sf_pdf_gvals[int(.84*len(sf_pdf_gvals))] - sf_pdf_gvals[int(.16*len(sf_pdf_gvals))])
            pdf_unc_Q = pdf_sigmaQ
            pdf_unc_G = pdf_sigmaG
            
            ### scale variation ### 
            
            for k in range(0,7):
                sv_jets =  ReadingDijetMCindex(sv,str(min),var,"sv1"+str(k),k)
            
                higher_quark_1 = sv_jets["hq"]
                lower_quark_1 = sv_jets["lq"]
                higher_gluon_1 = sv_jets["hg"]
                lower_gluon_1 = sv_jets["lg"]
                higher_other_1 = sv_jets["ho"]
                lower_other_1 = sv_jets["lo"]        
                higher_1 = sv_jets["hmc"]
                lower_1 = sv_jets["lmc"]  
                
                fg_1,cg_1,fq_1,cq_1 = fraction(lower_quark_1,lower_gluon_1,higher_quark_1,higher_gluon_1)


                if (higher_quark_1.Integral() != 0):
                    higher_quark_1.Scale(1./higher_quark_1.Integral())
                if(higher_gluon_1.Integral() != 0):
                    higher_gluon_1.Scale(1./higher_gluon_1.Integral())
                if(lower_quark_1.Integral() != 0):
                    lower_quark_1.Scale(1./lower_quark_1.Integral())
                if(lower_gluon_1.Integral() != 0):
                    lower_gluon_1.Scale(1./lower_gluon_1.Integral())
                if (higher_1.Integral() != 0):
                    higher_1.Scale(1./higher_1.Integral())                    
                if(lower_1.Integral() != 0):
                    lower_1.Scale(1./lower_1.Integral()) 
                    
                SV_extract_Q,SV_extract_G = mc_matrixmethod(lower_quark_1,lower_gluon_1,higher_quark_1,higher_gluon_1,fg_1,cg_1,fq_1,cq_1,higher_1,lower_1)
                extracted_data_SV_Q,extracted_data_SV_G = data_matrixmethod(lower_quark_1,lower_gluon_1,higher_quark_1,higher_gluon_1,higher_data,lower_data,fg_1,cg_1,fq_1,cq_1)
                sf_q_SV,sf_g_SV,q_eff_SV,g_rej_SV,q_eff_data_SV,g_rej_data_SV = wp_sf(mc_bin,SV_extract_Q,SV_extract_G,extracted_data_SV_Q,extracted_data_SV_G)
                
                if k == 0: 
                    sf_q_SV0 = sf_q_SV
                    sf_g_SV0 = sf_g_SV   
                    
                if k == 1:
                    SV_unc_Q = abs(sf_q_SV-sf_q_SV0)
                    SV_unc_G = abs(sf_g_SV-sf_g_SV0)

                if k > 1:
                    diffq = abs(sf_q_SV-sf_q_SV0)
                    diffg = abs(sf_g_SV-sf_g_SV0)

                    if(diffq>SV_unc_Q):
                        SV_unc_Q = diffq
                    if(diffg>SV_unc_G):
                        SV_unc_G = diffg
                        
            #mc closure
            sf_q_closure,sf_g_closure,q_eff_closure,g_rej_closure,q_eff_data_closure,g_rej_data_closure = wp_sf(mc_bin,higher_quark,higher_gluon,extracted_data_sherpa_Q,extracted_data_sherpa_G)
            closure_unc_Q = abs(sf_q_closure-sf_q)
            closure_unc_G = abs(sf_g_closure-sf_g)            
            
            #calculate the SF uncertainty
            qeff.append(q_eff_mc)
            grej.append(g_rej_mc)
            qsf_array.append(sf_q)
            gsf_array.append(sf_g)


            unc_q = np.sqrt(closure_unc_Q**2+had_unc_Q**2+show_unc_Q**2+pdf_unc_Q**2+SV_unc_Q**2+sigmaQ**2+RF_unc_Q**2+me_unc_Q**2)
            unc_g = np.sqrt(closure_unc_G**2+had_unc_G**2+show_unc_G**2+pdf_unc_G**2+SV_unc_G**2+sigmaG**2+RF_unc_G**2+me_unc_G**2)



            histq.SetBinContent(i+1,q_eff_mc)
            histg.SetBinContent(i+1,g_rej_mc)
            histq.SetBinError(i+1,Q_eff_unc)
            histg.SetBinError(i+1,G_rej_unc) 
            histqsf.SetBinContent(i+1,sf_q)
            histgsf.SetBinContent(i+1,sf_g)
            histqsf.SetBinError(i+1,unc_q)
            histgsf.SetBinError(i+1,unc_g)    

            pt_range = str(str(min)+"-"+str(max))
      

            unc_form_q_pt = [pt_range,round(sf_q,5),round(unc_q,5),round(sigmaQ,5),round(closure_unc_Q,5),round(had_unc_Q,5),round(me_unc_Q,5),round(show_unc_Q,5),round(pdf_unc_Q,5),round(SV_unc_Q,5),round(RF_unc_Q,5)]
            unc_form_g_pt = [pt_range,round(sf_g,5),round(unc_g,5),round(sigmaG,5),round(closure_unc_G,5),round(had_unc_G,5),round(me_unc_G,5),round(show_unc_G,5),round(pdf_unc_G,5),round(SV_unc_G,5),round(RF_unc_G,5)]


            unc_form_q.append(unc_form_q_pt)
            unc_form_g.append(unc_form_g_pt)

            

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
    
    c.Print("dijet-wp-" + str(wpoint) +"-"+ str(var)+ "-rej.pdf")
    
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
    
    c1.Print("dijet-sf-" + str(wpoint) +"-"+ str(var)+ "-rej.pdf")
    
    tbq = PrettyTable()

    tbq.field_names = ["Pt(GeV)","Quark Scale Factor","Uncertainty","Statistical","MC Closure","Hardonization","Matrix Element","Showering","PDF","Scale Variation","Reweight Factor"]
    
    for i in range(6):
        tbq.add_row(unc_form_q[i])
    print(tbq)
    
    tbg = PrettyTable()

    tbg.field_names = ["Pt(GeV)","Gluon Scale Factor","Uncertainty","Statistical","MC Closure","Hardonization","Matrix Element","Showering","PDF","Scale Variation","Reweight Factor"]
    
    for i in range(6):
        tbg.add_row(unc_form_g[i])
    print(tbg)
    foutput = open(str(wpoint) +"-"+ str(var)+".txt", "w")
    foutput.write(str(wpoint) + " " + str(var) +" "+ str(doreweight) +  "\n")
    foutput.write(tbq.get_string())
    foutput.write("\n")
    foutput.write(tbg.get_string())
    foutput.close()
        
    foutputq = open("Quark"+str(wpoint) +"-"+ str(var)+".csv", "w")
    

    foutputq.write("Pt(GeV),Quark Scale Factor,Uncertainty,Statistical,MC closure,Hardonization,Matrix Element,Showering,PDF,Scale Variation,Reweight Factor\n")
    
    for i in range(6):
        for j in (range(len(unc_form_q[i]))):
            foutputq.write(str(unc_form_q[i][j]))
            if j != (len(unc_form_q[i]))-1:
                        foutputq.write(",")
        foutputq.write("\n")
        
    foutputg = open("Gluon"+str(wpoint) +"-"+ str(var)+".csv", "w")
    

    foutputg.write("Pt(GeV),Gluon Scale Factor,Uncertainty,Statistical,MC closure,Hardonization,Matrix Element,Showering,,PDF,Scale Variation,Reweight Factor\n")
    
    for i in range(6):
        for j in (range(len(unc_form_g[i]))):
            foutputg.write(str(unc_form_g[i][j]))
            if j != (len(unc_form_g[i]))-1:
                        foutputg.write(",")
        foutputg.write("\n")
