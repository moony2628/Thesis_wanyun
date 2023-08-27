# uncertainties for statistics with pyroot ReadingTree output pdf
from ROOT import *
import numpy as np
import uncertainties as unc # propagate uncertainties
from uncertainties import unumpy as unp # array operations for type ufloat
import re

doreweight = 0   #decide if we want to do the reweighting process

var = "bdt"  #change the var name according to the inputvar you want to read
mc = "pythia"   #by setting it as "SF" or "MC", it will automatically making scale factor plots or MC closure plots
inputvar = var  #by setting it as bdt (or ntrk,width,c1..), it will read the corresponding histogram, but remember to change the TLine range according to X-axis of different variable, one can check it by browsing the histograms in root file.



def rebin(input_hist):
    hist = input_hist.Clone()
    a = np.array([])
    for j in range(1,hist.GetNbinsX()+2):
        if var == "ntrk":
            if j == 1 or (j >= 6 and j<=30 and (j+5)%4 == 0) or j==31 or j ==41 or j ==61: 
                a=np.append(a,hist.GetBinLowEdge(j))
        if var == "bdt":
            if j == 1 or (j >= 13 and j<=49 and (j-1)%4 == 0) or j==61: 
                a=np.append(a,hist.GetBinLowEdge(j))

    hist = hist.Rebin(len(a)-1,"",a)
    return(hist)

# input two eta region to be compared
#eta_bin = ["0-0.5_0.5-1","0.5-1_1-2.1"]
eta_bin = ["0-2.1_0-2.1","0-0.5_0.5-2.1"]


ntrackall = TFile("../roots/pythia_etas_all.root") #mc sample
ntrackall3  = TFile("../roots/data_etas_prescale.root") #data sample



def myText(x,y,text,color =1):
    l = TLatex()
    l.SetTextSize(0.025)
    l.SetNDC()
    l.SetTextColor(color)
    l.DrawLatex(x,y,text)
    pass


#convert histogram and error into unp.uarray
#if sample is pyroot input, GetBinError returns correct result
def unc_array(hist):
    value = np.zeros(hist.GetNbinsX())
    error = np.zeros(hist.GetNbinsX())
    for j in range(1,hist.GetNbinsX()+1):
        value[j-1] = hist.GetBinContent(j)
        error[j-1] = hist.GetBinError(j)
    result = unp.uarray(value,error)
    return(result)

#convert histogram and error into unp.uarray
#if sample is uproot input, err is the sumw2 of the corresponding histogram
def unc_array_err(hist,err):
    value = np.zeros(hist.GetNbinsX())
    error = np.zeros(hist.GetNbinsX())
    for j in range(1,hist.GetNbinsX()+1):
        value[j-1] = hist.GetBinContent(j)
        error[j-1] = np.sqrt(err.GetBinContent(j))
    result = unp.uarray(value,error)
    return(result)

def set_hist_error(hist,unc):
    for i in range(1,hist.GetNbinsX()+1):
        hist.SetBinError(i,unp.std_devs(unc[i-1]))
        
def etatovalue(etastring):
    pattern  = re.compile("^([0-9.]+)-([0-9.]+)_([0-9.]+)-([0-9.]+)$")
    result = pattern.match(etastring)
    return(result.group(1),result.group(2),result.group(3),result.group(4))

eta_start_central_1,eta_end_central_1,eta_start_forward_1,eta_end_forward_1 = etatovalue(eta_bin[0])
eta_start_central_2,eta_end_central_2,eta_start_forward_2,eta_end_forward_2 = etatovalue(eta_bin[1])


bin = [0,50,100,150,200,300,400,500,600,800,1000,1200,1500,2000]
for k in range(7,13):   #for only dijet event, start from jet pT>500 GeV
    quark_hist_array = []
    gluon_hist_array = []
    quark_hist_array_truth = []
    gluon_hist_array_truth = []
    
    quark_hist_array_ratio = []
    gluon_hist_array_ratio = []
    for e in range(len(eta_bin)):            

        min = bin[k]
        max = bin[k+1]
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

        c = TCanvas("c","c",800,800)

        #add leading and subleading jet from only dijet event together,
        #note that for gammajet+dijet event, we need to add leading jet from gammajet and leading jet from dijet sample together
        higher_data.Add(higher_data2)
        lower_data.Add(lower_data2)
        higher_quark.Add(higher_quark2)
        higher_gluon.Add(higher_gluon2)
        lower_quark.Add(lower_quark2)
        lower_gluon.Add(lower_gluon2)
        
        

        higher_mc = higher_quark.Clone()
        higher_mc.Add(higher_gluon)
        lower_mc = lower_quark.Clone()
        lower_mc.Add(lower_gluon)
    
        higher_quark2_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Forward_Quark_"+inputvar + "_err")
        higher_gluon2_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Forward_Gluon_"+inputvar + "_err")
        lower_quark2_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Central_Quark_"+inputvar + "_err")
        lower_gluon2_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_LeadingJet_Central_Gluon_"+inputvar + "_err")

        higher_quark_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Forward_Quark_"+inputvar + "_err")
        higher_gluon_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Forward_Gluon_"+inputvar + "_err")

        lower_quark_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Central_Quark_"+inputvar + "_err")
        lower_gluon_err = ntrackall.Get(str(eta_bin[e])+"_"+str(min)+"_SubJet_Central_Gluon_"+inputvar + "_err")



        #add leading and subleading jet from only dijet event together,
        #note that for gammajet+dijet event, we need to add leading jet from gammajet and leading jet from dijet sample together

        higher_quark_err.Add(higher_quark2_err)
        higher_gluon_err.Add(higher_gluon2_err)
        lower_quark_err.Add(lower_quark2_err)
        lower_gluon_err.Add(lower_gluon2_err)


 
        
        
        higher_mc_err = higher_quark_err.Clone()
        higher_mc_err.Add(higher_gluon_err)
        lower_mc_err = lower_quark_err.Clone()
        lower_mc_err.Add(lower_gluon_err)
        
        
        #uncertainty propagation
        higher_quark_unc = unc_array_err(higher_quark,higher_quark_err)
        higher_gluon_unc = unc_array_err(higher_gluon,higher_gluon_err)
        lower_quark_unc = unc_array_err(lower_quark,lower_quark_err)
        lower_gluon_unc = unc_array_err(lower_gluon,lower_gluon_err)
        higher_mc_unc = unc_array_err(higher_mc,higher_mc_err)
        lower_mc_unc = unc_array_err(lower_mc,lower_mc_err)
        
        higher_data_unc = unc_array(higher_data)
        lower_data_unc = unc_array(lower_data)
        
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
        higher_data_unc = higher_data_unc/higher_data_unc.sum()
        lower_data_unc = lower_data_unc/lower_data_unc.sum()
        higher_mc_unc = higher_mc_unc/higher_mc_unc.sum()
        lower_mc_unc = lower_mc_unc/lower_mc_unc.sum()        
        
        if (doreweight=="Quark"):
                for i in range(1,higher_quark.GetNbinsX()+1):
                        if (lower_quark.GetBinContent(i) > 0 and lower_gluon.GetBinContent(i) > 0):
                                factor_gluon_unc[i-1] = higher_gluon_unc[i-1]/lower_gluon_unc[i-1]
                                factor_quark_unc[i-1] = higher_quark_unc[i-1]/lower_quark_unc[i-1]
                        else:
                                factor_gluon_unc[i-1] = unc.ufloat(1, 0)
                                factor_quark_unc[i-1] = unc.ufloat(1, 0)
                lower_quark_unc=lower_quark_unc*factor_quark_unc
                lower_gluon_unc=lower_gluon_unc*factor_quark_unc
                lower_mc_unc=lower_mc_unc*factor_quark_unc
                lower_data_unc=lower_data_unc*factor_quark_unc

        if (doreweight=="Gluon"):
                for i in range(1,higher_quark.GetNbinsX()+1):
                        if (lower_quark.GetBinContent(i) > 0 and lower_gluon.GetBinContent(i) > 0):
                                #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                factor_gluon_unc[i-1] = higher_gluon_unc[i-1]/lower_gluon_unc[i-1]
                                factor_quark_unc[i-1] = higher_quark_unc[i-1]/lower_quark_unc[i-1]
                        else:
                                factor_gluon_unc[i-1] = unc.ufloat(1, 0)
                                factor_quark_unc[i-1] = unc.ufloat(1, 0)
                lower_quark_unc=lower_quark_unc*factor_gluon_unc
                lower_gluon_unc=lower_gluon_unc*factor_gluon_unc
                lower_mc_unc=lower_mc_unc*factor_quark_unc
                lower_data_unc=lower_data_unc*factor_gluon_unc
            
        higher_quark_unc = higher_quark_unc/higher_quark_unc.sum()
        higher_gluon_unc = higher_gluon_unc/higher_gluon_unc.sum()
        lower_quark_unc = lower_quark_unc/lower_quark_unc.sum()
        lower_gluon_unc = lower_gluon_unc/lower_gluon_unc.sum()
        higher_data_unc = higher_data_unc/higher_data_unc.sum()
        lower_data_unc = lower_data_unc/lower_data_unc.sum()
        higher_mc_unc = higher_mc_unc/higher_mc_unc.sum()
        lower_mc_unc = lower_mc_unc/lower_mc_unc.sum()   
        
        higher = higher_mc.Clone("")
        lower = lower_mc.Clone("")

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

        
        
        # calculate the value(not correct error)
        ToT_Fq2 = 0.
        ToT_Fg2 = 0.

        ToT_Cq2 =0.
        ToT_Cg2 = 0.

        for j in range(1,lower_quark.GetNbinsX()+1):
            ToT_Fq2+=higher_quark.GetBinContent(j)
            ToT_Cq2+=lower_quark.GetBinContent(j)
            ToT_Fg2+=higher_gluon.GetBinContent(j)
            ToT_Cg2+=lower_gluon.GetBinContent(j)

        # calculate the fraction of forward(higher) / central(lower) quark or gluon jet
        fg=ToT_Fg2/(ToT_Fg2+ToT_Fq2)
        cg=ToT_Cg2/(ToT_Cq2+ToT_Cg2)
        fq=1.-fg
        cq=1.-cg

        if (lower_quark.Integral() != 0):
                lower_quark.Scale(1./lower_quark.Integral())
        if(lower_gluon.Integral() != 0):
                lower_gluon.Scale(1./lower_gluon.Integral())
        if(higher_quark.Integral() != 0):
                higher_quark.Scale(1./higher_quark.Integral())
        if(higher_gluon.Integral() != 0):
                higher_gluon.Scale(1./higher_gluon.Integral())
        if(lower_data.Integral() != 0):
                lower_data.Scale(1./lower_data.Integral())
        if(higher_data.Integral() != 0):
                higher_data.Scale(1./higher_data.Integral())
        if(lower_mc.Integral() != 0):
                lower_mc.Scale(1./lower_mc.Integral())
        if(higher_mc.Integral() != 0):
                higher_mc.Scale(1./higher_mc.Integral())
                
        if (doreweight=="Quark"):
                for i in range(1,higher_quark.GetNbinsX()+1):
                        if (lower_quark.GetBinContent(i) > 0 and lower_gluon.GetBinContent(i) > 0):
                                #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                factor_gluon = higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                factor_quark = higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i)
                                lower_quark.SetBinContent(i,lower_quark.GetBinContent(i)*factor_quark)
                                lower_gluon.SetBinContent(i,lower_gluon.GetBinContent(i)*factor_quark)
                                lower_data.SetBinContent(i,lower_data.GetBinContent(i)*factor_quark)
                                lower_mc.SetBinContent(i,lower_mc.GetBinContent(i)*factor_quark)

            
        if (doreweight=="Gluon"):
                for i in range(1,higher_quark.GetNbinsX()+1):
                        if (lower_quark.GetBinContent(i) > 0 and lower_gluon.GetBinContent(i) > 0):
                                #print i,higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i),higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                factor_gluon = higher_gluon.GetBinContent(i)/lower_gluon.GetBinContent(i)
                                factor_quark = higher_quark.GetBinContent(i)/lower_quark.GetBinContent(i)
                                lower_quark.SetBinContent(i,lower_quark.GetBinContent(i)*factor_gluon)
                                lower_gluon.SetBinContent(i,lower_gluon.GetBinContent(i)*factor_gluon)
                                lower_data.SetBinContent(i,lower_data.GetBinContent(i)*factor_gluon)
                                lower_mc.SetBinContent(i,lower_mc.GetBinContent(i)*factor_gluon)


        if (lower_quark.Integral() != 0):
                lower_quark.Scale(1./lower_quark.Integral())
        if(lower_gluon.Integral() != 0):
                lower_gluon.Scale(1./lower_gluon.Integral())
        if(higher_quark.Integral() != 0):
                higher_quark.Scale(1./higher_quark.Integral())
        if(higher_gluon.Integral() != 0):
                higher_gluon.Scale(1./higher_gluon.Integral())
        if(lower_data.Integral() != 0):
                lower_data.Scale(1./lower_data.Integral())
        if(higher_data.Integral() != 0):
                higher_data.Scale(1./higher_data.Integral())
        if(lower_mc.Integral() != 0):
                lower_mc.Scale(1./lower_mc.Integral())
        if(higher_mc.Integral() != 0):
                higher_mc.Scale(1./higher_mc.Integral())

        higher = higher_mc.Clone("")
        lower = lower_mc.Clone("")

        #Now, let's solve.

        quark = higher_quark.Clone("quark")
        gluon = higher_gluon.Clone("gluon")
        quark_data = higher_data.Clone("")
        gluon_data = higher_data.Clone("") 
        for i in range(1,quark_data.GetNbinsX()+1):
            F = higher.GetBinContent(i)
            C = lower.GetBinContent(i)
            if((cg*fq-fg*cq) != 0 ):
                Q = -(C*fg-F*cg)/(cg*fq-fg*cq)
                G = (C*fq-F*cq)/(cg*fq-fg*cq)
                quark.SetBinContent(i,Q)
                gluon.SetBinContent(i,G)
                #print "   ",i,G,higher_gluon.GetBinContent(i),lower_gluon.GetBinContent(i)
        pass



        #lower_data.Scale(1./lower_data.Integral())
        #higher_data.Scale(1./higher_data.Integral())
        #quark_data = higher_data.Clone("")
        #gluon_data = higher_data.Clone("")

        for i in range(1,higher_data.GetNbinsX()+1):
                F = higher_data.GetBinContent(i)
                C = lower_data.GetBinContent(i)
                if((cg*fq-fg*cq) != 0):
                        Q = -(C*fg-F*cg)/(cg*fq-fg*cq)
                        G = (C*fq-F*cq)/(cg*fq-fg*cq)
                        quark_data.SetBinContent(i,Q)
                        gluon_data.SetBinContent(i,G)
                        #print "   ",i,"  ",G,"   ",Q
                pass

        # set error bar in histogram        
        set_hist_error(quark,Q_unc)
        set_hist_error(gluon,G_unc)
        set_hist_error(quark_data,Q_data_unc)
        set_hist_error(gluon_data,G_data_unc)
        set_hist_error(higher_quark,higher_quark_unc)
        set_hist_error(higher_gluon,higher_gluon_unc)

        ## below just do the ploting

        gPad.SetLeftMargin(0.15)
        gPad.SetTopMargin(0.05)
        gPad.SetBottomMargin(0.15)
 
        quark_ratio = quark.Clone("")
        gluon_ratio = gluon.Clone("")
        ex_quark_ratio = quark_data.Clone("")
        ex_gluon_ratio = gluon_data.Clone("")
        quark_ratio.GetYaxis().SetTitle("Extracted/Truth")
        gluon_ratio.GetYaxis().SetTitle("Extracted/Truth")

        ex_quark_ratio.Divide(quark)
        ex_gluon_ratio.Divide(gluon)
        quark_ratio.Divide(higher_quark)
        gluon_ratio.Divide(higher_gluon)
        
        #if k == 7:
        #    for j in range(1,quark.GetNbinsX()+1):
        #        print(j,gluon_ratio.GetBinContent(j),higher_gluon.GetBinContent(j),gluon.GetBinContent(j))
        gStyle.SetOptStat(0)
    ######################## for ratio plot

        quark_hist_array.append(quark)
        quark_hist_array_truth.append(higher_quark)
        gluon_hist_array.append(gluon)
        gluon_hist_array_truth.append(higher_gluon)        
        quark_hist_array_ratio.append(quark_ratio)
        gluon_hist_array_ratio.append(gluon_ratio)
    
    fg_quark_mc_ratio = quark_hist_array[0].Clone()
    fg_quark_mc_ratio.Divide(quark_hist_array[1])
    fg_quark_truth_ratio = quark_hist_array_truth[0].Clone()
    fg_quark_truth_ratio.Divide(quark_hist_array_truth[1])
    fg_gluon_mc_ratio = gluon_hist_array[0].Clone()
    fg_gluon_mc_ratio.Divide(gluon_hist_array[1])
    fg_gluon_truth_ratio = gluon_hist_array_truth[0].Clone()
    fg_gluon_truth_ratio.Divide(gluon_hist_array_truth[1])
    
    c.Divide(3,1)
    top = c.cd(1)
    top.SetPad(0.0,0.0,1.0,1.0)
    top.SetFillColor(0)
    top.SetBorderMode(0)
    top.SetBorderSize(2)
    top.SetTickx(1)
    top.SetTicky(1)
    top.SetLeftMargin(0.14)
    top.SetRightMargin(0.055)
    top.SetBottomMargin(0.4)#0.25
    top.SetFrameBorderMode(0)
    #top.SetLogy(1)
    middle = c.cd(2)
    middle.SetPad(0.0,0.15,1.0,0.4)
    middle.SetFillColor(0)
    middle.SetBorderMode(0)
    middle.SetBorderSize(2)
    middle.SetTickx(1)
    middle.SetTicky(1)
    middle.SetLeftMargin(0.14)
    middle.SetRightMargin(0.055)
    middle.SetTopMargin(0.045)
    middle.SetBottomMargin(0.4)
    middle.SetFrameBorderMode(0)
    
    bot = c.cd(3)
    bot.SetPad(0.0,0.0,1.0,0.25)
    bot.SetFillColor(0)
    bot.SetBorderMode(0)
    bot.SetBorderSize(2)
    bot.SetTickx(1)
    bot.SetTicky(1)
    bot.SetLeftMargin(0.14)
    bot.SetRightMargin(0.055)
    bot.SetTopMargin(0)
    bot.SetBottomMargin(0.4)
    bot.SetFrameBorderMode(0)
    
    quark_hist_array[0].SetTitle("")
    quark_hist_array[0].GetXaxis().SetTitle(var)
    quark_hist_array[0].GetYaxis().SetTitle("Normalized to unity")
    quark_hist_array[0].GetYaxis().SetNdivisions(505)
    quark_hist_array[0].GetYaxis().SetRangeUser(-0.01,quark.GetMaximum()*1.5)
    
    for j in range(2):
            quark_hist_array[j].SetMarkerColor(3*j+1)
            quark_hist_array[j].SetLineColor(3*j+1)
            quark_hist_array[j].SetMarkerSize(0.5)
            quark_hist_array[j].SetLineStyle(1)
  
    for j in range(2):
            quark_hist_array_truth[j].SetMarkerColor(3*j+1)
            quark_hist_array_truth[j].SetLineColor(3*j+1)
            quark_hist_array_truth[j].SetMarkerSize(0.5)
            quark_hist_array_truth[j].SetLineStyle(2)
      
    for j in range(2):
            quark_hist_array_ratio[j].SetMarkerColor(3*j+1)
            quark_hist_array_ratio[j].SetLineColor(3*j+1)
            quark_hist_array_ratio[j].SetMarkerSize(0.5)
            quark_hist_array_ratio[j].SetLineStyle(2)
   
            
    quark_hist_array_ratio[0].SetTitle("")
    quark_hist_array_ratio[0].GetYaxis().SetRangeUser(0.7,1.3)
    quark_hist_array_ratio[0].GetXaxis().SetTitleOffset(1)
    quark_hist_array_ratio[0].GetXaxis().SetTitleSize(0.11)
    quark_hist_array_ratio[0].GetXaxis().SetLabelSize(0.1)
    quark_hist_array_ratio[0].GetXaxis().SetLabelOffset(0.03)
    quark_hist_array_ratio[0].GetYaxis().SetTitleSize(0.1)
    quark_hist_array_ratio[0].GetYaxis().SetTitleOffset(0.5)
    quark_hist_array_ratio[0].GetYaxis().SetLabelOffset(0.01)
    quark_hist_array_ratio[0].GetYaxis().SetLabelSize(0.1)

    fg_quark_mc_ratio.GetYaxis().SetTitle("\eta_{1}/\eta_{2}")
    fg_quark_mc_ratio.GetYaxis().SetRangeUser(0.7,1.3)
    fg_quark_mc_ratio.GetXaxis().SetTitleOffset(1)
    fg_quark_mc_ratio.GetXaxis().SetTitleSize(0.11)
    fg_quark_mc_ratio.GetXaxis().SetLabelSize(0.1)
    fg_quark_mc_ratio.GetXaxis().SetLabelOffset(0.03)
    fg_quark_mc_ratio.GetYaxis().SetTitleSize(0.1)
    fg_quark_mc_ratio.GetYaxis().SetTitleOffset(0.5)
    fg_quark_mc_ratio.GetYaxis().SetLabelOffset(0.01) 
    fg_quark_mc_ratio.GetYaxis().SetLabelSize(0.1)

    fg_quark_mc_ratio.SetMarkerSize(0.5)
    fg_quark_mc_ratio.SetLineStyle(1)  
    fg_quark_mc_ratio.SetLineColor(1)

    fg_quark_truth_ratio.SetMarkerSize(0.5)
    fg_quark_truth_ratio.SetLineStyle(2)  
    fg_quark_truth_ratio.SetLineColor(1)
    
    top.cd()
    quark_hist_array[0].Draw("HIST e") #extracted mc
    quark_hist_array[1].Draw("HIST e same") 
    #quark_hist_array[2].Draw("HIST same") 
    quark_hist_array_truth[0].SetMarkerStyle(20)
    quark_hist_array_truth[1].SetMarkerStyle(20)

    quark_hist_array_truth[0].Draw("HIST same") #extracted mc
    quark_hist_array_truth[1].Draw("HIST same") 
    #quark_hist_array_truth[2].Draw("p same") 
    #higher_quark.Draw("HIST same") #truth mc
    if(inputvar == "ntrk"):
        line = TLine(0.,1,60,1)
        quark_hist_array_ratio[0].GetXaxis().SetTitle("N_{track}")
        leg = TLegend(0.50,0.7,0.9,0.90)
        leg.SetTextSize(0.018)
    if(inputvar == "bdt"):
        line = TLine(-0.8,1,0.7,1)
        quark_hist_array_ratio[0].GetXaxis().SetTitle("BDT score")
        leg = TLegend(0.6,0.7,0.9,0.9) ##0.6,0.5,0.9,0.7
        leg = TLegend(0.50,0.7,0.9,0.90)
        leg.SetTextSize(0.018)
    leg.SetTextFont(42)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(1)
    leg.AddEntry(quark_hist_array[0],"extracted quark(mc) \eta_{1}: "+eta_start_central_1+"<|\eta_{C}|<"+eta_end_central_1+","+eta_start_forward_1+"<|\eta_{C}|<"+eta_end_forward_1+"","l")
    leg.AddEntry(quark_hist_array[1],"extracted quark(mc) \eta_{2}: "+eta_start_central_2+"<|\eta_{C}|<"+eta_end_central_2+","+eta_start_forward_2+"<|\eta_{C}|<"+eta_end_forward_2+"","l")
    leg.AddEntry(quark_hist_array_truth[0],"quark(mc) \eta_{1}:"+eta_start_central_1+"<|\eta_{C}|<"+eta_end_central_1+","+eta_start_forward_1+"<|\eta_{C}|<"+eta_end_forward_1+"","l")
    leg.AddEntry(quark_hist_array_truth[1],"quark(mc) \eta_{2}:"+eta_start_central_2+"<|\eta_{C}|<"+eta_end_central_2+","+eta_start_forward_2+"<|\eta_{C}|<"+eta_end_forward_2+"","l")

    leg.Draw()
    myText(0.18,0.84,"#it{#bf{#scale[1.5]{#bf{ATLAS} Internal}}}")
    myText(0.18,0.80,"#bf{#scale[1.2]{#sqrt{s} = 13 TeV}}")
    myText(0.18,0.75,"#bf{#scale[1.2]{pT range: "+str(min)+" - "+str(max)+" GeV}}")
    
    middle.cd()
    fg_quark_mc_ratio.Draw("HIST e") #extracted mc
    fg_quark_truth_ratio.Draw("HIST e same")   
    line.Draw("same")
    bot.cd()
    #quark_ratio.Draw("HIST")
    
    quark_hist_array_ratio[0].Draw("HIST e") #extracted mc
    quark_hist_array_ratio[1].Draw("HIST e same") 
    #quark_hist_array_ratio[2].Draw("HIST same") 
    line.Draw("same")
    c.Print("./plots_"+var+"/"+"quark_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"_compare_"+eta_bin[0]+"_"+eta_bin[1]+".pdf")
    
    
    
    
    c1 = TCanvas("c","c",800,800)
    c1.Divide(3,1)
    top = c1.cd(1)
    top.SetPad(0.0,0.0,1.0,1.0)
    top.SetFillColor(0)
    top.SetBorderMode(0)
    top.SetBorderSize(2)
    top.SetTickx(1)
    top.SetTicky(1)
    top.SetLeftMargin(0.14)
    top.SetRightMargin(0.055)
    top.SetBottomMargin(0.4)#0.25
    top.SetFrameBorderMode(0)
    #top.SetLogy(1)
    middle = c1.cd(2)
    middle.SetPad(0.0,0.15,1.0,0.4)
    middle.SetFillColor(0)
    middle.SetBorderMode(0)
    middle.SetBorderSize(2)
    middle.SetTickx(1)
    middle.SetTicky(1)
    middle.SetLeftMargin(0.14)
    middle.SetRightMargin(0.055)
    middle.SetTopMargin(0.045)
    middle.SetBottomMargin(0.4)
    middle.SetFrameBorderMode(0)
    
    bot = c1.cd(3)
    bot.SetPad(0.0,0.0,1.0,0.25)
    bot.SetFillColor(0)
    bot.SetBorderMode(0)
    bot.SetBorderSize(2)
    bot.SetTickx(1)
    bot.SetTicky(1)
    bot.SetLeftMargin(0.14)
    bot.SetRightMargin(0.055)
    bot.SetTopMargin(0)
    bot.SetBottomMargin(0.4)
    bot.SetFrameBorderMode(0)
    
    gluon_hist_array[0].SetTitle("")
    gluon_hist_array[0].GetXaxis().SetTitle(var)
    gluon_hist_array[0].GetYaxis().SetTitle("Normalized to unity")
    gluon_hist_array[0].GetYaxis().SetNdivisions(505)
    gluon_hist_array[0].GetYaxis().SetRangeUser(-0.01,gluon.GetMaximum()*1.5)
    
    for j in range(2):
            gluon_hist_array[j].SetMarkerColor(3*j+1)
            gluon_hist_array[j].SetLineColor(3*j+1)
            gluon_hist_array[j].SetMarkerSize(0.5)
            gluon_hist_array[j].SetLineStyle(1)
  
    for j in range(2):
            gluon_hist_array_truth[j].SetMarkerColor(3*j+1)
            gluon_hist_array_truth[j].SetLineColor(3*j+1)
            gluon_hist_array_truth[j].SetMarkerSize(0.5)
            gluon_hist_array_truth[j].SetLineStyle(2)
      
    for j in range(2):
            gluon_hist_array_ratio[j].SetMarkerColor(3*j+1)
            gluon_hist_array_ratio[j].SetLineColor(3*j+1)
            gluon_hist_array_ratio[j].SetMarkerSize(0.5)
            gluon_hist_array_ratio[j].SetLineStyle(2)
   
            
    gluon_hist_array_ratio[0].SetTitle("")
    gluon_hist_array_ratio[0].GetYaxis().SetRangeUser(0.7,1.3)
    gluon_hist_array_ratio[0].GetXaxis().SetTitleOffset(1)
    gluon_hist_array_ratio[0].GetXaxis().SetTitleSize(0.11)
    gluon_hist_array_ratio[0].GetXaxis().SetLabelSize(0.1)
    gluon_hist_array_ratio[0].GetXaxis().SetLabelOffset(0.03)
    gluon_hist_array_ratio[0].GetYaxis().SetTitleSize(0.1)
    gluon_hist_array_ratio[0].GetYaxis().SetTitleOffset(0.5)
    gluon_hist_array_ratio[0].GetYaxis().SetLabelOffset(0.01)
    gluon_hist_array_ratio[0].GetYaxis().SetLabelSize(0.1)

    fg_gluon_mc_ratio.GetYaxis().SetTitle("\eta_{1}/\eta_{2}")
    fg_gluon_mc_ratio.GetYaxis().SetRangeUser(0.7,1.3)
    fg_gluon_mc_ratio.GetXaxis().SetTitleOffset(1)
    fg_gluon_mc_ratio.GetXaxis().SetTitleSize(0.11)
    fg_gluon_mc_ratio.GetXaxis().SetLabelSize(0.1)
    fg_gluon_mc_ratio.GetXaxis().SetLabelOffset(0.03)
    fg_gluon_mc_ratio.GetYaxis().SetTitleSize(0.1)
    fg_gluon_mc_ratio.GetYaxis().SetTitleOffset(0.5)
    fg_gluon_mc_ratio.GetYaxis().SetLabelOffset(0.01) 
    fg_gluon_mc_ratio.GetYaxis().SetLabelSize(0.1)

    fg_gluon_mc_ratio.SetMarkerSize(0.5)
    fg_gluon_mc_ratio.SetLineStyle(1)  
    fg_gluon_mc_ratio.SetLineColor(1)

    fg_gluon_truth_ratio.SetMarkerSize(0.5)
    fg_gluon_truth_ratio.SetLineStyle(2)  
    fg_gluon_truth_ratio.SetLineColor(1)
    
    top.cd()
    gluon_hist_array[0].Draw("HIST e") #extracted mc
    gluon_hist_array[1].Draw("HIST e same") 
    #gluon_hist_array[2].Draw("HIST same") 
    gluon_hist_array_truth[0].SetMarkerStyle(20)
    gluon_hist_array_truth[1].SetMarkerStyle(20)

    gluon_hist_array_truth[0].Draw("HIST same") #extracted mc
    gluon_hist_array_truth[1].Draw("HIST same") 
    #gluon_hist_array_truth[2].Draw("p same") 
    #higher_gluon.Draw("HIST same") #truth mc
    if(inputvar == "ntrk"):
        line = TLine(0.,1,60,1)
        gluon_hist_array_ratio[0].GetXaxis().SetTitle("N_{track}")
        leg = TLegend(0.50,0.7,0.9,0.90)
        leg.SetTextSize(0.018)
    if(inputvar == "bdt"):
        line = TLine(-0.8,1,0.7,1)
        gluon_hist_array_ratio[0].GetXaxis().SetTitle("BDT score")
        leg = TLegend(0.6,0.7,0.9,0.9) ##0.6,0.5,0.9,0.7
        leg = TLegend(0.50,0.7,0.9,0.90)
        leg.SetTextSize(0.018)
    leg.SetTextFont(42)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(1)
    leg.AddEntry(gluon_hist_array[0],"extracted gluon(mc) \eta_{1}: "+eta_start_central_1+"<|\eta_{C}|<"+eta_end_central_1+","+eta_start_forward_1+"<|\eta_{C}|<"+eta_end_forward_1+"","l")
    leg.AddEntry(gluon_hist_array[1],"extracted gluon(mc) \eta_{2}: "+eta_start_central_2+"<|\eta_{C}|<"+eta_end_central_2+","+eta_start_forward_2+"<|\eta_{C}|<"+eta_end_forward_2+"","l")
    leg.AddEntry(gluon_hist_array_truth[0],"gluon(mc) \eta_{1}:"+eta_start_central_1+"<|\eta_{C}|<"+eta_end_central_1+","+eta_start_forward_1+"<|\eta_{C}|<"+eta_end_forward_1+"","l")
    leg.AddEntry(gluon_hist_array_truth[1],"gluon(mc) \eta_{2}:"+eta_start_central_2+"<|\eta_{C}|<"+eta_end_central_2+","+eta_start_forward_2+"<|\eta_{C}|<"+eta_end_forward_1+"","l")
    leg.Draw()
    myText(0.18,0.84,"#it{#bf{#scale[1.5]{#bf{ATLAS} Internal}}}")
    myText(0.18,0.80,"#bf{#scale[1.2]{#sqrt{s} = 13 TeV}}")
    myText(0.18,0.75,"#bf{#scale[1.2]{pT range: "+str(min)+" - "+str(max)+" GeV}}")
    
    middle.cd()
    fg_gluon_mc_ratio.Draw("HIST e") #extracted mc
    fg_gluon_truth_ratio.Draw("HIST e same")   
    line.Draw("same")
    bot.cd()
    #gluon_ratio.Draw("HIST")
    
    gluon_hist_array_ratio[0].Draw("HIST e") #extracted mc
    gluon_hist_array_ratio[1].Draw("HIST e same") 
    #gluon_hist_array_ratio[2].Draw("HIST same") 
    line.Draw("same")
    c1.Print("./plots_"+var+"/"+"gluon_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+"_compare_"+eta_bin[0]+"_"+eta_bin[1]+".pdf")

