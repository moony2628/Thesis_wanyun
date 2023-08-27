from ROOT import *
import numpy as np
import uncertainties as unc # propagate uncertainties
from uncertainties import unumpy as unp # array operations for type ufloat

var = "ntrk"  #change the var name according to the inputvar you want to read
mc = "pythia"   #The mc sample in the output name
inputvar = var  #by setting it as bdt (or ntrk,width,c1..), it will read the corresponding histogram, but remember to change the TLine range according to X-axis of different variable, one can check it by browsing the histograms in root file.

#decide if we want to do the reweighting process
if var == "ntrk":
    doreweight = "Quark"  
if var == "bdt":
    doreweight = 0

ntrackall = TFile("../newroot/dijet_pythia_nominal.root")
ntrackall3  = TFile("../newroot/dijet_data_prescale.root")


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

bin = [0,50,100,150,200,300,400,500,600,800,1000,1200,1500,2000]
for k in range(7,13):   #for only dijet event, start from jet pT>500 GeV
#for i in range(13):	#for gamma+jet combined with dijet event, start from jet pT>0 GeV
        min = bin[k]
        max = bin[k+1]
        higher_quark2 = ntrackall.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar)
        higher_gluon2 = ntrackall.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar)
        higher_data2 = ntrackall3.Get(str(min)+"_LeadingJet_Forward_Data_"+inputvar)
        lower_quark2 = ntrackall.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar)
        lower_gluon2 = ntrackall.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar)
        lower_data2 = ntrackall3.Get(str(min)+"_LeadingJet_Central_Data_"+inputvar)

        higher_quark = ntrackall.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar)
        higher_gluon = ntrackall.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar)

        lower_quark = ntrackall.Get(str(min)+"_SubJet_Central_Quark_"+inputvar)
        lower_gluon = ntrackall.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar)

        higher_data = ntrackall3.Get(str(min)+"_SubJet_Forward_Data_"+inputvar)
        lower_data = ntrackall3.Get(str(min)+"_SubJet_Central_Data_"+inputvar)
        from ROOT import *    
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
        
        higher_quark2_err = ntrackall.Get(str(min)+"_LeadingJet_Forward_Quark_"+inputvar + "_err")
        higher_gluon2_err = ntrackall.Get(str(min)+"_LeadingJet_Forward_Gluon_"+inputvar + "_err")
        lower_quark2_err = ntrackall.Get(str(min)+"_LeadingJet_Central_Quark_"+inputvar + "_err")
        lower_gluon2_err = ntrackall.Get(str(min)+"_LeadingJet_Central_Gluon_"+inputvar + "_err")

        higher_quark_err = ntrackall.Get(str(min)+"_SubJet_Forward_Quark_"+inputvar + "_err")
        higher_gluon_err = ntrackall.Get(str(min)+"_SubJet_Forward_Gluon_"+inputvar + "_err")

        lower_quark_err = ntrackall.Get(str(min)+"_SubJet_Central_Quark_"+inputvar + "_err")
        lower_gluon_err = ntrackall.Get(str(min)+"_SubJet_Central_Gluon_"+inputvar + "_err")



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

        quark = higher_quark.Clone("")
        gluon = higher_gluon.Clone("")
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
        ex_quark_ratio.GetYaxis().SetTitle("Data/MC")
        ex_gluon_ratio.GetYaxis().SetTitle("Data/MC")

        ex_quark_ratio.Divide(quark)
        ex_gluon_ratio.Divide(gluon)
        quark_ratio.Divide(higher_quark)
        gluon_ratio.Divide(higher_gluon)
        

        gStyle.SetOptStat(0)
        ######################## for ratio plot
        c.Divide(2,1)

        top = c.cd(1)
        top.SetPad(0.0,0.0,1.0,1.0)
        top.SetFillColor(0)
        top.SetBorderMode(0)
        top.SetBorderSize(2)
        top.SetTickx(1)
        top.SetTicky(1)
        top.SetLeftMargin(0.14)
        top.SetRightMargin(0.055)
        top.SetBottomMargin(0.3)#0.25
        top.SetFrameBorderMode(0)
        #top.SetLogy(1)

        bot = c.cd(2)
        bot.SetPad(0.0,0.0,1.0,0.3)
        bot.SetFillColor(0)
        bot.SetBorderMode(0)
        bot.SetBorderSize(2)
        bot.SetTickx(1)
        bot.SetTicky(1)
        bot.SetLeftMargin(0.14)
        bot.SetRightMargin(0.055)
        bot.SetTopMargin(0.045)
        bot.SetBottomMargin(0.4)
        bot.SetFrameBorderMode(0)

        quark.SetTitle("")
        quark.GetXaxis().SetTitle(var)
        quark.GetYaxis().SetTitle("Normalized to unity")
        quark.GetYaxis().SetNdivisions(505)
        quark.GetYaxis().SetRangeUser(-0.01,quark.GetMaximum()*1.5)

        quark.SetMarkerColor(4)
        quark.SetLineColor(4)
        quark.SetMarkerSize(0.5)
        quark.SetMarkerStyle(20)

        quark_data.SetMarkerColor(1)
        quark_data.SetLineColor(1)
        quark_data.SetMarkerSize(0.5)
        quark_data.SetMarkerStyle(8)

        higher_quark.SetMarkerColor(4)
        higher_quark.SetLineColor(4)
        higher_quark.SetMarkerSize(0.5)
        higher_quark.SetLineStyle(2)

        ex_quark_ratio.SetTitle("")
        ex_quark_ratio.GetYaxis().SetRangeUser(0.7,1.3)
        ex_quark_ratio.GetXaxis().SetTitleOffset(1)
        ex_quark_ratio.GetXaxis().SetTitleSize(0.11)
        ex_quark_ratio.GetXaxis().SetLabelSize(0.1)
        ex_quark_ratio.GetXaxis().SetLabelOffset(0.03)
        ex_quark_ratio.GetYaxis().SetTitleSize(0.1)
        ex_quark_ratio.GetYaxis().SetTitleOffset(0.5)
        ex_quark_ratio.GetYaxis().SetLabelOffset(0.01)

        quark_ratio.SetMarkerColor(4)
        quark_ratio.SetLineColor(4)
        quark_ratio.SetMarkerSize(0.5)
        quark_ratio.SetMarkerStyle(20)

        ex_quark_ratio.SetMarkerColor(4)
        ex_quark_ratio.SetLineColor(4)
        ex_quark_ratio.SetMarkerSize(0.5)
        ex_quark_ratio.SetLineStyle(20)

        top.cd()

        quark.Draw("HIST") #extracted mc
        quark_data.Draw("same") #extracted mc
        #higher_quark.Draw("HIST same") #truth mc

        if(inputvar == "ntrk"):
            line = TLine(0.,1,60,1)
            ex_quark_ratio.GetXaxis().SetTitle("N_{track}")
            leg = TLegend(0.6,0.5,0.9,0.7)
        if(inputvar == "bdt"):
            line = TLine(-0.8,1,0.7,1)
            ex_quark_ratio.GetXaxis().SetTitle("BDT")
            leg = TLegend(0.6,0.65,0.9,0.85) ##0.6,0.5,0.9,0.7

        leg.SetTextFont(42)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetNColumns(1)
        leg.AddEntry(quark_data,"Extracted quark (data)","p")
        leg.AddEntry(quark,"Extracted quark (mc)","l")
        #leg.AddEntry(higher_quark,"Quark (mc)","l")
        leg.Draw()

        myText(0.18,0.84,"#it{#bf{#scale[1.8]{#bf{ATLAS} Internal}}}")
        myText(0.18,0.80,"#bf{#scale[1.5]{#sqrt{s} = 13 TeV}}")
        myText(0.18,0.75,"#bf{#scale[1.5]{pT range: "+str(min)+" - "+str(max)+" GeV}}")

        bot.cd()
        #quark_ratio.Draw("HIST")
        ex_quark_ratio.SetFillColor(4)
        ex_quark_ratio.SetFillColor(4)
        ex_quark_ratio.SetFillStyle(3005)
        ex_quark_ratio.Draw("e2")
        line.Draw("same")
        c.Print("./plots_"+var+"/quark_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+".pdf")

        gluon.SetTitle("")
        gluon.GetXaxis().SetTitle(var)
        gluon.GetYaxis().SetTitle("Normalized to unity")
        gluon.GetYaxis().SetNdivisions(505)
        gluon.GetYaxis().SetRangeUser(-0.01,gluon.GetMaximum()*1.5)

        gluon.SetMarkerColor(2)
        gluon.SetLineColor(2)
        gluon.SetMarkerSize(0.5)
        gluon.SetMarkerStyle(20)

        gluon_data.SetMarkerColor(1)
        gluon_data.SetLineColor(1)
        gluon_data.SetMarkerSize(0.5)
        gluon_data.SetMarkerStyle(8)

        higher_gluon.SetMarkerColor(2)
        higher_gluon.SetLineColor(2)
        higher_gluon.SetMarkerSize(0.5)
        higher_gluon.SetLineStyle(2)

        ex_gluon_ratio.SetTitle("")
        ex_gluon_ratio.GetYaxis().SetRangeUser(0.7,1.3)
        ex_gluon_ratio.GetXaxis().SetTitleOffset(1)
        ex_gluon_ratio.GetXaxis().SetTitleSize(0.11)
        ex_gluon_ratio.GetXaxis().SetLabelSize(0.1)
        ex_gluon_ratio.GetXaxis().SetLabelOffset(0.03)
        ex_gluon_ratio.GetYaxis().SetTitleSize(0.1)
        ex_gluon_ratio.GetYaxis().SetTitleOffset(0.5)
        ex_gluon_ratio.GetYaxis().SetLabelOffset(0.01)


        ex_gluon_ratio.SetMarkerColor(2)
        ex_gluon_ratio.SetLineColor(2)
        ex_gluon_ratio.SetMarkerSize(0.5)
        ex_gluon_ratio.SetLineStyle(20)

        top.cd()
        gluon.Draw("HIST")
        gluon_data.Draw("same")
        #higher_gluon.Draw("HIST same")
        
        if(inputvar == "ntrk"):
            line = TLine(0.,1,60,1)
            ex_gluon_ratio.GetXaxis().SetTitle("N_{track}")
            leg = TLegend(0.6,0.5,0.9,0.7)
        if(inputvar == "bdt"):
            line = TLine(-0.8,1,0.7,1)
            ex_gluon_ratio.GetXaxis().SetTitle("BDT")
            leg = TLegend(0.2,0.5,0.5,0.7) ##0.6,0.5,0.9,0.7

        leg.SetTextFont(42)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetNColumns(1)
        leg.AddEntry(gluon_data,"Extracted gluon (data)","p")
        leg.AddEntry(gluon,"Extracted gluon (mc)","l")
        #leg.AddEntry(higher_gluon,"Gluon (mc)","l")
        leg.Draw()

        myText(0.18,0.84,"#it{#bf{#scale[1.8]{#bf{ATLAS} Internal}}}")
        myText(0.18,0.80,"#bf{#scale[1.5]{#sqrt{s} = 13 TeV}}")
        myText(0.18,0.75,"#bf{#scale[1.5]{pT range: "+str(min)+" - "+str(max)+" GeV}}")
        bot = c.cd(2)
        bot.SetPad(0.0,0.0,1.0,0.3)
        bot.SetFillColor(0)
        bot.SetBorderMode(0)
        bot.SetBorderSize(2)
        bot.SetTickx(1)
        bot.SetTicky(1)
        bot.SetLeftMargin(0.14)
        bot.SetRightMargin(0.055)
        bot.SetTopMargin(0.045)
        bot.SetBottomMargin(0.4)
        bot.SetFrameBorderMode(0)
        bot.cd()
        #gluon_ratio.Draw("HIST")
        ex_gluon_ratio.SetMarkerSize(2)
        ex_gluon_ratio.SetFillColor(2)
        ex_gluon_ratio.SetFillStyle(3005)
        ex_gluon_ratio.Draw("e2")
        line.Draw("same")
        c.Print("./plots_"+var+"/gluon_"+str(min)+"_"+str(doreweight)+"_"+mc+"_"+var+".pdf")
