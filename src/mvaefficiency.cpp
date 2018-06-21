#include <iostream>
#include <iomanip>
#include "mvaefficiency.h"
#include "mva_methods.h"
#include "./root_mva/tmvaglob.C"

using namespace std;

MvaEfficiency::MvaEfficiency(int sigN, int bgdN, string *plotdir)
{
   methodName = "";
   methodTitle = "";
   sig = 0;
   bgd = 0;
   origSigE = 0;
   origBgdE = 0;
   sigE = 0;
   bgdE = 0;
   purS = 0;
   sSig = 0;
   effpurS = 0;
   effcanvas = 0;
   line1 = 0;
   line2 = 0;
   line3 = 0;
   line4 = 0;
   plotXmin = -0.1;
   plotXmax = 1.1;
   rightAxis = 0;
   maxSignificance = 0;
   maxSignificanceErr = 0;
   fNSignal = sigN;
   fNBackground = bgdN;
   fFormula = "S/sqrt(S+B)";
   maxLenTitle = 0;
   sigpurCutBest = -1;
   sigbgdCutBest = -1;
   found[0] = false;
   found[1] = false;
   optimalCut = -1;
   sigbgdCut = -1;
   sigpurCut = -1;
   maxbin = -1;
   plotlocation = *plotdir;
}

MvaEfficiency::~MvaEfficiency()
{
   delete c;
   delete legend1;
   delete legend2;
   delete effline;

   delete rightAxis;
   delete purS;
   delete sSig;
   delete effpurS;
}

void MvaEfficiency::RunMvaEfficiency(TString *inname)
{
   TMVAGlob::Initialize(kTRUE);
  
   TFile *globFile = TMVAGlob::OpenFile(*inname);
  
   ReadHistogram(globFile);
   SetFormula(fFormula);
   UpdateSignificanceHist();
   DrawHistogram();
}

void MvaEfficiency::SetPlotLimits()
{
   cout << endl << "# SetPlotLimits()" << endl;

   plotXmin = GetMethodMin(methodTitle);
   plotXmax = GetMethodMax(methodTitle);

/*   // setup the ploting limits, depending on the used mva method
   if(methodTitle.Contains("CutsD"))
   {
      plotXmin = 0.;
      plotXmax = 1.;
   }
   else if(methodTitle.Contains("Cuts"))
   {
      plotXmin = 0.;
      plotXmax = 1.;
   }
   else if(methodTitle.Contains("LikelihoodPCA"))
   {
      plotXmin = 0.;
      plotXmax = 1.;
   }
   else if(methodTitle.Contains("Likelihood"))
   {
      plotXmin = 0.;
      plotXmax = 1.;
   }
   else if(methodTitle.Contains("PDERS"))
   {
      plotXmin = 0.;
      plotXmax = 1.;
   }
   else if(methodTitle.Contains("KNN"))
   {
      plotXmin = 0.;
      plotXmax = 1.;
   }
   else if(methodTitle.Contains("LD"))
   {
      plotXmin = -1.;
      plotXmax = 1.;
   }
   else if(methodTitle.Contains("Fisher"))
   {
      plotXmin = -1.6;
      plotXmax = 1.6;
   }
   else if(methodTitle.Contains("FDA_GA"))
   {
      plotXmin = 0.;
      plotXmax = 1.;
   }
   else if(methodTitle.Contains("MLPBNN"))
   {
      plotXmin = -0.1;
      plotXmax = 1.1;
   }
   else if(methodTitle.Contains("SVM"))
   {
      plotXmin = 0.;
      plotXmax = 1.;
   }
   else if(methodTitle.Contains("BDT"))
   {
      plotXmin = -1.;
      plotXmax = 1.;
   }
   else if(methodTitle.Contains("RuleFit"))
   {
      plotXmin = -1.6;
      plotXmax = 1.6;
   }
*/
}

void MvaEfficiency::ReadHistogram(TFile *file)
{
   cout << endl << "# ReadHistogram()" << endl;

   // search for the right histograms in full list of keys
   TIter next(file->GetListOfKeys());
   TKey *key(0);
   while( key = (TKey*)next() )
   {
      if(!TString(key->GetName()).BeginsWith("Method_")) continue;
      if(!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) continue;

      cout << "--- Found directory: " << ((TDirectory*)key->ReadObj())->GetName() << endl;

      TDirectory *mDir = (TDirectory*)key->ReadObj();

      TIter keyIt(mDir->GetListOfKeys());
      TKey *titkey;
      while( titkey = (TKey*)keyIt() )
      {
         if(!gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory")) continue;
	 TDirectory *titDir = (TDirectory*)titkey->ReadObj();

	 TMVAGlob::GetMethodName(methodName, key);
	 cout << "--- Got method name: " << methodName << endl;
         TMVAGlob::GetMethodTitle(methodTitle, titDir);        
         if (methodTitle.Length() > maxLenTitle)
            maxLenTitle = methodTitle.Length();
	 TString hname = "MVA_" + methodTitle;
         cout << "--- Classifier: " << methodTitle << endl;

	 SetPlotLimits();

	 sig = dynamic_cast<TH1*>(titDir->Get(hname + "_S"));
	 bgd = dynamic_cast<TH1*>(titDir->Get(hname + "_B"));
	 origSigE = dynamic_cast<TH1*>(titDir->Get(hname + "_effS"));
	 origBgdE = dynamic_cast<TH1*>(titDir->Get(hname + "_effB"));
         if(origSigE == 0 || origBgdE == 0) continue;

	 SetResultHist();
      }
   }
}

void MvaEfficiency::SetFormula(const TString& f)
{
   cout << endl << "# SetFormula()" << endl;
   fFormula = f;
}

void MvaEfficiency::UpdateSignificanceHist()
{
   cout << endl << "# UpdateSignificanceHist()" << endl;

   TFormula f("sigf", GetFormula());
   TString cname = "Classifier";
   if(cname.Length() > maxLenTitle)
      maxLenTitle = cname.Length();
   TString str = Form("%*s   (  #signal, #backgr.)  Optimal-cut  %s      NSig      NBkg   EffSig   EffBkg   SigPur", maxLenTitle, cname.Data(), fFormula.Data());
   cout << "--- " << setfill('=') << setw(str.Length()) << "" << setfill(' ') << endl;
   cout << "--- " << str << endl;
   cout << "--- " << setfill('-') << setw(str.Length()) << "" << setfill(' ') << endl;

   double maxSig = -1;
   double maxSigErr = -1;

   for(int i = 1; i < origSigE->GetNbinsX(); i++)
   {
      float eS = origSigE->GetBinContent(i);
      float S = eS*fNSignal;
      float B = origBgdE->GetBinContent(i)*fNBackground;
      purS->SetBinContent(i, (S+B == 0) ? 0 : S/(S+B));
//cout << "# UpdateSignificance    #: " << "Bin " << i << ", SigE = " << eS << ", BgdE = " << origBgdE->GetBinContent(i) << ", Signal = " << S << ", Background = " << B << ", S/(S+B) = " << S/(S+B) << endl;

      double dsig = f.Eval(S, B);
      if(dsig > maxSig)
      {
         maxSig = dsig;
	 if(fFormula == "S/sqrt(B)")
	 {
            maxSigErr = dsig*sqrt(1./S + 1./(2.*B));
	 }
      }
      sSig->SetBinContent(i, dsig);
      effpurS->SetBinContent(i, eS*purS->GetBinContent(i));
   }
   maxSignificance = sSig->GetMaximum();
   maxSignificanceErr = (maxSigErr > 0) ? maxSigErr : 0;
   sSig->Scale(1/maxSignificance);

   PrintResults();
   cout << "--- " << setfill('-') << setw(str.Length()) << "" << setfill(' ') << endl << endl;
}

void MvaEfficiency::DrawHistogram()
{
   cout << endl << "# DrawHistogram()" << endl;

   // define Canvas layout here!
   const int width = 1000;   // size of canvas
   int signifColor = TColor::GetColor("#00aa00");

   // create new canvas
   c = new TCanvas("canvas1", Form("Cut efficiencies for %s classifier", methodTitle.Data()), 200, 0, width, (int)(width*0.78)); 
   effcanvas = c;
   
   // draw grid
   c->SetGrid(1);
   c->SetTickx(0);
   c->SetTicky(0);
   
   TStyle *TMVAStyle = gROOT->GetStyle("Plain"); // our style is based on Plain
   TMVAStyle->SetLineStyleString(5, "[32 22]");
   TMVAStyle->SetLineStyleString(6, "[12 22]");
      
   c->SetTopMargin(.2);
   
   // and the signal purity and quality
   effpurS->SetTitle("Cut efficiencies and optimal cut value");
   if (methodTitle.Contains("Cuts"))
      effpurS->GetXaxis()->SetTitle("Signal Efficiency");
   else
      effpurS->GetXaxis()->SetTitle(TString("Cut value applied on ") + methodTitle + " output");

   effpurS->GetYaxis()->SetTitle("Efficiency (Purity)");
   TMVAGlob::SetFrameStyle(effpurS);
   
   c->SetTicks(0,0);
   c->SetRightMargin(2.0);
   
   effpurS->SetMaximum(1.1);
   effpurS->Draw("histl");
   
   purS->Draw("samehistl");      
   
   // overlay signal and background histograms
   sigE->Draw("samehistl");
   bgdE->Draw("samehistl");
   
   sSig->SetLineColor(signifColor);
   sSig->Draw("samehistl");
   
   // redraw axes
   effpurS->Draw( "sameaxis" );
   
   // Draw legend               
   legend1= new TLegend(c->GetLeftMargin(), 1 - c->GetTopMargin(), c->GetLeftMargin() + 0.4, 1 - c->GetTopMargin() + 0.12);
   legend1->SetFillStyle(1);
   legend1->AddEntry(sigE, "Signal efficiency", "L");
   legend1->AddEntry(bgdE, "Background efficiency", "L");
   legend1->Draw("same");
   legend1->SetBorderSize(1);
   legend1->SetMargin(0.3);
   
   legend2= new TLegend(c->GetLeftMargin() + 0.4, 1 - c->GetTopMargin(), 1 - c->GetRightMargin(), 1 - c->GetTopMargin() + 0.12);
   legend2->SetFillStyle(1);
   legend2->AddEntry(purS, "Signal purity", "L");
   legend2->AddEntry(effpurS, "Signal efficiency*purity", "L");
   legend2->AddEntry(sSig, GetLatexFormula().Data(), "L");
   legend2->Draw("same");
   legend2->SetBorderSize(1);
   legend2->SetMargin(0.3);
      
   // line to indicate maximum efficiency
   effline = new TLine(sSig->GetXaxis()->GetXmin(), 1, sSig->GetXaxis()->GetXmax(), 1);
   effline->SetLineWidth(1);
   effline->SetLineColor(1);
   effline->Draw();
   
   // print comments
   TLatex tl;
   tl.SetNDC();
//   tl.SetTextSize(0.033);
   tl.SetTextSize(0.03);
   maxbin = sSig->GetMaximumBin();
   
   if(maxSignificanceErr > 0)
   {
      line1 = tl.DrawLatex(0.15, 0.23, Form("For %1.0f signal and %1.0f background", fNSignal, fNBackground));
      tl.DrawLatex(0.15, 0.19, "events the maximum " + GetLatexFormula() + " is");
      line2 = tl.DrawLatex(0.15, 0.15, Form("%5.2f +- %4.2f when cutting at %5.4f", maxSignificance, maxSignificanceErr, sSig->GetXaxis()->GetBinCenter(maxbin)));
   }
   else
   {
      line1 = tl.DrawLatex(0.15, 0.31, Form("For %1.0f signal and %1.0f background", fNSignal, fNBackground));
      tl.DrawLatex(0.15, 0.27, "events the maximum " + GetLatexFormula() + " is");
      line2 = tl.DrawLatex(0.15, 0.23, Form("%4.2f when cutting at %5.4f (optimal)", GetHistValue(0, 3), optimalCut));
      line3 = tl.DrawLatex(0.15, 0.19, Form("%4.2f when cutting at %5.4f (equal sig and bgd)", GetHistValue(1, 3), sigbgdCut));
      line4 = tl.DrawLatex(0.15, 0.15, Form("%4.2f when cutting at %5.4f (equal sig and sig purity)", GetHistValue(2, 3), sigpurCut));
   }
   
   // add comment for Method cuts
   if(methodTitle.Contains("Cuts"))
   {
      tl.DrawLatex(0.13, 0.77, "Method Cuts provides a bundle of cut selections, each tuned to a");
      tl.DrawLatex(0.13, 0.74, "different signal efficiency. Shown is the purity for each cut selection.");
   }
   // save canvas to file
   c->Update();

   effpurS->GetXaxis()->SetRangeUser(plotXmin, plotXmax);
   
   // Draw second axes
   rightAxis = new TGaxis(plotXmax, c->GetUymin(), plotXmax, c->GetUymax(), 0, 1.1*maxSignificance, 510, "+L");
   rightAxis->SetLineColor(signifColor);
   rightAxis->SetLabelColor(signifColor);
   rightAxis->SetTitleColor(signifColor);
   
   rightAxis->SetTitleSize(sSig->GetXaxis()->GetTitleSize());
   rightAxis->SetTitle("Significance");
   rightAxis->Draw();
   
   c->Update();
   
   const bool Save_Images = kTRUE;
   if(Save_Images)
   {
      plotlocation = plotlocation + string(Form("/mvaeffs_%s", methodTitle.Data()));
      TMVAGlob::imgconv(c, plotlocation.c_str());
   }
}

void MvaEfficiency::SetResultHist()
{
   cout << endl << "# SetResultHist()" << endl;

   TString pname    = "purS_"         + methodTitle;
   TString epname   = "effpurS_"      + methodTitle;
   TString ssigname = "significance_" + methodTitle;
   
   sigE = (TH1*)origSigE->Clone("sigEffi");
   bgdE = (TH1*)origBgdE->Clone("bgdEffi");
   
   Int_t nbins = sigE->GetNbinsX();
   Double_t low = sigE->GetBinLowEdge(1);
   Double_t high = sigE->GetBinLowEdge(nbins+1);
   purS    = new TH1F(pname, pname, nbins, low, high);
   sSig    = new TH1F(ssigname, ssigname, nbins, low, high);
   effpurS = new TH1F(epname, epname, nbins, low, high);        
   
   // chop off useless stuff
   sigE->SetTitle( Form("Cut efficiencies for %s classifier", methodTitle.Data()) );
      
   // set the histogram style
   TMVAGlob::SetSignalAndBackgroundStyle( sigE, bgdE );
   TMVAGlob::SetSignalAndBackgroundStyle( purS, bgdE );
   TMVAGlob::SetSignalAndBackgroundStyle( effpurS, bgdE );
   sigE->SetFillStyle(0);
   bgdE->SetFillStyle(0);
   sSig->SetFillStyle(0);
   sigE->SetLineWidth(3);
   bgdE->SetLineWidth(3);
   sSig->SetLineWidth(3);
   
   // the purity and quality
   purS->SetFillStyle(0);
   purS->SetLineWidth(2);
   purS->SetLineStyle(5);
   effpurS->SetFillStyle(0);
   effpurS->SetLineWidth(2);
   effpurS->SetLineStyle(6);
}

TString MvaEfficiency::GetFormula()
{
   cout << endl << "# GetFormula()" << endl;
   TString f = fFormula;
   f.ReplaceAll("S", "x");
   f.ReplaceAll("B", "y");
   return f;
}

TString MvaEfficiency::GetLatexFormula()
{
   cout << endl << "# GetLatexFormula()" << endl;
   TString f = fFormula;
   f.ReplaceAll("(","{");
   f.ReplaceAll(")","}");
   f.ReplaceAll("sqrt","#sqrt");
   return f;
}

void MvaEfficiency::PrintResults()
{
   maxbin = sSig->GetMaximumBin();
   if(line1 != 0)
      line1->SetText(0.15, 0.23, Form("For %1.0f signal and %1.0f background", fNSignal, fNBackground));
   
   if(line2 != 0)
   {
      if(maxSignificanceErr > 0)
         line2->SetText(0.15, 0.15, Form("%3.2g +- %3.2g when cutting at %3.4g", maxSignificance, maxSignificanceErr, sSig->GetXaxis()->GetBinCenter(maxbin)));
      else
      {
         line2->SetText(0.15, 0.15, Form("%3.4f when cutting at %3.4f", maxSignificance, sSig->GetXaxis()->GetBinCenter(maxbin)));
      }
   }

   if(maxSignificanceErr <= 0)
   {
      TString opt = Form("%%%is:  (%%9.8g,%%9.8g)    %%9.4f   %%10.6g  %%8.7g  %%8.7g %%8.4g %%8.4g %%8.4g", maxLenTitle);
      cout << "--- " << Form(opt.Data(), methodTitle.Data(), fNSignal, fNBackground, sSig->GetXaxis()->GetBinCenter(maxbin), maxSignificance, origSigE->GetBinContent(maxbin)*fNSignal, origBgdE->GetBinContent(maxbin)*fNBackground, origSigE->GetBinContent(maxbin), origBgdE->GetBinContent(maxbin), purS->GetBinContent(maxbin)) << endl;
   }
   else
   {
      TString opt = Form("%%%is:  (%%9.8g,%%9.8g)    %%9.4f   (%%8.3g  +-%%6.3g)  %%8.7g  %%8.7g %%8.4g %%8.4g %%8.4g", maxLenTitle);
      cout << "--- " << Form(opt.Data(), methodTitle.Data(), fNSignal, fNBackground, sSig->GetXaxis()->GetBinCenter(maxbin), maxSignificance, maxSignificanceErr, origSigE->GetBinContent(maxbin)*fNSignal, origBgdE->GetBinContent(maxbin)*fNBackground, origSigE->GetBinContent(maxbin), origBgdE->GetBinContent(maxbin), purS->GetBinContent(maxbin)) << endl;
   }

/*   cout << "Number of bins in each histogram:" << endl;
   cout << "- sig = " << sig->GetNbinsX() << endl;
   cout << "- bgd = " << bgd->GetNbinsX() << endl;
   cout << "- origSigE = " << origSigE->GetNbinsX() << endl;
   cout << "- origBgdE = " << origBgdE->GetNbinsX() << endl;
   cout << "- sigE = " << sigE->GetNbinsX() << endl;
   cout << "- bgdE = " << bgdE->GetNbinsX() << endl;
   cout << "- purS = " << purS->GetNbinsX() << endl;
   cout << "- sSig = " << sSig->GetNbinsX() << endl;
   cout << "- effpurS = " << effpurS->GetNbinsX() << endl;*/

   for(int i = 1; i <= sigE->GetNbinsX(); i++)
   {
//      cout << i << "\t" << sigE->GetXaxis()->GetBinCenter(i) << "\t" << sigE->GetBinContent(i) << "\t" << bgdE->GetBinContent(i) << "\t" << purS->GetBinContent(i) << endl;

      if( (purS->GetBinContent(i) >= sigE->GetBinContent(i)) && (!found[0]) )
      {
	 sigpurCutBest = i;
	 found[0] = true;
      }

      if( (1.-bgdE->GetBinContent(i) >= sigE->GetBinContent(i)) && (!found[1]) )
      {
	 sigbgdCutBest = i;
	 found[1] = true;
      }

      if(found[0] && found[1]) break;
   }

   if(purS->GetBinContent(sigpurCutBest) > sigE->GetBinContent(sigpurCutBest))
      sigpurCutBest = sigpurCutBest - 1;

   if(1.-bgdE->GetBinContent(sigbgdCutBest) > sigE->GetBinContent(sigbgdCutBest))
      sigbgdCutBest = sigbgdCutBest - 1;

   cout << endl << "Optimal cut (" << maxbin << ", " << sigE->GetXaxis()->GetBinCenter(maxbin) << "), sig = " << sigE->GetBinContent(maxbin) << ", back = " << 1.-bgdE->GetBinContent(maxbin) << ", purity = " << purS->GetBinContent(maxbin) << endl;
   cout << "Best cut (" << sigpurCutBest << ", " << sigE->GetXaxis()->GetBinCenter(sigpurCutBest) << "), sig = " << sigE->GetBinContent(sigpurCutBest) << ", back = " << 1.-bgdE->GetBinContent(sigpurCutBest) << ", purity = " << purS->GetBinContent(sigpurCutBest) << endl;
   cout << "Equal cut (" << sigbgdCutBest << ", " << sigE->GetXaxis()->GetBinCenter(sigbgdCutBest) << "), sig = " << sigE->GetBinContent(sigbgdCutBest) << ", back = " << 1.-bgdE->GetBinContent(sigbgdCutBest) << ", purity = " << purS->GetBinContent(sigbgdCutBest) << endl;

   optimalCut = sigE->GetXaxis()->GetBinCenter(maxbin);
   sigpurCut = sigE->GetXaxis()->GetBinCenter(sigpurCutBest);
   sigbgdCut = sigE->GetXaxis()->GetBinCenter(sigbgdCutBest);

   cout << endl;
}

double MvaEfficiency::GetHistValue(int type, int sigbgdpur)
{
   if(type == 0)		// optimal cut
   {
      if(sigbgdpur == 0)	// signal
         return sigE->GetBinContent(maxbin);
      else if(sigbgdpur == 1)	// background
         return 1.-bgdE->GetBinContent(maxbin);
      else if(sigbgdpur == 2)	// purity
         return purS->GetBinContent(maxbin);
      else if(sigbgdpur == 3)	// signal-to-noise ratio
         return (sSig->GetBinContent(maxbin))*maxSignificance;
      else
         return -1.;
   }
   else if(type == 1)		// signal/background best cut
   {
      if(sigbgdpur == 0)	// signal
         return sigE->GetBinContent(sigbgdCutBest);
      else if(sigbgdpur == 1)	// background
         return 1.-bgdE->GetBinContent(sigbgdCutBest);
      else if(sigbgdpur == 2)	// purity
         return purS->GetBinContent(sigbgdCutBest);
      else if(sigbgdpur == 3)	// signal-to-noise ratio
         return (sSig->GetBinContent(sigbgdCutBest))*maxSignificance;
      else
         return -1.;
   }
   else if(type == 2)		// signal/purity best cut
   {
      if(sigbgdpur == 0)	// signal
         return sigE->GetBinContent(sigpurCutBest);
      else if(sigbgdpur == 1)	// background
         return 1.-bgdE->GetBinContent(sigpurCutBest);
      else if(sigbgdpur == 2)	// purity
         return purS->GetBinContent(sigpurCutBest);
      else if(sigbgdpur == 3)	// signal-to-noise ratio
         return (sSig->GetBinContent(sigpurCutBest))*maxSignificance;
      else
         return -1.;
   }
   else
      return -1.;
}
