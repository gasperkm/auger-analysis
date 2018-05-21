void print_fractions(char *filename)
{
   string stemp;
   double dtemp[3];
   vector<double> xbin;
   vector<double> fracSig, fracBack, fracData, fracNormData;

   TCanvas *c1 = new TCanvas("c1","",1200,900);
   c1->SetGrid();
   gStyle->SetEndErrorSize(6);

   TFile *infile = TFile::Open(filename);
   TGraphAsymmErrors *ingr;
   TList *tlist = (TList*)infile->GetListOfKeys();
   int nrkeys = infile->GetNkeys();

   for(int i = 0; i < nrkeys; i++)
   {
      stemp = tlist->At(i)->GetName();
      cout << "Printing out fractions for " << stemp << endl;
      ingr = (TGraphAsymmErrors*)infile->Get(stemp.c_str());
      for(int j = 0; j < ingr->GetN(); j++)
      {
         ingr->GetPoint(j, dtemp[0], dtemp[1]);
	 xbin.push_back(dtemp[0]);
	 xbin.push_back(ingr->GetErrorXlow(j));
	 xbin.push_back(ingr->GetErrorXhigh(j));

	 if(stemp.find("sig_") != string::npos)
	 {
            fracSig.push_back(dtemp[1]);
            fracSig.push_back(ingr->GetErrorYlow(j));
            fracSig.push_back(ingr->GetErrorYhigh(j));
	    printf("  %d: %6.2lf\t%6.6lf\t%6.6lf\t%6.6lf\t%6.6lf\t%6.6lf\n", j, xbin[0], xbin[1], xbin[2], dtemp[1], ingr->GetErrorYlow(j), ingr->GetErrorYhigh(j));
	 }

	 if(stemp.find("back_") != string::npos)
	 {
            fracBack.push_back(100.-dtemp[1]);
            fracBack.push_back(ingr->GetErrorYlow(j));
            fracBack.push_back(ingr->GetErrorYhigh(j));
	    printf("  %d: %6.2lf\t%6.6lf\t%6.6lf\t%6.6lf\t%6.6lf\t%6.6lf\n", j, xbin[0], xbin[1], xbin[2], dtemp[1], ingr->GetErrorYlow(j), ingr->GetErrorYhigh(j));
	 }

	 if(stemp.find("data_") != string::npos)
	 {
            fracData.push_back(dtemp[1]);
            fracData.push_back(ingr->GetErrorYlow(j));
            fracData.push_back(ingr->GetErrorYhigh(j));
	    printf("  %d: %6.2lf\t%6.6lf\t%6.6lf\t%6.6lf\t%6.6lf\t%6.6lf\n", j, xbin[0], xbin[1], xbin[2], dtemp[1], ingr->GetErrorYlow(j), ingr->GetErrorYhigh(j));
	 }

	 if(stemp.find("datanorm_") != string::npos)
	 {
            fracNormData.push_back(dtemp[1]);
            fracNormData.push_back(ingr->GetErrorYlow(j));
            fracNormData.push_back(ingr->GetErrorYhigh(j));
	    printf("  %d: %6.2lf\t%6.6lf\t%6.6lf\t%6.6lf\t%6.6lf\t%6.6lf\n", j, xbin[0], xbin[1], xbin[2], dtemp[1], ingr->GetErrorYlow(j), ingr->GetErrorYhigh(j));
	 }
      }

      cout << endl;
   }

   int nrpoints = fracSig.size()/3;

   cout << "Signal fractions:" << endl;
   for(int i = 0; i < nrpoints; i++)
      cout << "  " << fracSig[3*i] << "\t" << fracSig[3*i+1] << "\t" << fracSig[3*i+2] << endl;
   cout << endl;

   cout << "Background fractions:" << endl;
   for(int i = 0; i < nrpoints; i++)
      cout << "  " << fracBack[3*i] << "\t" << fracBack[3*i+1] << "\t" << fracBack[3*i+2] << endl;
   cout << endl;

   cout << "Data fractions:" << endl;
   for(int i = 0; i < nrpoints; i++)
      cout << "  " << fracData[3*i] << "\t" << fracData[3*i+1] << "\t" << fracData[3*i+2] << endl;
   cout << endl;

   cout << "Normated data fractions:" << endl;
   for(int i = 0; i < nrpoints; i++)
      cout << "  " << fracNormData[3*i] << "\t" << fracNormData[3*i+1] << "\t" << fracNormData[3*i+2] << endl;
   cout << endl;

   // Calculate lnA
/*   cout << "Calculating lnA with: lnA = r_p*ln(1) + r_Fe*ln(56)" << endl;
   for(int i = 0; i < nrpoints; i++)
   {
      dtemp[0] = (fracData[3*i]/100.)*TMath::Log(1) + (1. - fracData[3*i]/100.)*TMath::Log(56);
      cout << "  lnA = " << dtemp[0] << ",\tA = " << TMath::Exp(dtemp[0]) << endl;
   }
   cout << endl;*/
   cout << "Calculating lnA with (no normalization): lnA = r_p*ln(1) + r_Fe*ln(56)" << endl;
   for(int i = 0; i < nrpoints; i++)
   {
      dtemp[0] = (fracData[3*i]/100.)*TMath::Log(1) + (1. - fracData[3*i]/100.)*TMath::Log(56);
      cout << "  lnA = " << dtemp[0] << ",\tA = " << TMath::Exp(dtemp[0]) << endl;
   }
   cout << endl;

   cout << "Calculating lnA with (with normalization): lnA = r_p*ln(1) + r_Fe*ln(56)" << endl;
   for(int i = 0; i < nrpoints; i++)
   {
      dtemp[0] = (fracNormData[3*i]/100.)*TMath::Log(1) + (1. - fracNormData[3*i]/100.)*TMath::Log(56);
      cout << "  lnA = " << dtemp[0] << ",\tA = " << TMath::Exp(dtemp[0]) << endl;
   }
   cout << endl;

   cout << "Calculating lnA with (no normalization): lnA = [ln(56)*(r_psig - r_pdata)/(r_psig - r_pback)]" << endl;
   for(int i = 0; i < nrpoints; i++)
   {
      dtemp[0] = TMath::Log(56)*(fracSig[3*i] - fracData[3*i])/(fracSig[3*i] - fracBack[3*i]);
      cout << "  lnA = " << dtemp[0] << ",\tA = " << TMath::Exp(dtemp[0]) << endl;
   }
   cout << endl;

   cout << "Calculating lnA with (with normalization): lnA = [ln(56)*(r_psig - r_pdatanorm)/(r_psig - r_pback)]" << endl;
   for(int i = 0; i < nrpoints; i++)
   {
      dtemp[0] = TMath::Log(56)*(fracSig[3*i] - fracNormData[3*i])/(fracSig[3*i] - fracBack[3*i]);
      cout << "  lnA = " << dtemp[0] << ",\tA = " << TMath::Exp(dtemp[0]) << endl;
   }
   cout << endl;

   cout << "Calculating lnA with (no normalization): lnA = [ln(56)*(100. - r_pdata)/(100. - 0.)]" << endl;
   for(int i = 0; i < nrpoints; i++)
   {
      dtemp[0] = TMath::Log(56)*(100. - fracData[3*i])/100.;
      cout << "  lnA = " << dtemp[0] << ",\tA = " << TMath::Exp(dtemp[0]) << endl;
   }
   cout << endl;

   cout << "Calculating lnA with (with normalization): lnA = [ln(56)*(100. - r_pdatanorm)/(100. - 0.)]" << endl;
   for(int i = 0; i < nrpoints; i++)
   {
      dtemp[0] = TMath::Log(56)*(100. - fracNormData[3*i])/100.;
      cout << "  lnA = " << dtemp[0] << ",\tA = " << TMath::Exp(dtemp[0]) << endl;
   }
   cout << endl;

   // Plot lnA
   TGraphAsymmErrors *outgr = new TGraphAsymmErrors();
   outgr->Set(nrpoints);

   cout << "Plotting lnA with (no rebalancing): lnA = [ln(56)*(r_psig - r_pdata)/(r_psig - r_pback)]" << endl;
   for(int i = 0; i < nrpoints; i++)
   {
      dtemp[0] = TMath::Log(56)*(fracSig[3*i] - fracData[3*i])/(fracSig[3*i] - fracBack[3*i]);
      dtemp[1] = TMath::Log(56)*((fracSig[3*i]-fracSig[3*i+1]) - (fracData[3*i]-fracData[3*i+1]))/((fracSig[3*i]-fracSig[3*i+1]) - (fracBack[3*i]-fracBack[3*i+1]));
      dtemp[1] = TMath::Abs(dtemp[0] - dtemp[1]);
      dtemp[2] = TMath::Log(56)*((fracSig[3*i]+fracSig[3*i+2]) - (fracData[3*i]+fracData[3*i+2]))/((fracSig[3*i]+fracSig[3*i+2]) - (fracBack[3*i]+fracBack[3*i+2]));
      dtemp[2] = TMath::Abs(dtemp[0] - dtemp[2]);

      outgr->SetPoint(i, xbin[3*i], dtemp[0]);
      outgr->SetPointError(i, xbin[3*i+1], xbin[3*i+2], dtemp[1], dtemp[2]);
   }

   outgr->SetLineColor(2);
   outgr->SetLineWidth(2);
   outgr->GetYaxis()->SetRangeUser(-0.5,4.5);
   outgr->SetTitle(";FD energy [log(E/eV)];#left(lnA#right)");
   outgr->GetXaxis()->CenterTitle();
   outgr->GetYaxis()->CenterTitle();
   outgr->Draw("AP");

   // Draw axis on the right side
   c1->Update();
   TText *t = new TText();
   t->SetTextAlign(12);
   t->SetTextSize(0.035);
   t->SetTextFont(72);
   t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(1),  "p");
   t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(4),  "He");
   t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(14), "N");
   t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(56), "Fe");

   // Draw lines marking the particle types
   TLine *l = new TLine();
   l->SetLineWidth(2);
   l->SetLineStyle(7);
   l->DrawLine(gPad->GetUxmin(), TMath::Log(1), gPad->GetUxmax(), TMath::Log(1));
   l->DrawLine(gPad->GetUxmin(), TMath::Log(4), gPad->GetUxmax(), TMath::Log(4));
   l->DrawLine(gPad->GetUxmin(), TMath::Log(14), gPad->GetUxmax(), TMath::Log(14));
   l->DrawLine(gPad->GetUxmin(), TMath::Log(56), gPad->GetUxmax(), TMath::Log(56));

   c1->SaveAs("lnA_mva.pdf");
   delete outgr;

   // Plot lnA - normalized
   outgr = new TGraphAsymmErrors();
   outgr->Set(nrpoints);

   cout << "Plotting lnA with (with rebalancing): lnA = [ln(56)*(r_psig - r_pdatanorm)/(r_psig - r_pback)]" << endl;
   for(int i = 0; i < nrpoints; i++)
   {
      dtemp[0] = TMath::Log(56)*(100. - fracNormData[3*i])/100.;
      dtemp[1] = TMath::Log(56)*(100. - (fracNormData[3*i]-fracData[3*i+1]))/100.;
      dtemp[1] = TMath::Abs(dtemp[0] - dtemp[1]);
      dtemp[2] = TMath::Log(56)*(100. - (fracNormData[3*i]+fracData[3*i+2]))/100.;
      dtemp[2] = TMath::Abs(dtemp[0] - dtemp[2]);

      outgr->SetPoint(i, xbin[3*i], dtemp[0]);
      outgr->SetPointError(i, xbin[3*i+1], xbin[3*i+2], dtemp[1], dtemp[2]);
   }

   outgr->SetLineColor(2);
   outgr->SetLineWidth(2);
   outgr->GetYaxis()->SetRangeUser(-0.5,4.5);
   outgr->SetTitle(";FD energy [log(E/eV)];#left(lnA#right)");
   outgr->GetXaxis()->CenterTitle();
   outgr->GetYaxis()->CenterTitle();
   outgr->Draw("AP");

   // Draw axis on the right side
   c1->Update();
   t = new TText();
   t->SetTextAlign(12);
   t->SetTextSize(0.035);
   t->SetTextFont(72);
   t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(1),  "p");
   t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(4),  "He");
   t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(14), "N");
   t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(56), "Fe");

   // Draw lines marking the particle types
   l = new TLine();
   l->SetLineWidth(2);
   l->SetLineStyle(7);
   l->DrawLine(gPad->GetUxmin(), TMath::Log(1), gPad->GetUxmax(), TMath::Log(1));
   l->DrawLine(gPad->GetUxmin(), TMath::Log(4), gPad->GetUxmax(), TMath::Log(4));
   l->DrawLine(gPad->GetUxmin(), TMath::Log(14), gPad->GetUxmax(), TMath::Log(14));
   l->DrawLine(gPad->GetUxmin(), TMath::Log(56), gPad->GetUxmax(), TMath::Log(56));

   c1->SaveAs("lnA_mva_normalized.pdf");
}
