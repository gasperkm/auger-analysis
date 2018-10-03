#define _STANDALONE_ 1
#include "workstation.h"
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include "separate_functions.h"
#include "mva_methods.h"
#include "mva_result_read.h"
#include "root_style.h"
#include "primary_type.h"
#if OFFVER == 0
   #include "OfflineIncludeOld.h"
#elif OFFVER == 1
   #include "OfflineIncludeNew.h"
#endif

using namespace std;

class AdstFile
{
private:
   string *stemp;
   int *itemp;
   double *dtemp;

   RecEventFile *fFile;
   RecEvent *fRecEvent;
   DetectorGeometry *fDetGeo;
   SdRecShower *sdrecshw;
   vector<SdRecStation> stationVector;
   vector<FDEvent> acteyes;

   double *wQuantity, *quantitySum, *wQuantitySum, *quantitySumErr;

   // Temporary variables and holders for number of stations and events
   int nrstations, nevents, nreyes;
   bool goodrec;
   bool isHeco;
public:
   AdstFile();
   virtual ~AdstFile();

   void ReadAdstFile(string inname, double *energylimitFull, double *zenithlimitFull, vector<double> *energy, vector<double> *zenith, vector<double> *shwsize, int *allevts);
};

AdstFile::AdstFile()
{
   fRecEvent = new RecEvent();
   fDetGeo = new DetectorGeometry();
   sdrecshw = new SdRecShower();

   stemp = new string[3];
   itemp = new int[4];
   dtemp = new double[4];

   wQuantity = new double[2];
   quantitySum = new double[2];
   wQuantitySum = new double[2];
   quantitySumErr = new double[2];
}

AdstFile::~AdstFile()
{
   delete[] wQuantity;
   delete[] quantitySum;
   delete[] wQuantitySum;
   delete[] quantitySumErr;

   delete[] dtemp;
   delete[] itemp;
   delete[] stemp;

   delete sdrecshw;
   delete fDetGeo;
   delete fRecEvent;
}

void AdstFile::ReadAdstFile(string inname, double *energylimitFull, double *zenithlimitFull, vector<double> *energy, vector<double> *zenith, vector<double> *shwsize, int *allevts)
{
   cout << endl << "Opening file: " << inname << " -------------------------------------------------------" << endl;
   cerr << endl << "Opening file: " << inname << " -------------------------------------------------------" << endl;

   // Open and prepare the ADST files for reading
   fFile = new RecEventFile(inname.c_str(), RecEventFile::eRead);
   nevents = fFile->GetNEvents();
   cout << "Number of events: " << nevents << endl;

   fFile->SetBuffers(&fRecEvent);
   fFile->ReadDetectorGeometry(*fDetGeo);

   for(int j = 0; j < nevents; j++)
   {
      if(nevents < 20)
         cerr << "Currently at " << j << "/" << nevents << endl;
      else if(j%((int)(nevents*0.05)) == 0)
         cerr << "Currently at " << j << "/" << nevents << endl;
      else
      {
         if(j%((int)(nevents*0.1)) == 0)
            cerr << "Currently at " << j << "/" << nevents << endl;
      }
   
      fFile->ReadEvent(j);
      cout << "New event (" << j+1 << ", ID = " << fRecEvent->GetEventId() << ", Time = [" << fRecEvent->GetYYMMDD() << "," << fRecEvent->GetHHMMSS() << "], E = " << TMath::Log10(fRecEvent->GetSDEvent().GetSdRecShower().GetEnergy()) << ", sec(theta) = " <<SecTheta(fRecEvent->GetSDEvent().GetSdRecShower().GetZenith(),false)  << ") -------" << endl;
      goodrec = true;
      isHeco = false;
            
      // Prepare SD station events -------------------------------------------
      *sdrecshw = fRecEvent->GetSDEvent().GetSdRecShower();
  
      // Check if there are triggered SD stations
      if(!(fRecEvent->GetSDEvent().HasTriggeredStations()))
         goodrec = false;
      // Check if there are any SD stations in the event
      if(!(fRecEvent->GetSDEvent().HasStations()))
         goodrec = false;
      // Check if SD stations have a VEM trace
      if(!(fRecEvent->GetSDEvent().HasVEMTraces()))
         goodrec = false;

      // Prepare FD eye events -------------------------------------------
      if(fRecEvent->GetNEyes() == 0)
         goodrec = false;
      else
      {
         acteyes.clear();
	 acteyes = fRecEvent->GetFDEvents();
	 nreyes = acteyes.size();

         wQuantity[0] = 0;
         wQuantity[1] = 0;
         quantitySum[0] = 0;
         quantitySum[1] = 0;
         quantitySumErr[0] = 0;
         quantitySumErr[1] = 0;
         wQuantitySum[0] = 0;
         wQuantitySum[1] = 0;

	 for(int i = 0; i < nreyes; i++)
	 {
            if(!acteyes[i].IsHybridEvent())
               goodrec = false;
	    else
	    {
               goodrec = true;
	       break;
	    }
	 }

	 if(goodrec)
	 {
	    for(int i = 0; i < nreyes; i++)
	    {
               if(acteyes[i].GetEyeId() == 6) // for now don't use HECO
                  isHeco = true;
	    }

	    if(isHeco)
               nreyes--;

	    for(int i = 0; i < nreyes; i++)
	    {
               if(acteyes[i].IsHybridEvent() && (acteyes[i].GetEyeId() < 6)) // for now don't use HECO
	       {
//                  cout << "Eye ID = " << acteyes[i].GetEyeId() << endl;
	          dtemp[0] = (acteyes[i].GetFdRecShower()).GetEnergy();
	          dtemp[1] = (acteyes[i].GetFdRecShower()).GetEnergyError();
	          dtemp[2] = (acteyes[i].GetFdRecShower()).GetZenith();
	          dtemp[3] = (acteyes[i].GetFdRecShower()).GetZenithError();
//		  cout << "EYE " << i << ": energy = " << dtemp[0] << ", " << dtemp[1] << ", zenith = " << dtemp[2] << ", " << dtemp[3] << endl;
//		  cout << "EYE " << i << ": cos zenith = " << (acteyes[i].GetFdRecShower()).GetCosZenith() << ", " << (acteyes[i].GetFdRecShower()).GetCosZenithError() << endl;

		  wQuantity[0] = 1./TMath::Power(dtemp[1],2);
		  wQuantity[1] = 1./TMath::Power(dtemp[3],2);
		  quantitySum[0] += dtemp[0]*wQuantity[0];
		  quantitySum[1] += dtemp[2]*wQuantity[1];
		  wQuantitySum[0] += wQuantity[0];
		  wQuantitySum[1] += wQuantity[1];
		  quantitySumErr[0] += wQuantity[0];
		  quantitySumErr[1] += wQuantity[1];
	       }
	    }

	    if((nreyes == 0) && isHeco)
               goodrec = false;
	    else
            {
	       quantitySumErr[0] = TMath::Sqrt(1./quantitySumErr[0]);
	       quantitySumErr[1] = TMath::Sqrt(1./quantitySumErr[1]);
	       quantitySum[0] = quantitySum[0]/wQuantitySum[0];
	       quantitySum[1] = quantitySum[1]/wQuantitySum[1];
	    }

//	    cout << "nr eyes = " << nreyes << ": final energy = " << quantitySum[0] << ", " << quantitySumErr[0] << ", final zenith = " << quantitySum[1] << ", " << quantitySumErr[1] << endl;

            // Limit in energy
            if(energylimitFull[0] != -1)
            {
               if(TMath::Log10(quantitySum[0]) < energylimitFull[0])
                  goodrec = false;
            }
            if(energylimitFull[1] != -1)
            {
               if(TMath::Log10(quantitySum[0]) > energylimitFull[1])
                  goodrec = false;
            }

            // Full limit in zenith angle (will be further binned later)
            if(zenithlimitFull[0] != -1)
            {
               if(SecTheta(quantitySum[1],false) < zenithlimitFull[0])
                  goodrec = false;
            }
            if(zenithlimitFull[1] != -1)
            {
               if(SecTheta(quantitySum[1],false) > zenithlimitFull[1])
                  goodrec = false;
            }

	    if(goodrec)
	    {
	       energy->push_back(quantitySum[0]);
	       energy->push_back(quantitySumErr[0]);
	       zenith->push_back(quantitySum[1]);
	       zenith->push_back(quantitySumErr[1]);

	       cout << "energy = " << quantitySum[0] << ", " << quantitySumErr[0] << endl;
	       cout << "zenith = " << quantitySum[1] << ", " << quantitySumErr[1] << endl;
	    }
	 }
      }

      if(goodrec)
      {
         shwsize->push_back(sdrecshw->GetShowerSize());
         shwsize->push_back(sdrecshw->GetShowerSizeError());

	 cout << "s1000 = " << sdrecshw->GetShowerSize() << ", " << sdrecshw->GetShowerSizeError() << endl;

	 (*allevts)++;
      }
   }
   
   fFile->Close();
   delete fFile;
}

void BinNaming(string *instring, double *energy, double *zenith)
{
   if( (energy[0] != -1) || (energy[1] != -1) )
      *instring += "_en_";
   if(energy[0] != -1)
      *instring += ToString(energy[0],2);
   if(energy[1] != -1)
      *instring += "-" + ToString(energy[1],2);
   if( (zenith[0] != -1) || (zenith[1] != -1) )
      *instring += "_zen_";
   if(zenith[0] != -1)
      *instring += ToString(zenith[0],2);
   if(zenith[1] != -1)
      *instring += "-" + ToString(zenith[1],2);
}

int ReadCustomBinning(string infile, vector<double> *outBins, double *maxRange)
{
   outBins->clear();

   ifstream *ifs = new ifstream;
   double *dtemp = new double[4];
   int itemp;

   ifs->open(infile.c_str(), ifstream::in);

   dtemp[2] = 1e+20;
   dtemp[3] = -1;

   if(ifs->is_open())
   {
      itemp = 0;
      cout << "Writing out custom binning:" << endl;
      while(1)
      {
         *ifs >> dtemp[0] >> dtemp[1];
   
         if(dtemp[0] < dtemp[2])
            dtemp[2] = dtemp[0];
         if(dtemp[1] > dtemp[3])
            dtemp[3] = dtemp[1];
   
         outBins->push_back(dtemp[0]);
         outBins->push_back(dtemp[1]);
   
         itemp++;
         ifs->ignore(1,' ');
         if(ifs->eof()) break;
      }
   
      cout << "Number of bins = " << itemp << ", Minimum bin value = " << dtemp[2] << ", Maximum bin value = " << dtemp[3] << endl;
      maxRange[0] = dtemp[2];
      maxRange[1] = dtemp[3];
   }
   
   ifs->close();

   delete ifs;
   delete[] dtemp;

   return itemp;
}

// Calculate the S38 value and it's uncertainty
//   S38 = S1000/fCIC
//   dS38 = sqrt((dS1000/fCIC)^2 + (-S1000*da*x/fCIC^2)^2 + (-S1000*db*x^2/fCIC^2)^2 + (-S1000*dc*x^3/fCIC^2)^2 + (-S1000*(a+2b*x+3c*x^2)/fCIC^2)^2)
void CalculateS38(vector<double> *shwsize, vector<double> *zenith, double *fitpar, double *fitparErr, vector<double> *outVect)
{
   double *dtemp = new double[4];
   double *seczenith = new double[2];

   for(int i = 0; i < zenith->size()/2; i++)
   {
      // Calculate sec(theta) from zenith angle
      seczenith[0] = SecTheta(zenith->at(2*i),false);
      seczenith[1] = (TMath::Sin(zenith->at(2*i))*zenith->at(2*i+1))/TMath::Power(TMath::Cos(zenith->at(2*i)),2);

      cout << i << ": sec(theta) = " << seczenith[0] << " ± " << seczenith[1] << ", S1000 = " << shwsize->at(2*i) << " ± " << shwsize->at(2*i+1);
      // x = (1/sec(theta))^2 - (cos(thetaref))^2
      dtemp[0] = TMath::Power(1./seczenith[0],2) - TMath::Power(TMath::Cos(DegToRad(38.)),2);
      cout << ", x = " << dtemp[0];
      // fCIC = 1 + a*x + b*x^2 + c*x^3
      dtemp[1] = 1.+fitpar[1]*dtemp[0]+fitpar[2]*TMath::Power(dtemp[0],2)+fitpar[3]*TMath::Power(dtemp[0],3);
      cout << ", fCIC = " << dtemp[1];
      // S38 = S1000/fCIC
      dtemp[2] = shwsize->at(2*i)/dtemp[1];
      cout << ", S38 = " << dtemp[2] << endl;
   
      outVect->push_back(dtemp[2]);
   
      // (dS1000/fCIC)^2
      dtemp[2] = shwsize->at(2*i+1)/dtemp[1];
      dtemp[3] = TMath::Power(dtemp[2],2);
      cout << "  S1000 err = " << dtemp[2];
      // (-S1000*da*x/fCIC^2)^2
      dtemp[2] = -(shwsize->at(2*i)*fitparErr[1]*dtemp[0])/TMath::Power(dtemp[1],2);
      dtemp[3] += TMath::Power(dtemp[2],2);
      cout << ", a err = " << dtemp[2];
      // (-S1000*db*x^2/fCIC^2)^2
      dtemp[2] = -(shwsize->at(2*i)*fitparErr[2]*TMath::Power(dtemp[0],2))/TMath::Power(dtemp[1],2);
      dtemp[3] += TMath::Power(dtemp[2],2);
      cout << ", b err = " << dtemp[2];
      // (-S1000*dc*x^3/fCIC^2)^2
      dtemp[2] = -(shwsize->at(2*i)*fitparErr[3]*TMath::Power(dtemp[0],3))/TMath::Power(dtemp[1],2);
      dtemp[3] += TMath::Power(dtemp[2],2);
      cout << ", c err = " << dtemp[2];
      // (-S1000*dx*(a+2b*x+3c*x^2)/fCIC^2)^2
      // dx = 2*cos(theta)*sin(theta)*dtheta
      dtemp[2] = 2*seczenith[1]/TMath::Power(seczenith[0],3);
      cout << ", dx = " << dtemp[2];
      dtemp[2] = -(shwsize->at(2*i)*dtemp[2]*(fitpar[1]+2*fitpar[2]*dtemp[0]+3*fitpar[3]*TMath::Power(dtemp[0],2)))/TMath::Power(dtemp[1],2);
      dtemp[3] += TMath::Power(dtemp[2],2);
      cout << ", x err = " << dtemp[2];
   
      dtemp[3] = TMath::Sqrt(dtemp[3]);
      cout << ", dS38 = " << dtemp[3] << endl;
      outVect->push_back(dtemp[3]);
   }

   delete[] dtemp;
   delete[] seczenith;
}

int main(int argc, char **argv)
{
   gSystem->Load("libTree.so");

   string inname;
   string *stemp = new string[3];
   int *itemp = new int[3];
   double *dtemp = new double[7];
   TGraphErrors *result[20];
   RootStyle *mystyle;
   TCanvas *c1;
   TLine *line = new TLine();
   line->SetLineWidth(2);
   line->SetLineStyle(7);
   line->SetLineColor(28);
   TF1 *fitfunc[6];
   TF1 *tempfunc;
   double fitparam[4];
   double fitparamErr[4];
   double tempparam[4];
   double tempparamErr[4];

   TGraphErrors *allS38, *alldeltaS38;
   int s38count = 0;

   TGraphErrors *allFitparam[3];
   int fitparamcount = 0;

   vector<double> *energy = new vector<double>;
   vector<double> *zenith = new vector<double>;
   vector<double> *shwsize = new vector<double>;
   vector<double> *s38 = new vector<double>;
   vector<double> *deltas38 = new vector<double>;
   vector<double> energyVect;
   vector<double> zenithVect;
   vector<double> seczenithVect;
   vector<double> shwsizeVect;
   vector<double> s38Vect;

   int allevts = 0;
   bool goodrec;

   // Variables for setting energy and zenith angle limits
   double energylimit[2] = {-1,-1};
   double energylimitFull[2] = {-1,-1};
   int energynr = 12;
   vector<double> *energybins = new vector<double>;
   int energyref = -1;
   double zenithlimit[2] = {-1,-1};
   double zenithlimitFull[2] = {-1,-1};
   int zenithnr = 9;
   vector<double> *zenithbins = new vector<double>;

   stemp[0] = string(rootdir) + "/input/custom_energy_bins.txt";
   energynr = ReadCustomBinning(stemp[0], energybins, energylimitFull);
   stemp[0] = string(rootdir) + "/input/custom_zenith_bins.txt";
   zenithnr = ReadCustomBinning(stemp[0], zenithbins, zenithlimitFull);

   AdstFile *adfile = new AdstFile();

   if(argc > 1)
   {
      // Prepare plotting style
      mystyle = new RootStyle();
      mystyle->SetBaseStyle();
      c1 = new TCanvas("c1","",1200,900);
      gStyle->SetEndErrorSize(3);

      // Read all FD energy, FD zenith and S1000 values inside the energy and zenith angle ranges
      for(int k = 0; k < argc-1; k++)
      {
         inname = string(argv[k+1]);
         adfile->ReadAdstFile(inname, energylimitFull, zenithlimitFull, energy, zenith, shwsize, &allevts);
      }

      cout << "Surviving a total of " << allevts << " valid events." << endl;
      cerr << "Surviving a total of " << allevts << " valid events." << endl;

      allS38 = new TGraphErrors();
      allFitparam[0] = new TGraphErrors();
      allFitparam[1] = new TGraphErrors();
      allFitparam[2] = new TGraphErrors();
      s38count = 0;
      fitparamcount = 0;
      s38->clear();

      // Go through all energy zenith bins
      for(int iEn = 0; iEn < energynr; iEn++)
      {
         energylimit[0] = energybins->at(2*iEn);
         energylimit[1] = energybins->at(2*iEn+1);

         cout << "Chosen energy bin = " << energylimit[0] << ", " << energylimit[1] << endl;

	 // Vectors for energy, zenith angle and S1000 of separate energy bins
	 energyVect.clear();
	 zenithVect.clear();
	 seczenithVect.clear();
	 shwsizeVect.clear();
	 s38Vect.clear();

	 // Go through all selected events
	 for(int i = 0; i < shwsize->size()/2; i++)
	 {
            goodrec = true;

	    // Select only events in the energy bin
            if(TMath::Log10(energy->at(2*i)) < energylimit[0])
               goodrec = false;
            if(TMath::Log10(energy->at(2*i)) > energylimit[1])
               goodrec = false;

	    if(goodrec)
	    {
               cout << "Event " << i << ", energy = " << energy->at(2*i) << " (" << energy->at(2*i+1)  << "), zenith = " << SecTheta(zenith->at(2*i),false) << " (" << (TMath::Sin(zenith->at(2*i))*zenith->at(2*i+1))/TMath::Power(TMath::Cos(zenith->at(2*i)),2) << "), shwsize = " << shwsize->at(2*i) << " (" << shwsize->at(2*i+1) << ")" << endl;

	       // Push back all three values and their uncertainties
	       energyVect.push_back(energy->at(2*i));
	       energyVect.push_back(energy->at(2*i+1));
	       zenithVect.push_back(zenith->at(2*i));
	       zenithVect.push_back(zenith->at(2*i+1));
	       shwsizeVect.push_back(shwsize->at(2*i));
	       shwsizeVect.push_back(shwsize->at(2*i+1));
	       // Convert theta to sec(theta):
	       //   sec(theta) = 1/cos(theta)
	       //   dsec(theta) = sin(theta)*dtheta/(cos(theta))^2
	       seczenithVect.push_back(SecTheta(zenith->at(2*i),false));
	       seczenithVect.push_back((TMath::Sin(zenith->at(2*i))*zenith->at(2*i+1))/TMath::Power(TMath::Cos(zenith->at(2*i)),2));
	    }
	 }

	 allevts = 0;

         result[0] = new TGraphErrors();
	 dtemp[4] = 0.;

	 // Save values for sec(theta)/S1000 to graph
	 for(int i = 0; i < zenithVect.size()/2; i++)
         {
            result[0]->SetPoint(i, seczenithVect[2*i], shwsizeVect[2*i]);
            result[0]->SetPointError(i, seczenithVect[2*i+1], shwsizeVect[2*i+1]);

            dtemp[4] += shwsizeVect[2*i];

	    allevts++;
	 }

	 cout << "Total of " << allevts << " events in chosen energy bin." << endl;
	 
	 if(allevts > 0)
         {
	    // Mean value of all S1000 event signals for initial fitting parameter S
	    dtemp[4] = dtemp[4]/allevts;
	    cout << "Mean value of S1000 = " << dtemp[4] << endl;

            // Plotting S1000 vs sec(theta)
            mystyle->SetGraphColor(result[0], 1);
	    result[0]->SetMarkerStyle(20);
	    result[0]->SetMarkerSize(0.8);
            mystyle->SetAxisTitles(result[0], "FD zenith angle (sec#theta)", "S_{1000} (VEM)");

	    c1->SetLogx(kFALSE);
	    c1->SetLogy(kFALSE);
	    result[0]->Draw("AP");

	    // Fitting function fCIC: x = sec(theta) -> fCIC = S + A*x + B*x^2 + C*x^3
            fitfunc[0] = new TF1("fitfunc0", "[0]+[1]*(TMath::Power(1./x,2)-TMath::Power(cos(0.663225115),2))+[2]*TMath::Power((TMath::Power(1./x,2)-TMath::Power(cos(0.663225115),2)),2)+[3]*TMath::Power((TMath::Power(1./x,2)-TMath::Power(cos(0.663225115),2)),3)", 1., 2.);
	    fitfunc[0]->SetParameters(dtemp[4],dtemp[4]*1.,-dtemp[4]*1.5,-dtemp[4]*1.3);
	    result[0]->Fit("fitfunc0");
	    c1->Update();
	    line->DrawLine(SecTheta(38.,true), gPad->GetUymin(), SecTheta(38.,true), gPad->GetUymax());

	    tempfunc = (TF1*)result[0]->GetFunction("fitfunc0");
	    tempfunc->Draw("L;SAME");

	    // Convert fitting parameters to: a = A/S, b = B/S, c = C/S
	    //   in order to convert the function to fCIC = 1 + a*x + b*x^2 + c*x^3
	    //   uncertainties are: da = sqrt((dA/S)^2 + (-A*dS/S^2)^2)
	    //                      db = sqrt((dB/S)^2 + (-B*dS/S^2)^2)
	    //                      dc = sqrt((dC/S)^2 + (-C*dS/S^2)^2)
	    for(int j = 0; j < 4; j++)
	    {
               tempparam[j] = fitfunc[0]->GetParameter(j);
               tempparamErr[j] = fitfunc[0]->GetParError(j);

	       if(j == 0)
	       {
                  fitparam[0] = tempparam[0];
                  fitparamErr[0] = tempparamErr[0];
	       }
	       else
	       {
                  fitparam[j] = tempparam[j]/tempparam[0];
//	          fitparamErr[j] = TMath::Sqrt(TMath::Power(fitparamErr[j]/fitparam[0],2) + TMath::Power(-(fitparam[j]*fitparamErr[0])/TMath::Power(fitparam[0],2),2));
	          fitparamErr[j] = TMath::Sqrt(TMath::Power(tempparamErr[j]/fitparam[0],2) + TMath::Power(-(tempparam[j]*fitparamErr[0])/TMath::Power(fitparam[0],2),2));
	       }
	    }

	    cout << endl << "Fitting parameters of fCIC are (chi2/ndf = " << fitfunc[0]->GetChisquare() << "/" << fitfunc[0]->GetNDF() << " = " << (fitfunc[0]->GetChisquare())/(fitfunc[0]->GetNDF()) << "):" << endl;
	    cout << "- S1000 at 38 deg = " << fitfunc[0]->Eval(SecTheta(38.,true)) << endl;
	    cout << "- S = " << fitparam[0] << " ± " << fitparamErr[0] << endl;
	    cout << "- a = " << fitparam[1] << " ± " << fitparamErr[1] << endl;
	    cout << "- b = " << fitparam[2] << " ± " << fitparamErr[2] << endl;
	    cout << "- c = " << fitparam[3] << " ± " << fitparamErr[3] << endl;
	    cout << endl;

	    stemp[0] = "rm ./plots/s1000_vs_zenithFD";
            BinNaming(&stemp[0], energylimit, zenithlimitFull);
            stemp[0] += "*";
            system(stemp[0].c_str());

	    stemp[0] = "./plots/s1000_vs_zenithFD";
            BinNaming(&stemp[0], energylimit, zenithlimitFull);
            stemp[1] = stemp[0] + ".pdf";
	    c1->SaveAs(stemp[1].c_str());

	    // Save fitting parameters for plotting
	    for(int i = 0; i < 3; i++)
	    {
	       allFitparam[i]->SetPoint(fitparamcount, (energylimit[0]+energylimit[1])/2., fitparam[i+1]);
	       allFitparam[i]->SetPointError(fitparamcount, 0., fitparamErr[i+1]);
	    }
	    fitparamcount++;

	    // Calculate the S38 value and it's uncertainty
	    //   S38 = S1000/fCIC
	    //   dS38 = sqrt((dS1000/fCIC)^2 + (-S1000*da*x/fCIC^2)^2 + (-S1000*db*x^2/fCIC^2)^2 + (-S1000*dc*x^3/fCIC^2)^2 + (-S1000*(a+2b*x+3c*x^2)/fCIC^2)^2)
	    cout << "Calculating S38:" << endl;
	    CalculateS38(&shwsizeVect, &zenithVect, fitparam, fitparamErr, &s38Vect);

            result[1] = new TGraphErrors();
            result[2] = new TGraphErrors();

	    for(int i = 0; i < zenithVect.size()/2; i++)
	    {
               result[1]->SetPoint(i, seczenithVect[2*i], s38Vect[2*i]);
               result[1]->SetPointError(i, seczenithVect[2*i+1], s38Vect[2*i+1]);

               result[2]->SetPoint(i, energyVect[2*i]/1.e+18, s38Vect[2*i]);
               result[2]->SetPointError(i, energyVect[2*i+1]/1.e+18, s38Vect[2*i+1]);

               allS38->SetPoint(i, energyVect[2*i]/1.e+18, s38Vect[2*i]);
               allS38->SetPointError(i, energyVect[2*i+1]/1.e+18, s38Vect[2*i+1]);
	       s38count++;
	    }

            // Plotting S38 vs sec(theta)
            mystyle->SetGraphColor(result[1], 1);
	    result[1]->SetMarkerStyle(20);
	    result[1]->SetMarkerSize(0.8);
            mystyle->SetAxisTitles(result[1], "FD zenith angle (sec#theta)", "S_{38} (VEM)");

	    c1->SetLogx(kFALSE);
	    c1->SetLogy(kFALSE);
	    result[1]->Draw("AP");

	    stemp[0] = "rm ./plots/s38_vs_zenithFD";
            BinNaming(&stemp[0], energylimit, zenithlimitFull);
            stemp[0] += "*";
            system(stemp[0].c_str());

	    stemp[0] = "./plots/s38_vs_zenithFD";
            BinNaming(&stemp[0], energylimit, zenithlimitFull);
            stemp[1] = stemp[0] + ".pdf";
	    c1->SaveAs(stemp[1].c_str());

            // Plotting S38 vs energy
            mystyle->SetGraphColor(result[2], 1);
	    result[2]->SetMarkerStyle(20);
	    result[2]->SetMarkerSize(0.8);
            mystyle->SetAxisTitles(result[2], "FD energy (EeV)", "S_{38} (VEM)");

	    c1->SetLogx(kFALSE);
	    c1->SetLogy(kFALSE);
	    result[2]->Draw("AP");

	    stemp[0] = "rm ./plots/s38_vs_energyFD";
            BinNaming(&stemp[0], energylimit, zenithlimitFull);
            stemp[0] += "*";
            system(stemp[0].c_str());

	    stemp[0] = "./plots/s38_vs_energyFD";
            BinNaming(&stemp[0], energylimit, zenithlimitFull);
            stemp[1] = stemp[0] + ".pdf";
	    c1->SaveAs(stemp[1].c_str());

	    delete result[1];
	    delete result[2];
	    delete fitfunc[0];
	 }

         delete result[0];
      }

      // Plot the energy dependence of fCIC fitting parameters
      stemp[0] = "rm ./plots/energyFD_vs_fitparam*";
      system(stemp[0].c_str());

      for(int i = 0; i < 3; i++)
      {
         mystyle->SetGraphColor(allFitparam[i], 1);
         allFitparam[i]->SetMarkerStyle(20);
         allFitparam[i]->SetMarkerSize(0.8);

	 if(i == 0)
            mystyle->SetAxisTitles(allFitparam[i], "FD energy (log(E/eV))", "f_{CIC} fitting parameter a [f_{CIC} = 1 + a*x + b*x^{2} + c*x^{3}]");
	 else if(i == 1)
            mystyle->SetAxisTitles(allFitparam[i], "FD energy (log(E/eV))", "f_{CIC} fitting parameter b [f_{CIC} = 1 + a*x + b*x^{2} + c*x^{3}]");
	 else if(i == 2)
            mystyle->SetAxisTitles(allFitparam[i], "FD energy (log(E/eV))", "f_{CIC} fitting parameter c [f_{CIC} = 1 + a*x + b*x^{2} + c*x^{3}]");
         
         c1->SetLogx(kFALSE);
         c1->SetLogy(kFALSE);

         allFitparam[i]->Draw("AP");
	 if(i == 0)
            stemp[1] = "./plots/energyFD_vs_fitparamA.pdf";
	 else if(i == 1)
            stemp[1] = "./plots/energyFD_vs_fitparamB.pdf";
	 else if(i == 2)
            stemp[1] = "./plots/energyFD_vs_fitparamC.pdf";
         c1->SaveAs(stemp[1].c_str());
      }

      // Plot the S38 values for all energies
      mystyle->SetGraphColor(allS38, 1);
      allS38->SetMarkerStyle(20);
      allS38->SetMarkerSize(0.8);
      mystyle->SetAxisTitles(allS38, "FD energy (EeV)", "S_{38} (VEM)");
     
      c1->SetLogx(kTRUE);
      c1->SetLogy(kTRUE);
/*      allS38->GetYaxis()->SetRange(4.,500.);
      allS38->GetYaxis()->SetRangeUser(4.,500.);*/
      allS38->GetXaxis()->SetMoreLogLabels(kTRUE);
      allS38->Draw("AP");
      fitfunc[1] = new TF1("fitfunc1", "TMath::Power(x/[0],1./[1])", 3., 100.);
      fitfunc[1]->SetParameters(0.1,1.);
      allS38->Fit("fitfunc1");
      tempfunc = (TF1*)allS38->GetFunction("fitfunc1");
      tempfunc->Draw("L;SAME");

      stemp[0] = "rm ./plots/s38_vs_energyFD_all*";
      system(stemp[0].c_str());
      stemp[1] = "./plots/s38_vs_energyFD_all.pdf";
      c1->SaveAs(stemp[1].c_str());
     
      cout << endl << "Fitting parameters for linear dependence between FD energy and S38:" << endl;
      cout << "- A = " << fitfunc[1]->GetParameter(0) << " ± " << fitfunc[1]->GetParError(0) << endl;
      cout << "- B = " << fitfunc[1]->GetParameter(1) << " ± " << fitfunc[1]->GetParError(1) << endl;
      cout << endl;

      // Prepare reference value of S38 (at energy 10^19)
      alldeltaS38 = new TGraphErrors();
      // Reference value of S38 from the fitting function at 10^19
      dtemp[0] = fitfunc[1]->Eval(10.);
      cout << "S38 at reference energy of 10 EeV = " << dtemp[0] << endl;
      deltas38->clear();
      cout << "Number of points = " << allS38->GetN() << endl;
      cout << "Energy points = " << energy->size()/2 << endl;
      dtemp[5] = 0.;
      dtemp[6] = 0.;
      for(int i = 0; i < allS38->GetN(); i++)
      {
         allS38->GetPoint(i, dtemp[1], dtemp[2]);
         dtemp[3] = dtemp[2] - fitfunc[1]->Eval(dtemp[1]);
	 deltas38->push_back(dtemp[3]);
         dtemp[3] = allS38->GetErrorX(i);
	 dtemp[4] = allS38->GetErrorY(i);
	 deltas38->push_back(dtemp[4]);

	 cout << "Point " << i << ": energy = " << dtemp[1] << " ± " << dtemp[3] << ", s38 = " << dtemp[2] << " ± " << dtemp[4] << ", deltas38 = " << deltas38->at(2*i) << " ± " << deltas38->at(2*i+1) << endl; 

         if(TMath::Abs(deltas38->at(2*i)) > dtemp[5])
            dtemp[5] = TMath::Abs(deltas38->at(2*i));
         if(deltas38->at(2*i+1) > dtemp[6])
	    dtemp[6] = deltas38->at(2*i+1);

	 alldeltaS38->SetPoint(i, dtemp[1], deltas38->at(2*i));
	 alldeltaS38->SetPointError(i, dtemp[3], deltas38->at(2*i+1));
      }

      // Plot the deltaS38 values for all energies
      mystyle->SetGraphColor(alldeltaS38, 1);
      alldeltaS38->SetMarkerStyle(20);
      alldeltaS38->SetMarkerSize(0.8);
      mystyle->SetAxisTitles(alldeltaS38, "FD energy (EeV)", "#deltaS_{38} (VEM)");
     
      c1->SetLogx(kTRUE);
      c1->SetLogy(kFALSE);
      cout << "Max value = " << dtemp[5] << " + " << dtemp[6] << endl;
      dtemp[5] += dtemp[6];
      dtemp[5] += dtemp[5]*0.05;
      alldeltaS38->GetYaxis()->SetRange(-dtemp[5],dtemp[5]);
      alldeltaS38->GetYaxis()->SetRangeUser(-dtemp[5],dtemp[5]);
      alldeltaS38->GetXaxis()->SetMoreLogLabels(kTRUE);
      alldeltaS38->Draw("AP");

      stemp[0] = "rm ./plots/deltas38_vs_energyFD_all*";
      system(stemp[0].c_str());
      stemp[1] = "./plots/deltas38_vs_energyFD_all.pdf";
      c1->SaveAs(stemp[1].c_str());

      delete allS38;
      delete alldeltaS38;
      delete fitfunc[1];
      delete allFitparam[0];
      delete allFitparam[1];
      delete allFitparam[2];

      delete mystyle;
      delete c1;
   }
   else
   {
      cerr << "Error! No input files supplied. Rerun program and add input files as arguments (mvatree_file.root)." << endl;
      return 1;
   }

   delete adfile;

   delete energy;
   delete zenith;
   delete shwsize;
   delete s38;
   delete deltas38;

   delete[] stemp;
   delete[] itemp;
   delete[] dtemp;

   return 0;
}
