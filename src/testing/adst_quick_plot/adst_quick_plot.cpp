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

/*class AdstFile
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

   // Temporary variables and holders for number of stations and events
   int nrstations, nevents;
   bool goodrec;

   // Limits for risetime calculations
   double limitTankDistance[2];
   double minSignal;
   bool includeSaturated;
   int minPoints;

   // Vectors and variables needed for calculation of risetime
   vector<double> tempVect;
   vector<float> time;
   vector<float> vemtrace;
   vector<float> yvem;
   vector<float> yintvem;
   vector<float> yvalue;
   vector<float> tempRise;
   int start_bin, stop_bin;
   double maxval;
   double xp, yp;
   double risemean, riseerr;
   double *byrange;
   double *bzrange;
   TFormula *fRTWeights;

   double eventThetaRec;
   double secZenith;
   double alpha;
   double gamma;
   double g;
   double zeta;

   RootStyle *mystyle;

   const int c_SignalLine = TColor::GetColor("#0000ee");
   const int c_SignalFill = TColor::GetColor("#7d99d1");
   const int c_BackgroundLine = TColor::GetColor("#ff0000");
   const int c_BackgroundFill = TColor::GetColor("#ff0000");
   const int c_DataLine = TColor::GetColor("#000000");
   const int c_DataFill = TColor::GetColor("#808080");
   const int c_DataNormLine = TColor::GetColor("#00bb00");
   const int c_DataNormFill = TColor::GetColor("#00bb00");
   const int c_ResidLine = TColor::GetColor("#000000");
   const int c_ResidFill = TColor::GetColor("#808080");

public:
   AdstFile();
   virtual ~AdstFile();

   void ReadAdstFile(string inname, vector<double> *sdidVect, vector<bool> *HGsat, vector<double> *distVect, vector<double> *riseVect, vector<double> *energy, vector<double> *zenith, vector<double> *shwsize, vector<int> *eventVect, int *allcount, int *allevts, int *nrrun);

   void PlotCurrentVemTrace(int eventNr, int stationID, int pmtID, double max);
};

AdstFile::AdstFile()
{
   fRecEvent = new RecEvent();
   fDetGeo = new DetectorGeometry();
   sdrecshw = new SdRecShower();

   stemp = new string[3];
   itemp = new int[4];
   dtemp = new double[3];

   limitTankDistance[0] = 300.;
   limitTankDistance[1] = 1400.;
   minSignal = 5.0;
   includeSaturated = false;
   minPoints = 3;

   byrange = new double[2];
   bzrange = new double[2];
   fRTWeights = new TFormula("RiseTimeWeights", "(80.0+(5.071e-7+6.48e-4*y-3.051e-4*y*y)*x*x)/z-16.46*y+36.16");

   byrange[0] = 1.e+40;
   byrange[1] = -1.e+40;
   bzrange[0] = 1.e+40;
   bzrange[1] = -1.e+40;
}

AdstFile::~AdstFile()
{
   delete fRTWeights;
   delete[] byrange;
   delete[] bzrange;

   delete[] dtemp;
   delete[] itemp;
   delete[] stemp;

   delete sdrecshw;
   delete fDetGeo;
   delete fRecEvent;
}*/

/*double AdstFile::GetDistanceLimit(int type)
{
   return limitTankDistance[type];
}*/

/*void AdstFile::PlotCurrentVemTrace(int eventNr, int stationID, int pmtID, double max)
{
   mystyle = new RootStyle();
   mystyle->SetBaseStyle();
   TCanvas *c1 = new TCanvas("c1", "", 1200, 900);

   TGraph *grVem = new TGraph();
   TGraph *grIntVemNorm = new TGraph();
   TGraph *grIntVem = new TGraph();

   TH1F *histVem = new TH1F("histVem", "", yvem.size()-1, 0, time[yvem.size()-1]-time[0]);
   TH1F *histIntVemNorm = new TH1F("histIntVemNorm", "", yintvem.size()-1, 0, time[yintvem.size()-1]-time[0]);
   TH1F *histIntVem = new TH1F("histIntVem", "", yintvem.size()-1, 0, time[yintvem.size()-1]-time[0]);

   dtemp[2] = 0;
   for(int i = 0; i < yvem.size(); i++)
   {
      grVem->SetPoint(i, time[i]-time[0], yvem[i]);

      if(yvem[i] < 0)
         histVem->SetBinContent(i, 0);
      else
         histVem->SetBinContent(i, yvem[i]);

      if(yvem[i] > dtemp[2])
         dtemp[2] = yvem[i];
   }

   for(int i = 0; i < yintvem.size(); i++)
   {
      grIntVemNorm->SetPoint(i, time[i]-time[0], (yintvem[i]*dtemp[2])/max);
      grIntVem->SetPoint(i, time[i]-time[0], yintvem[i]/max);

      if(yintvem[i] < 0)
      {
         histIntVemNorm->SetBinContent(i, 0);
         histIntVem->SetBinContent(i, 0);
      }
      else
      {
         histIntVemNorm->SetBinContent(i, (yintvem[i]*dtemp[2])/max);
         histIntVem->SetBinContent(i, yintvem[i]/max);
      }
   }

   stemp[0] = "mkdir -p ./plots";
   system(stemp[0].c_str());

   mystyle->SetSinglePlot(0, -1, c1);

   // Create the PMT signal vs. integrated PMT signal for whole range
   mystyle->SetHistColor(histVem, 2);
   histVem->GetXaxis()->SetRangeUser(0, time[yvem.size()-1]-time[0]);
   histVem->GetYaxis()->SetRangeUser(0, dtemp[2]*1.2);
   histVem->GetYaxis()->SetTitleOffset(mystyle->GetSingleYoffset(c1));
   histVem->GetXaxis()->SetTitleOffset(mystyle->GetSingleXoffset(c1));
   mystyle->SetAxisTitles((TH1*)histVem, "Time (ns)", "PMT signal (VEM)");
   histVem->Draw();

   mystyle->SetHistColor(histIntVemNorm, 1);
   histIntVemNorm->SetFillColorAlpha(c_SignalFill, 0.5);
   histIntVemNorm->GetXaxis()->SetRangeUser(0, time[yintvem.size()-1]-time[0]);
   histIntVemNorm->GetYaxis()->SetRangeUser(0, dtemp[2]*1.2);
   histIntVemNorm->GetYaxis()->SetTitleOffset(mystyle->GetSingleYoffset(c1));
   histIntVemNorm->GetXaxis()->SetTitleOffset(mystyle->GetSingleXoffset(c1));
   mystyle->SetAxisTitles((TH1*)histIntVemNorm, "Time (ns)", "PMT signal (VEM)");
   histIntVemNorm->Draw("SAME");

   c1->Update();

   TLine *line = new TLine();
   line->SetLineWidth(2);
   line->SetLineStyle(7);
   line->SetLineColor(2);
   line->DrawLine(tempRise[0]-time[0], gPad->GetUymin(), tempRise[0]-time[0], gPad->GetUymax());
   line->DrawLine(tempRise[1]-time[0], gPad->GetUymin(), tempRise[1]-time[0], gPad->GetUymax());

   TLatex *risetext = new TLatex();
   risetext->SetTextAlign(31);
   stemp[0] = "t_{i} = " + ToString(tempRise[1]-tempRise[0], 2) + " ns";
   risetext->DrawLatex(gPad->GetUxmax()-(0.03*(gPad->GetUxmax()-gPad->GetUxmin())), gPad->GetUymax()-(0.065*(gPad->GetUymax()-gPad->GetUymin())), stemp[0].c_str());
   
   stemp[0] = "./plots/current_vem_trace_event-" + ToString(eventNr) + "_pmt-" + ToString(pmtID) + "_station-" + ToString(stationID) + ".pdf";
   c1->SaveAs(stemp[0].c_str());

   // Create integrated PMT signal close to risetime counts
   mystyle->SetSinglePlot(0, -1, c1);

   mystyle->SetHistColor(histIntVem, 1);
   histIntVem->SetFillColorAlpha(c_SignalFill, 0.5);
   histIntVem->GetXaxis()->SetRangeUser(tempRise[0]-time[0]-100., tempRise[1]-time[0]+100.);
   histIntVem->GetYaxis()->SetRangeUser(0, 0.7);
   histIntVem->GetYaxis()->SetTitleOffset(mystyle->GetSingleYoffset(c1));
   histIntVem->GetXaxis()->SetTitleOffset(mystyle->GetSingleXoffset(c1));
   mystyle->SetAxisTitles((TH1*)histIntVem, "Time (ns)", "Integrated PMT signal (normalized)");
   histIntVem->Draw();

   mystyle->SetGraphColor(grIntVem, 0);
   grIntVem->GetXaxis()->SetRangeUser(tempRise[0]-time[0]-100., tempRise[1]-time[0]+100.);
   grIntVem->GetYaxis()->SetRangeUser(0, 0.7);
   grIntVem->GetYaxis()->SetTitleOffset(mystyle->GetSingleYoffset(c1));
   grIntVem->GetXaxis()->SetTitleOffset(mystyle->GetSingleXoffset(c1));
   mystyle->SetAxisTitles(grIntVem, "Time (ns)", "Integrated PMT signal (normalized)");
   grIntVem->Draw("L;SAME");

   c1->Update();

   line->DrawLine(tempRise[0]-time[0], gPad->GetUymin(), tempRise[0]-time[0], tempRise[2]);
   line->DrawLine(tempRise[1]-time[0], gPad->GetUymin(), tempRise[1]-time[0], tempRise[3]);

   TMarker *mark = new TMarker();
   mark->SetMarkerSize(1);
   mark->SetMarkerColor(2);
   mark->SetMarkerStyle(20);
   mark->DrawMarker(tempRise[0]-time[0], tempRise[2]);
   mark->DrawMarker(tempRise[1]-time[0], tempRise[3]);

   stemp[0] = "t_{i} = " + ToString(tempRise[1]-tempRise[0], 2) + " ns";
   risetext->DrawLatex(gPad->GetUxmax()-(0.03*(gPad->GetUxmax()-gPad->GetUxmin())), gPad->GetUymax()-(0.065*(gPad->GetUymax()-gPad->GetUymin())), stemp[0].c_str());
   
   stemp[0] = "./plots/current_vem_tracezoom_event-" + ToString(eventNr) + "_pmt-" + ToString(pmtID) + "_station-" + ToString(stationID) + ".pdf";
   c1->SaveAs(stemp[0].c_str());

   delete risetext;
   delete line;
   delete mark;
   delete histVem;
   delete histIntVemNorm;
   delete histIntVem;
   delete grVem;
   delete grIntVemNorm;
   delete grIntVem;
   delete c1;
   delete mystyle;
}*/

/*void AdstFile::ReadAdstFile(string inname, vector<double> *sdidVect, vector<bool> *HGsat, vector<double> *distVect, vector<double> *riseVect, vector<double> *energy, vector<double> *zenith, vector<double> *shwsize, vector<int> *eventVect, int *allcount, int *allevts, int *nrrun)
{
   cout << endl << "Opening file: " << inname << " -------------------------------------------------------" << endl;
   cerr << endl << "Opening file: " << inname << " -------------------------------------------------------" << endl;

   // Open and prepare the ADST files for reading
   fFile = new RecEventFile(inname.c_str(), RecEventFile::eRead);
   nevents = fFile->GetNEvents();
   cout << "Number of events: " << nevents << endl;

   fFile->SetBuffers(&fRecEvent);
   fFile->ReadDetectorGeometry(*fDetGeo);

   for(int j = 0; j < *nrrun; j++)
   {
      fFile->ReadEvent(j);
      cout << "New event (" << j+1 << ", ID = " << fRecEvent->GetEventId() << ", Time = [" << fRecEvent->GetYYMMDD() << "," << fRecEvent->GetHHMMSS() << "], E = " << TMath::Log10(fRecEvent->GetSDEvent().GetSdRecShower().GetEnergy()) << ", sec(theta) = " <<SecTheta(fRecEvent->GetSDEvent().GetSdRecShower().GetZenith(),false)  << ") -------" << endl;
      goodrec = true;
            
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

      if(goodrec)
      {
         dtemp[0] = 0;
         itemp[0] = 0;

         // Loop over all triggered SD stations
         stationVector = fRecEvent->GetSDEvent().GetStationVector();
         itemp[3] = 0;
         tempVect.clear();
         for(int i = 0; i < stationVector.size(); i++)
         {
            // Only use stations that are valid candidates
            if(stationVector[i].IsCandidate())
            {
               cout << "New station (" << stationVector[i].GetId() << ")" << endl;	// DEBUG
               start_bin = stationVector[i].GetSignalStartSlot() - 4;
               stop_bin = stationVector[i].GetSignalEndSlot();
   
               if( (start_bin >= stop_bin) || (start_bin < 0) || (start_bin > 5000) )
                  start_bin = 0;
   
               dtemp[1] = 0;
               itemp[1] = 0;
   
               // Check all PMTs
               for(int iPMT = 1; iPMT <= 3; iPMT++)
               {
                  time.clear();
                  yvem.clear();
                  yintvem.clear();
                  yvalue.clear();
                  vemtrace.clear();
   
                  yp = 0;
                  maxval = -1.e40;
   
                  vemtrace = stationVector[i].GetVEMTrace(iPMT);
   
                  //cout << "PMT " << iPMT << ": Number of points in the VEM trace: " << vemtrace.size() << " --------------------------------------------------------" << endl;	// DEBUG
   
                  // Continue if there is a VEM trace
                  if(vemtrace.size() > 0)
                  {
                     itemp[0]++;
                     itemp[1]++;
   
                     dtemp[2] = 0;
   
                     // Prepare the time vector (each point is multiplied by 25 to get nanoseconds)
		     cout << "VEM trace printout + integrated VEM trace for PMT " << iPMT << " and station ID " << stationVector[i].GetId() << endl;
		     cout << "time\tVEM trace\tIntegrated VEM trace" << endl;
                     for(int iVEM = 0; iVEM < vemtrace.size(); iVEM++)
                     {
                        if( (iVEM >= start_bin) && (iVEM <= stop_bin) )
                        {
                           time.push_back((float)iVEM*25.);
   
                           yp += vemtrace[iVEM];
                           if(yp > maxval)
                              maxval = yp;
                        
                           yvem.push_back(vemtrace[iVEM]);
                           yintvem.push_back(yp);
                           yvalue.push_back(yp);
                           dtemp[2] += yp;

			   cout << iVEM*25. << "\t" << vemtrace[iVEM] << "\t" << yp << endl;
                        }
                     }
		     cout << endl;
   
                     //cout << "Number of points in the signal slot: " << yvalue.size() << ", Maxval: " << maxval << endl;	// DEBUG
   
                     if(dtemp[2] < 0)
                     {
                        cout << "Rejected PMT " << iPMT << " in tank " << stationVector[i].GetId() << ": Negative signal integral value = " << dtemp[2] << endl;
                        itemp[0]--;
                        itemp[1]--;
                     }
                     else
                     {
                        for(int iy = 0; iy < yvalue.size(); iy++)
                        {
                           //cout << time[iy]/25. << "\t" << yvalue[iy]/maxval << endl;	// DEBUG
   
                           if(yvalue[iy]/maxval > 0.95)
                              break;
   
                           if(yvalue[iy]/maxval <= 0.10)
                           {
                              byrange[0] = yvalue[iy]/maxval;
                              byrange[1] = yvalue[iy+1]/maxval;
   
                              yp = 0.1;
                              //cout << "yp = " << yp << ", byrange = " << byrange[0] << ", " << byrange[1] << ", time = " << time[iy] << ", " << time[iy+1] << endl;	// DEBUG
                              // Find the x value of point with y value = yp = 0.1, that lies on a line between two points
                              // y = k*x + a
                              //    k = (y2 - y1)/(x2 - x1)
                              //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                              // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                              xp = ((time[iy+1] - time[iy])*(yp - byrange[1]))/(byrange[1] - byrange[0]) + time[iy+1];
   
                              byrange[0] = xp;
                              byrange[1] = yp;
                           }
   
                           if(yvalue[iy]/maxval <= 0.50)
                           {
                              bzrange[0] = yvalue[iy]/maxval;
                              bzrange[1] = yvalue[iy+1]/maxval;
   
                              yp = 0.5;
                              // Find the x value of point with y value = yp = 0.5, that lies on a line between two points
                              // y = k*x + a
                              //    k = (y2 - y1)/(x2 - x1)
                              //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                              // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                              xp = ((time[iy+1] - time[iy])*(yp - bzrange[1]))/(bzrange[1] - bzrange[0]) + time[iy+1];
   
                              bzrange[0] = xp;
                              bzrange[1] = yp;
                           }
                        }
   
                        cout << "Calculated risetime (" << byrange[0]/25. << "," << bzrange[0]/25. << ") = " << bzrange[0] - byrange[0] << endl;	// DEBUG

			tempRise.clear();
			tempRise.push_back(byrange[0]);
			tempRise.push_back(bzrange[0]);
			tempRise.push_back(byrange[1]);
			tempRise.push_back(bzrange[1]);
		        PlotCurrentVemTrace(j+1, stationVector[i].GetId(), iPMT, maxval);
   
                        dtemp[0] += bzrange[0] - byrange[0];
                        dtemp[1] += bzrange[0] - byrange[0];
             	     }
                  }
               }
   
               dtemp[1] = dtemp[1]/itemp[1];
   
               cout << "Station " << stationVector[i].GetId() << ", " << stationVector[i].GetSPDistance() << " m: Calculated average risetime (for " << itemp[1] << " PMTs in the tank) = " << dtemp[1] << ", Total signal = " << stationVector[i].GetTotalSignal() << endl;	// DEBUG
	    }
	 }
      }
   }
   
   fFile->Close();
   delete fFile;
}*/

int main(int argc, char **argv)
{
   gSystem->Load("libTree.so");

   string *stemp = new string[3];
   int *itemp = new int[4];
   double *dtemp = new double[7];
   bool goodrec;
   int printEvent = -1;

   string inname;

   RecEventFile *fFile;
   RecEvent *fRecEvent;
   DetectorGeometry *fDetGeo;
   SdRecShower *sdrecshw;
   int nevents;

   RootStyle *mystyle = new RootStyle();

   if(argc > 1)
   {
      mystyle->SetBaseStyle();
      TCanvas *c1 = new TCanvas("c1", "", 1200, 900);

      TGraph *gr = new TGraph();
      TH1F *hist = new TH1F("hist", "", 60, 0., 30.);

      itemp[0] = 0;

      for(int i = 0; i < argc-1; i++)
      {
         inname = argv[i+1];

         // Open and prepare the ADST files for reading
         fFile = new RecEventFile(inname.c_str(), RecEventFile::eRead);
         fRecEvent = new RecEvent();
         fDetGeo = new DetectorGeometry();
         sdrecshw = new SdRecShower();
         nevents = fFile->GetNEvents();
         cout << "Number of events: " << nevents << endl;

         fFile->SetBuffers(&fRecEvent);
         fFile->ReadDetectorGeometry(*fDetGeo);

         for(int j = 0; j < nevents; j++)
         {
            fFile->ReadEvent(j);

	    dtemp[0] = TMath::Log10(fRecEvent->GetSDEvent().GetSdRecShower().GetEnergy());
	    dtemp[1] = SecTheta(fRecEvent->GetSDEvent().GetSdRecShower().GetZenith(),false);

	    if((dtemp[0] >= 19.4) && (dtemp[0] < 19.5))
	    {
	       if((dtemp[1] >= 1.) && (dtemp[1] < 2.))
	       {
                  cout << "New event (" << j+1 << ", ID = " << fRecEvent->GetEventId() << ", Time = [" << fRecEvent->GetYYMMDD() << "," << fRecEvent->GetHHMMSS() << "], E = " << TMath::Log10(fRecEvent->GetSDEvent().GetSdRecShower().GetEnergy()) << ", sec(theta) = " << SecTheta(fRecEvent->GetSDEvent().GetSdRecShower().GetZenith(),false) << ") -------" << endl;

                  // Energy SD error
//	          dtemp[1] = fRecEvent->GetSDEvent().GetSdRecShower().GetEnergyError()/1.e+18;
                  // Energy SD total error
//	          dtemp[1] = fRecEvent->GetSDEvent().GetSdRecShower().GetEnergyTotalError()/1.e+18;
                  // Statistic S1000 error
//	          dtemp[1] = fRecEvent->GetSDEvent().GetSdRecShower().GetShowerSizeError();
	          // Systematic S1000 error
	          dtemp[1] = fRecEvent->GetSDEvent().GetSdRecShower().GetShowerSizeSys();
//
	          gr->SetPoint(itemp[0], dtemp[0], dtemp[1]);
	          hist->Fill(dtemp[1]);
	          itemp[0]++;
	       }
	    }
	 }
   
         fFile->Close();
         delete fRecEvent;
         delete fDetGeo;
         delete sdrecshw;
         delete fFile;
      }

      gr->SetMarkerColor(2);
      gr->SetMarkerSize(0.6);
      gr->SetMarkerStyle(20);
      gr->Draw("AP");

      cout << "X-axis mean = " << gr->GetMean(1) << endl;
      cout << "Y-axis mean = " << gr->GetMean(2) << endl;

      c1->SaveAs("output_graph.pdf");

      hist->SetLineColor(2);
      hist->SetLineWidth(2);
      hist->Draw();

      c1->SaveAs("output_distribution.pdf");

      delete hist;
      delete gr;
      delete c1;
      delete mystyle;
   }
   else
   {
      cerr << "Error! No input files supplied. Rerun program and add input files as arguments (ADST files)." << endl;
      return 1;
   }

   delete[] stemp;
   delete[] itemp;
   delete[] dtemp;

   return 0;
}
