vector<string> methods;
methods.push_back()
methods.push_back("Cuts");
methods.push_back("CutsD");
methods.push_back("CutsPCA");
methods.push_back("CutsGA");
methods.push_back("CutsSA");
methods.push_back("Likelihood");
methods.push_back("LikelihoodD");
methods.push_back("LikelihoodPCA");
methods.push_back("LikelihoodKDE");
methods.push_back("LikelihoodMIX");
methods.push_back("PDERS");
methods.push_back("PDERSD");
methods.push_back("PDERSPCA");
methods.push_back("PDEFoam");
methods.push_back("PDEFoamBoost");
methods.push_back("KNN");
methods.push_back("LD");
methods.push_back("Fisher");
methods.push_back("FisherG");
methods.push_back("BoostedFisher");
methods.push_back("HMatrix");
methods.push_back("FDA_GA");
methods.push_back("FDA_SA");
methods.push_back("FDA_MC");
methods.push_back("FDA_MT");
methods.push_back("FDA_GAMT");
methods.push_back("FDA_MCMT");
methods.push_back("MLP");
methods.push_back("MLPBFGS");
methods.push_back("MLPBNN");
methods.push_back("CFMlpANN");
methods.push_back("TMlpANN");
methods.push_back("SVM");
methods.push_back("BDT");
methods.push_back("BDTG");
methods.push_back("BDTB");
methods.push_back("BDTD");
methods.push_back("BDTF");
methods.push_back("RuleFit");
int nrmethods = methods.size();

void BookTheMethod(int iMethod, TMVA::Factory *factory)
{
   for(int i = 0; i < nrmethods; i++)
   {
      if(methods[iMethod] == "Cuts")
         factory->BookMethod(TMVA::Types::kCuts, methods[iMethod], "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart");

      if(methods[iMethod] == "CutsD")
         factory->BookMethod(TMVA::Types::kCuts, methods[iMethod], "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate");

      if(methods[iMethod] == "CutsPCA")
         factory->BookMethod(TMVA::Types::kCuts, methods[iMethod], "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA");

      if(methods[iMethod] == "CutsGA")
         factory->BookMethod(TMVA::Types::kCuts, methods[iMethod], "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95");

      if(methods[iMethod] == "CutsSA")
         factory->BookMethod(TMVA::Types::kCuts, methods[iMethod], "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");

      if(methods[iMethod] == "Likelihood")
         factory->BookMethod(TMVA::Types::kLikelihood, methods[iMethod], "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50");

      if(methods[iMethod] == "LikelihoodD")
         factory->BookMethod(TMVA::Types::kLikelihood, methods[iMethod], "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate");

      if(methods[iMethod] == "LikelihoodPCA")
         factory->BookMethod(TMVA::Types::kLikelihood, methods[iMethod], "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA"); 

      if(methods[iMethod] == "LikelihoodKDE")
         factory->BookMethod(TMVA::Types::kLikelihood, methods[iMethod], "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50"); 

      if(methods[iMethod] == "LikelihoodMIX")
         factory->BookMethod(TMVA::Types::kLikelihood, methods[iMethod], "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50"); 

      if(methods[iMethod] == "PDERS")
         factory->BookMethod(TMVA::Types::kPDERS, methods[iMethod], "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600");

      if(methods[iMethod] == "PDERSD")
         factory->BookMethod(TMVA::Types::kPDERS, methods[iMethod], "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate");

      if(methods[iMethod] == "PDERSPCA")
         factory->BookMethod(TMVA::Types::kPDERS, methods[iMethod], "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA");

      if(methods[iMethod] == "PDEFoam")
         factory->BookMethod(TMVA::Types::kPDEFoam, methods[iMethod], "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T");

      if(methods[iMethod] == "PDEFoamBoost")
         factory->BookMethod(TMVA::Types::kPDEFoam, methods[iMethod], "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T");

      if(methods[iMethod] == "KNN")
         factory->BookMethod(TMVA::Types::kKNN, methods[iMethod], "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim");

      if(methods[iMethod] == "LD")
         factory->BookMethod(TMVA::Types::kLD, methods[iMethod], "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");

      if(methods[iMethod] == "Fisher")
         factory->BookMethod(TMVA::Types::kFisher, methods[iMethod], "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");

      if(methods[iMethod] == "FisherG")
         factory->BookMethod(TMVA::Types::kFisher, methods[iMethod], "H:!V:VarTransform=Gauss");

      if(methods[iMethod] == "BoostedFisher")
         factory->BookMethod(TMVA::Types::kFisher, methods[iMethod], "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring");

      if(methods[iMethod] == "HMatrix")
         factory->BookMethod(TMVA::Types::kHMatrix, methods[iMethod], "!H:!V:VarTransform=None");

      if(methods[iMethod] == "FDA_GA")
         factory->BookMethod(TMVA::Types::kFDA, methods[iMethod], "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1");

      if(methods[iMethod] == "FDA_SA")
         factory->BookMethod(TMVA::Types::kFDA, methods[iMethod], "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");

      if(methods[iMethod] == "FDA_MC")
         factory->BookMethod(TMVA::Types::kFDA, methods[iMethod], "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1");

      if(methods[iMethod] == "FDA_MT")
         factory->BookMethod(TMVA::Types::kFDA, methods[iMethod], "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch");

      if(methods[iMethod] == "FDA_GAMT")
         factory->BookMethod(TMVA::Types::kFDA, methods[iMethod], "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim");

      if(methods[iMethod] == "FDA_MCMT")
         factory->BookMethod(TMVA::Types::kFDA, methods[iMethod], "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20");

      if(methods[iMethod] == "MLP")
         factory->BookMethod(TMVA::Types::kMLP, methods[iMethod], "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator");

      if(methods[iMethod] == "MLPBFGS")
         factory->BookMethod(TMVA::Types::kMLP, methods[iMethod], "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator");

      if(methods[iMethod] == "MLPBNN")
         factory->BookMethod(TMVA::Types::kMLP, methods[iMethod], "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator");

      if(methods[iMethod] == "CFMlpANN")
         factory->BookMethod(TMVA::Types::kCFMlpANN, methods[iMethod], "!H:!V:NCycles=2000:HiddenLayers=N+1,N");

      if(methods[iMethod] == "TMlpANN")
         factory->BookMethod(TMVA::Types::kTMlpANN, methods[iMethod], "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3");

      if(methods[iMethod] == "SVM")
         factory->BookMethod(TMVA::Types::kSVM, methods[iMethod], "Gamma=0.25:Tol=0.001:VarTransform=Norm");

      if(methods[iMethod] == "BDT")
         factory->BookMethod(TMVA::Types::kBDT, methods[iMethod], "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

      if(methods[iMethod] == "BDTG")
         factory->BookMethod(TMVA::Types::kBDT, methods[iMethod], "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");

      if(methods[iMethod] == "BDTB")
         factory->BookMethod(TMVA::Types::kBDT, methods[iMethod], "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20");

      if(methods[iMethod] == "BDTD")
         factory->BookMethod(TMVA::Types::kBDT, methods[iMethod], "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate");

      if(methods[iMethod] == "BDTF")
         factory->BookMethod(TMVA::Types::kBDT, "BDTMitFisher", "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");

      if(methods[iMethod] == "RuleFit")
         factory->BookMethod(TMVA::Types::kRuleFit, methods[iMethod], "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02");
   }
}
