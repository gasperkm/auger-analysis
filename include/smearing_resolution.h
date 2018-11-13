#ifndef _Resolution_h_
#define _Resolution_h_
#include "smearing_selection.h"
#include <iostream>
#include <cmath>

namespace erAna {
  namespace reso {

    enum EResoSource {
      eMolAtm = 0,   //
      eStatVAOD,     //
      eCorrVAOD,      //
      eVAODUniformity,     //
      eMultipleScattering, //
      eAlignment,          //
      eDetector,
      eNReso,
      eFirst = 0,
      eLast = eNReso - 1,
      eMolFirst = eMolAtm,
      eMolLast = eMolAtm,
      eVAODFirst = eStatVAOD,
      eVAODLast = eMultipleScattering,
      eDetFirst = eAlignment,
      eDetLast = eDetector,
      eNonMCFirst = eMolAtm,
      eNonMCLast = eAlignment
    };

    enum EResoSys {
      eVariance,
      eVarUp,
      eVarLow
    };

    enum ESigma {
      eMinus = -1,
      eZero = 0,
      ePlus = +1
    };

    inline double XmaxMolAtmVariance(const double lgE, const EResoSys sys)
    {
      // Eq. (E.5)
      if (sys == eVariance)
        return pow(2 + 0.75 * (lgE-18), 2);
      else
        return 0;
    }

    inline double XmaxStatVAODVariance(const double lgE, const EResoSys sys)
    {
      // Eq. (E.6)
      if (sys == eVariance)
        return 12. / (exp((17.9-lgE)/0.28)+1);
      else
        return 0;
    }

    inline double XmaxCorrVAODVariance(const double lgE, const EResoSys /*sys*/)
    {
      // Eq. (E.11) and (E.12)
      const double sigmaHalf = 2.7 / (exp((17.4-lgE)/0.6)+1) / 2;

      // for all three cases: sigma = sigmaHalf +/- sigmaHalf
      return pow(sigmaHalf, 2);

    }

    inline double XmaxVAODUniformityVariance(const double lgE, const EResoSys sys)
    {
      // Eq. (E.8)
      const double varLAUnif = pow(14. / (exp((17.8-lgE)/0.65)+1), 2);

      // Eq. (E.9)
      const double varVAODStat =  XmaxStatVAODVariance(lgE, sys);
      const double varUnif = (varLAUnif - 2 * varVAODStat) / 2;
      const double sigmaUnif = sqrt(varUnif);

      // Eq. (E.10): +0 - sigma/2
      if (sys == eVariance)
        return pow(0.75 * sigmaUnif, 2);
      else
        return pow(0.25 * sigmaUnif, 2);
    }

    inline double XmaxAlignmentVariance(const double lgE, const EResoSys /*sys*/)
    {
      // Eq. (E.16) and (E.17)
      const double sigmaHalf = (5 + 1.1 * (lgE - 18)) / 2;

      // for all three cases: sigma = sigmaHalf +/- sigmaHalf
      return pow(sigmaHalf, 2);
    }

    inline double XmaxMultScattVariance(const double /*lgE*/, const EResoSys sys)
    {
      // Eq. (E.22)
      if (sys == eVariance)
        return 1*1;
      else
        return 0;
    }

    inline double XmaxDetectorWidthRatio(const double lgE, const ESelection selType)
    {
      // Eq. (E.25)
      const double z18 = lgE-18;
      switch (selType)
      {
      case eDefault:
      case eDefaultProton:
      case eDefaultIron:
      case eDefaultMixed:
        return 2.1 + 0.71 * z18;
      case eNoFOV1:
      case eNoFOV1Up:
      case eNoFOV1Low:
        return 2.05302 + 7.34562e-01 * z18;
      case eNoFOV2:
        return 2.26589 + 8.21584e-01 * z18;
      case eVertical:
        return 2.13629e+00 + 5.63503e-01 * z18;
      case eInclined:
        return 2.10425e+00 + 6.98234e-01 * z18;
      default:
        std::cerr << " unknown selection type " << selType << std::endl;
        return 0;
      }
    }

    inline double XmaxDetectorGaussFraction(const double lgE, const ESelection selType)
    {
      // Eq. (E.27)
      const double z18 = lgE-18;
      switch (selType)
      {
      case eDefault:
      case eDefaultProton:
      case eDefaultIron:
      case eDefaultMixed:
        return 0.63 + 0.088 * z18;
      case eNoFOV1:
      case eNoFOV1Up:
      case eNoFOV1Low:
        return 5.12966e-01 + 1.05651e-01 * z18;
      case eNoFOV2:
        return 5.34273e-01 + 6.16180e-02 * z18;
      case eVertical:
        return 6.95688e-01 - 2.59457e-02 * z18;
      case eInclined:
        return 5.96234e-01 + 1.38184e-01 * z18;
      default:
        std::cerr << " unknown selection type " << selType << std::endl;
        return 0;
      }
    }

    inline double detFunc(const double z18, const double* p)
    {
      return pow(p[0],2) + pow(p[1]*exp(-z18/(p[2]-p[3]*z18)), 2);
    }

    inline double XmaxDetectorVariance(const double lgE, const EResoSys sys, const ESelection selType)
    {
      // Eq. (E.26)
      double pDefault[4] = { 12, 19.415, 1.537, 0.621 };
      double pDefaultProton[4] = { 6.93778e-03, 2.24665e+01, 1.87603e+00, -8.40818e-01 };
      double pDefaultIron[4] = { 9.68743e+00, 1.90939e+01, 1.62752e+00, 7.75500e-01 };
      double pDefaultMixed[4] = { 1.16179e+01, 1.86866e+01, 1.49197e+00, 4.45729e-01 };
      double pNoFOV1[4] = { 1.36653e+01, 2.27797e+01, 2.13796e+00, 6.86524e-01 };
      double pNoFOV2[4] = { 1.76932e+01, 2.40393e+01, 2.5597, 1.18866 };
      double pVertical[4] = { 8.90458e+00, 1.86387e+01, 1.23990, -3.92643e-01 };
      double pHorizontal[4] = { 3.43324e-04, 2.47195e+01, 1.86651, 9.48260e-02 };

      double* p = NULL;
      switch (selType)
      {
      case eDefault:
        p = pDefault;
        break;
      case eDefaultProton:
        p = pDefaultProton;
        break;
      case eDefaultIron:
        p = pDefaultIron;
        break;
      case eDefaultMixed:
        p = pDefaultMixed;
        break;
      case eNoFOV1:
      case eNoFOV1Up:
      case eNoFOV1Low:
        p = pNoFOV1;
        break;
      case eNoFOV2:
        p = pNoFOV2;
        break;
      case eVertical:
        p = pVertical;
        break;
      case eInclined:
        p = pHorizontal;
        break;
      default:
        std::cerr << " unknown selection type " << selType << std::endl;
        break;
      }

      const double z18 = lgE-18;
      const double vDet = detFunc(z18, p);
      if (sys == eVariance)
        return vDet;
      else if (sys == eVarUp)
      {
        p[1] += 1.4;
        p[2] += 0.09;
        p[3] -=  0.17;
        const double vUp = detFunc(z18, p);
        return pow(sqrt(vUp) - sqrt(vDet), 2);
      }
      else
      {
        p[1] -= 0.9;
        p[2] -= 0.05;
        p[3] +=  0.07;
        const double vLo = detFunc(z18, p);
        return pow(sqrt(vDet) - sqrt(vLo), 2);
      }
    }

    // first sigma of detector variance double-Gauss
    inline double XmaxDetectorVariance1(const double lgE, const EResoSys sys, const ESelection selType)
    {
      const double varianceFull = pow(sqrt(XmaxDetectorVariance(lgE, eVariance, selType)) + (sys == eVariance ? 0 : (sys == eVarUp ? 1 : -1) * sqrt(XmaxDetectorVariance(lgE, sys, selType))), 2);
      const double alpha2 = pow(XmaxDetectorWidthRatio(lgE, selType), 2);
      const double f = XmaxDetectorGaussFraction(lgE, selType);
      return varianceFull/(f + (1 - f) * alpha2);
    }

    // second sigma of detector variance double-Gauss
    inline double XmaxDetectorVariance2(const double lgE, const EResoSys sys, const ESelection selType)
    {
      const double varianceFull = pow(sqrt(XmaxDetectorVariance(lgE, eVariance, selType)) + (sys == eVariance ? 0 : (sys == eVarUp ? 1 : -1) * sqrt(XmaxDetectorVariance(lgE, sys, selType))), 2);
      const double variance1 = XmaxDetectorVariance1(lgE, sys, selType);
      const double f = XmaxDetectorGaussFraction(lgE, selType);
      return (varianceFull - f * variance1)/(1 - f);
    }

    inline double XmaxVarianceInternal(const double lgE, const EResoSource source, const ESelection selType, const EResoSys sys)
    {
      switch (source)
      {
      case eMolAtm:
        return XmaxMolAtmVariance(lgE, sys);
      case eStatVAOD:
        return XmaxStatVAODVariance(lgE, sys);
      case eCorrVAOD:
        return XmaxCorrVAODVariance(lgE, sys);
      case eVAODUniformity:
        return XmaxVAODUniformityVariance(lgE, sys);
      case eAlignment:
        return XmaxAlignmentVariance(lgE, sys);
      case eMultipleScattering:
        return XmaxMultScattVariance(lgE, sys);
      case eDetector:
        return XmaxDetectorVariance(lgE, sys, selType);
      default:
        std::cerr << " this should never happen " << std::endl;
        return 0;
      }
    }

    inline double XmaxVariance(const double lgE, const EResoSource source1, const EResoSource source2, const ESelection selType, const ESigma sigma)
    {
      double varSum = 0;
      double varVarSum = 0;
      for (int i = source1; i <= source2; ++i)
      {
        EResoSource source = (EResoSource) i;
        const double centralVar = XmaxVarianceInternal(lgE, source, selType, eVariance);
        const double varVar = XmaxVarianceInternal(lgE, source, selType, sigma == ePlus ? eVarUp : eVarLow);
        varSum += centralVar;
        varVarSum += std::abs(double(sigma)) * centralVar * varVar;

      }
      const double varianceOfVariance = varVarSum/varSum;
      return pow(sqrt(varSum) + sigma * sqrt(varianceOfVariance), 2);
    }

    inline double XmaxVariance(const double lgE, const ESelection selType, const ESigma sigma, const bool atmosphere)
    {
      if (atmosphere)
        return XmaxVariance(lgE, eFirst, eLast, selType, sigma);
      else
        return XmaxVariance(lgE, eDetFirst, eDetLast, selType, sigma);
    }

    // first sigma of full double-Gauss
    inline double XmaxVariance1(const double lgE, const ESigma nSigma, const ESelection selType)
    {
      // full variance including syst up/down
      const double xmaxVarTot = XmaxVariance(lgE, eFirst, eLast, selType, nSigma);

      // default variances (no syst)
      const double xmaxVarNonDet = XmaxVariance(lgE, eNonMCFirst, eNonMCLast, selType, eZero);
      const double xmaxVarDet1 = XmaxDetectorVariance1(lgE, eVariance, selType);
      const double xmaxVarDet2 = XmaxDetectorVariance2(lgE, eVariance, selType);
      const double fraction = XmaxDetectorGaussFraction(lgE, selType);
      const double xmaxVar1NoSys = xmaxVarDet1 + xmaxVarNonDet;
      const double xmaxVar2NoSys = xmaxVarDet2 + xmaxVarNonDet;
      const double xmaxVarTotNoSys = fraction * xmaxVar1NoSys + (1- fraction) * xmaxVar2NoSys;
      // rescaling factor to obtain full width
      const double factor = xmaxVarTot/xmaxVarTotNoSys;
      return xmaxVar1NoSys * factor;
    }

    // second sigma of full double-Gauss
    inline double XmaxVariance2(const double lgE, const ESigma nSigma, const ESelection selType)
    {
      // full variance including syst up/down
      const double xmaxVarTot = XmaxVariance(lgE, eFirst, eLast, selType, nSigma);

      // default variances (no syst)
      const double xmaxVarNonDet = XmaxVariance(lgE, eNonMCFirst, eNonMCLast, selType, eZero);
      const double xmaxVarDet1 = XmaxDetectorVariance1(lgE, eVariance, selType);
      const double xmaxVarDet2 = XmaxDetectorVariance2(lgE, eVariance, selType);
      const double fraction = XmaxDetectorGaussFraction(lgE, selType);
      const double xmaxVar1NoSys = xmaxVarDet1 + xmaxVarNonDet;
      const double xmaxVar2NoSys = xmaxVarDet2 + xmaxVarNonDet;
      const double xmaxVarTotNoSys = fraction * xmaxVar1NoSys + (1- fraction) * xmaxVar2NoSys;
      // rescaling factor to obtain full width
      const double factor = xmaxVarTot/xmaxVarTotNoSys;
      return xmaxVar2NoSys * factor;
    }

    // Xmax^\prime = Xmax + XmaxCorrection
    inline double XmaxCorrection(const double lgE, const ESelection selType, const bool isMC = false)
    {
      // Eq. (E.24)
      const double z18 = lgE-18;
      double mu = 0;
      switch (selType)
      {
      case eDefault:
        mu = -3.4 + 0.93  * z18;
        break;
      case eDefaultMixed:
        mu = -3.38264e+00 + 1.59347e+00  * z18;
        break;
      case eDefaultIron:
        mu = -3.19945e+00 + 1.02833e+00 * z18;
        break;
      case eDefaultProton:
        mu = -3.59148e+00 + 2.17117e+00 * z18;
        break;
      case eNoFOV1:
      case eNoFOV1Up:
      case eNoFOV1Low:
        mu = -3.1 + 0.98  * z18;
        break;
      case eNoFOV2:
        mu = -2.7 + 0.82  * z18;
        break;
      case eVertical:
        mu = -2.49675 + 1.057 * z18;
        break;
      case eInclined:
        mu = -4.26646 + 9.71897e-01 * z18;
        break;
      default:
        std::cerr << " unknown selection type " << selType << std::endl;
        mu = 0;
      }

      // Eq.(E.31)
      const double lwCorrBias = isMC ? 0 : 6.5 / (exp((lgE-18.23)/0.41)+1);

      return -lwCorrBias - mu;
    }
  }
}

#endif
