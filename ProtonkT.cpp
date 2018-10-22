/*
 * ProtonkT.cpp
 *
 *  Created on: Sep 10, 2018
 *      Author: hohlweger
 */

#include "ProtonkT.h"

#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TString.h"
#include <vector>
#include "TLatex.h"
#include "TStyle.h"
#include "TGraphErrors.h"
void FitkTBins(const char* filename, const char* OutputDir, int DataSet) {
  std::vector<int> fFillColors;
  fFillColors.push_back(kGray + 1);
  fFillColors.push_back(kRed - 10);
  fFillColors.push_back(kBlue - 9);
  fFillColors.push_back(kGreen - 8);
  fFillColors.push_back(kMagenta - 9);
  fFillColors.push_back(kOrange - 9);
  fFillColors.push_back(kCyan - 3);
  fFillColors.push_back(kYellow - 7);

  std::vector<int> fColors;
  fColors.push_back(kBlack);
  fColors.push_back(kRed + 1);
  fColors.push_back(kBlue + 2);
  fColors.push_back(kGreen + 3);
  fColors.push_back(kMagenta + 1);
  fColors.push_back(kOrange - 1);
  fColors.push_back(kCyan + 2);
  fColors.push_back(kYellow + 2);

  std::vector<int> fMarkers;
  fMarkers.push_back(kFullCircle);
  fMarkers.push_back(kFullSquare);
  fMarkers.push_back(kOpenCircle);
  fMarkers.push_back(kOpenSquare);
  fMarkers.push_back(kOpenDiamond);
  fMarkers.push_back(kOpenCross);
  fMarkers.push_back(kFullCross);
  fMarkers.push_back(kFullDiamond);
  fMarkers.push_back(kFullStar);
  fMarkers.push_back(kOpenStar);

  const unsigned NumMomBins_pp = 50;
  const double kMin_pp = 4;
  const double kMax_pp = kMin_pp + 8 * NumMomBins_pp;
  unsigned NumkTBins = 4;
  const float right = 0.025;
  const float top = 0.025;
  const float left = 0.15;

  TCanvas** cankTBins = new TCanvas*[NumkTBins];
  TCanvas* canRadius = new TCanvas("canRadius", "canRadius", 0, 0, 350, 500);

  canRadius->SetRightMargin(right);
  canRadius->SetLeftMargin(left);
  canRadius->SetTopMargin(top);

  TH1F** CFMeV_kT = new TH1F*[NumkTBins];

  TFile* input = TFile::Open(filename, "READ");
  for (int ikT = 0; ikT < NumkTBins; ++ikT) {
//    TString CFkTHistName = Form("hCkTotNormWeightMeV_kTBin_%i", ikT);
    TString CFkTHistName = Form("hCk_RebinnedMeV_0_kTBin_%i", ikT);
    CFMeV_kT[ikT] = (TH1F*) input->Get(CFkTHistName.Data());
    if (!CFMeV_kT[ikT]) {
      std::cout << CFkTHistName.Data() << " missing \n";
      return;
    }
  }
  TGraphErrors *meankT = (TGraphErrors*)input->Get("AveragekT");
  if (!meankT) {
    std::cout << "No mean kT" << std::endl;
    return;
  }
  TString CalibBaseDir = "";
  if (DataSet == 0) {  // pPb MB
    CalibBaseDir +=
        "/home/gu74req/Analysis/CATS_Input/SystematicsAndCalib/pPbRun2_MB";
  } else if (DataSet == 1) {  // pp MB
    CalibBaseDir +=
        "/home/gu74req/Analysis/CATS_Input/SystematicsAndCalib/ppRun2_MB";
  } else if (DataSet == 2) {  // pp HM
    CalibBaseDir +=
        "/home/gu74req/Analysis/CATS_Input/SystematicsAndCalib/ppRun2_HM";
  }
  TString ResMatrixFileName = TString::Format("%s/run2_decay_matrices_old.root",
                                              CalibBaseDir.Data());
  TString ResMatrixHistName = "hRes_pp_pL";
  TString SigmaMatrixFileName = TString::Format("%s/Sample3_MeV_compact.root",
                                                CalibBaseDir.Data());
  TString SigmaMatrixHistName = "hSigmaMeV_Proton_Proton";

  int Fraction_Res = 1;
  int Fraction_Sig = 1;
  double UnitConv_Res = 1;
  double UnitConv_Sig = 1;

  TFile* GraphFile = new TFile(
      TString::Format("%sGraphFileMULT.root", OutputDir), "recreate");
  TGraph* FitResult_pp = new TGraph[NumkTBins];
  TGraphErrors* RadiusResults = new TGraphErrors();
  TH2F* hRes_pp_pL;
  TH2F* hSigma_pp;

  TFile* FileRes = new TFile(ResMatrixFileName, "read");
  TFile* FileSigma = new TFile(SigmaMatrixFileName, "read");

  //for the Xi we will make the assumption that all residuals are flat since we do not know better
  FileRes->cd();
  hRes_pp_pL = (TH2F*) FileRes->Get(ResMatrixHistName);

  FileSigma->cd();
  hSigma_pp = (TH2F*) FileSigma->Get(SigmaMatrixHistName);

  printf("hSigma_pp=%p\n", hSigma_pp);

  TH2F* hRes_pp_pL_MeV = new TH2F(
      "hRes_pp_pL_MeV",
      "hRes_pp_pL_MeV",
      hRes_pp_pL->GetNbinsX() / Fraction_Res,
      hRes_pp_pL->GetXaxis()->GetBinLowEdge(1) * UnitConv_Res,
      hRes_pp_pL->GetXaxis()->GetBinUpEdge(
          hRes_pp_pL->GetNbinsX() / Fraction_Res) * UnitConv_Res,
      hRes_pp_pL->GetNbinsY() / Fraction_Res,
      hRes_pp_pL->GetYaxis()->GetBinLowEdge(1) * UnitConv_Res,
      hRes_pp_pL->GetXaxis()->GetBinUpEdge(
          hRes_pp_pL->GetNbinsY() / Fraction_Res) * UnitConv_Res);
  TH2F* hSigma_pp_MeV = new TH2F(
      "hSigma_pp_MeV",
      "hSigma_pp_MeV",
      hSigma_pp->GetNbinsX() / Fraction_Sig,
      hSigma_pp->GetXaxis()->GetBinLowEdge(1) * UnitConv_Sig,
      hSigma_pp->GetXaxis()->GetBinUpEdge(hSigma_pp->GetNbinsX() / Fraction_Sig)
          * UnitConv_Sig,
      hSigma_pp->GetNbinsY() / Fraction_Sig,
      hSigma_pp->GetYaxis()->GetBinLowEdge(1) * UnitConv_Sig,
      hSigma_pp->GetXaxis()->GetBinUpEdge(hSigma_pp->GetNbinsY() / Fraction_Sig)
          * UnitConv_Sig);

  for (int iBinX = 1; iBinX <= hRes_pp_pL->GetNbinsX() / Fraction_Res;
      iBinX++) {
    for (int iBinY = 1; iBinY <= hRes_pp_pL->GetNbinsY() / Fraction_Res;
        iBinY++) {
      hRes_pp_pL_MeV->SetBinContent(iBinX, iBinY,
                                    hRes_pp_pL->GetBinContent(iBinX, iBinY));
    }
  }
  for (int iBinX = 1; iBinX <= hSigma_pp->GetNbinsX() / Fraction_Sig; iBinX++) {
    for (int iBinY = 1; iBinY <= hSigma_pp->GetNbinsY() / Fraction_Sig;
        iBinY++) {
      hSigma_pp_MeV->SetBinContent(iBinX, iBinY,
                                   hSigma_pp->GetBinContent(iBinX, iBinY));
    }
  }

  int vSource;
  int vFit;

  TRandom3 rangen(1);

  vSource = 0;  //0 Gauss, 1 Cauchy, 2 DoubleGauss, 3 EPOS
  vFit = 1;  //1 = with BL

  double FemtoRegion_pp[2];
  FemtoRegion_pp[0] = 4;
  FemtoRegion_pp[1] = 200;

  double BlRegion[2];
  BlRegion[0] = 300;
  BlRegion[1] = 500;

  const double GaussSourceSize = 1.2;
  const double Mass_p = 938.272;

  double PotPars1S0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0, 0 };
  double PotPars3P0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 0 };
  double PotPars3P1[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 1 };
  double PotPars3P2[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 2 };

  const double Weight1S0 = 3. / 12.;
  const double Weight3P0 = 1. / 12.;
  const double Weight3P1 = 3. / 12.;
  const double Weight3P2 = 5. / 12.;

  CATS AB_pp[NumkTBins];
  double Pars_pp[NumkTBins][6];

  for (unsigned ikT = 0; ikT < NumkTBins; ikT++) {
    Pars_pp[ikT][0] = 0;
    Pars_pp[ikT][1] = 0;
    Pars_pp[ikT][2] = 0;
    Pars_pp[ikT][3] = GaussSourceSize * 1.2;
    Pars_pp[ikT][4] = GaussSourceSize / 1.2;
    Pars_pp[ikT][5] = 0.5;

    if (vSource == 0) {
      AB_pp[ikT].SetAnaSource(GaussSource, Pars_pp[ikT]);
      AB_pp[ikT].SetUseAnalyticSource(true);
      AB_pp[ikT].SetThetaDependentSource(false);
    } else if (vSource == 1) {
      AB_pp[ikT].SetAnaSource(CauchySource, Pars_pp[ikT]);
      AB_pp[ikT].SetUseAnalyticSource(true);
      AB_pp[ikT].SetThetaDependentSource(false);
    } else if (vSource == 2) {
      AB_pp[ikT].SetAnaSource(DoubleGaussSource, Pars_pp[ikT]);
      AB_pp[ikT].SetUseAnalyticSource(true);
      AB_pp[ikT].SetThetaDependentSource(false);
    }

    AB_pp[ikT].SetExcludeFailedBins(false);
    AB_pp[ikT].SetMomBins(NumMomBins_pp, kMin_pp, kMax_pp);

    AB_pp[ikT].SetNumChannels(4);
    AB_pp[ikT].SetNumPW(0, 2);
    AB_pp[ikT].SetNumPW(1, 2);
    AB_pp[ikT].SetNumPW(2, 2);
    AB_pp[ikT].SetNumPW(3, 2);
    AB_pp[ikT].SetSpin(0, 0);
    AB_pp[ikT].SetSpin(1, 1);
    AB_pp[ikT].SetSpin(2, 1);
    AB_pp[ikT].SetSpin(3, 1);
    AB_pp[ikT].SetChannelWeight(0, Weight1S0);
    AB_pp[ikT].SetChannelWeight(1, Weight3P0);
    AB_pp[ikT].SetChannelWeight(2, Weight3P1);
    AB_pp[ikT].SetChannelWeight(3, Weight3P2);

    AB_pp[ikT].SetQ1Q2(1);
    AB_pp[ikT].SetPdgId(2212, 2212);
    AB_pp[ikT].SetRedMass(0.5 * Mass_p);

    AB_pp[ikT].SetShortRangePotential(0, 0, fDlmPot, PotPars1S0);
    AB_pp[ikT].SetShortRangePotential(1, 1, fDlmPot, PotPars3P0);
    AB_pp[ikT].SetShortRangePotential(2, 1, fDlmPot, PotPars3P1);
    AB_pp[ikT].SetShortRangePotential(3, 1, fDlmPot, PotPars3P2);

    AB_pp[ikT].KillTheCat();
  }

  //!IMPORTANT SETTINGS
  const unsigned NumSourcePars = vSource == 0 ? 1 : vSource == 1 ? 1 :
                                 vSource == 2 ? 3 : 1;

  DLM_Ck** Ck_pp = new DLM_Ck*[NumkTBins];

  DLM_CkDecomposition** CkDec_pp = new DLM_CkDecomposition*[NumkTBins];

  double lam_pp;
  double lam_pp_pL;
  double lam_pp_fake;

  double lam_pL;
  double lam_pL_Feed;
  double lam_pL_fake;

  if (DataSet == 0) {
    lam_pp = 0.720;
    lam_pp_pL = 0.161;
    lam_pp_fake = 1. - lam_pp - lam_pp_pL;

    lam_pL = 0.589;
    lam_pL_fake = 0.083;
    lam_pL_Feed = 1. - lam_pL - lam_pL_fake;
  } else if (DataSet == 1) {
    lam_pp = 0.752;
    lam_pp_pL = 0.151;
    lam_pp_fake = 1. - lam_pp - lam_pp_pL;

    lam_pL = 0.519;
    lam_pL_fake = 0.043;
    lam_pL_Feed = 1. - lam_pL - lam_pL_fake;
  }
  DLM_Fitter1** fTest1 = new DLM_Fitter1*[NumkTBins];
  double CkVal;
  double Momentum;
  double ExpRad = 0;
  for (unsigned ikT = 0; ikT < NumkTBins; ikT++) {
    Ck_pp[ikT] = new DLM_Ck(NumSourcePars, 0, AB_pp[ikT]);
    Ck_pp[ikT]->Update();
    CkDec_pp[ikT] = new DLM_CkDecomposition("pp", 2, *Ck_pp[ikT],
                                            hSigma_pp_MeV);
    CkDec_pp[ikT]->AddContribution(0, 1. - lam_pp - lam_pp_fake,
                                   DLM_CkDecomposition::cFeedDown);
    CkDec_pp[ikT]->AddContribution(1, lam_pp_fake, DLM_CkDecomposition::cFake);  //0.02
    fTest1[ikT] = new DLM_Fitter1(1);
    fTest1[ikT]->SetOutputDir(OutputDir);

    if (vFit == 0) {
      fTest1[ikT]->SetSystem(0, *CFMeV_kT[ikT], 1, *CkDec_pp[ikT],
                             FemtoRegion_pp[0], FemtoRegion_pp[1],
                             FemtoRegion_pp[1], FemtoRegion_pp[1]);
    } else {
      fTest1[ikT]->SetSystem(0, *CFMeV_kT[ikT], 1, *CkDec_pp[ikT],
                             FemtoRegion_pp[0], FemtoRegion_pp[1], BlRegion[0],
                             BlRegion[1]);
    }

    fTest1[ikT]->SetParameter("pp", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
    if (vFit < 1)
      fTest1[ikT]->FixParameter("pp", DLM_Fitter1::p_b, 0);
    else
      fTest1[ikT]->SetParameter("pp", DLM_Fitter1::p_b, 1e-4, 0, 2e-3);
    if (vFit < 2)
      fTest1[ikT]->FixParameter("pp", DLM_Fitter1::p_c, 0);
    else
      fTest1[ikT]->SetParameter("pp", DLM_Fitter1::p_c, 0, -2e-4, 2e-4);
    if (vSource == 0) {
      fTest1[ikT]->SetParameter("pp", DLM_Fitter1::p_sor0, GaussSourceSize, 0.,
                                2.0);
      //fTest1->FixParameter("pp",DLM_Fitter1::p_sor0,1.2);
    } else if (vSource == 1) {
      fTest1[ikT]->SetParameter("pp", DLM_Fitter1::p_sor0,
                                GaussSourceSize / 1.4, 0.5, 1.6);
    } else if (vSource == 2) {
      fTest1[ikT]->SetParameter("pp", DLM_Fitter1::p_sor0, 0.8, 0.4, 3.6);
      fTest1[ikT]->SetParameter("pp", DLM_Fitter1::p_sor1, 2.4, 0.4, 3.6);
      fTest1[ikT]->SetParameter("pp", DLM_Fitter1::p_sor2, 0.5, 0, 1);
    }

    fTest1[ikT]->SetParameter("pp", DLM_Fitter1::p_Cl, -0.9, -1.2, -0.8);
    fTest1[ikT]->SetSeparateBL(0, false);

    CkDec_pp[ikT]->Update();
    fTest1[ikT]->GoBabyGo();

    printf("For ikT=%u I got r = %.3f fm \n", ikT,
           fTest1[ikT]->GetParameter("pp", DLM_Fitter1::p_sor0));

    FitResult_pp[ikT].SetName(TString::Format("FitResult_pp_%u", ikT));
    fTest1[ikT]->GetFitGraph(0, FitResult_pp[ikT]);
    FitResult_pp[ikT].SetMarkerColor(fColors[ikT]);
    FitResult_pp[ikT].SetMarkerStyle(20);
    FitResult_pp[ikT].SetMarkerSize(0);
    FitResult_pp[ikT].SetLineColor(fColors[ikT]);
    FitResult_pp[ikT].SetLineWidth(2);
    FitResult_pp[ikT].SetFillColor(kWhite);

    GraphFile->cd();
    FitResult_pp[ikT].Write();
    TString canName = Form("c%i", ikT);
    cankTBins[ikT] = new TCanvas(canName.Data(), canName.Data(), 0, 0, 650,
                                 550);
    cankTBins[ikT]->SetRightMargin(right);
    cankTBins[ikT]->SetTopMargin(top);
    cankTBins[ikT]->cd();
    CFMeV_kT[ikT]->GetXaxis()->SetLabelSize(0.045);
    CFMeV_kT[ikT]->GetXaxis()->SetTitleSize(0.05);
    CFMeV_kT[ikT]->GetXaxis()->SetLabelOffset(0.01);
    CFMeV_kT[ikT]->GetXaxis()->SetTitleOffset(1.);
    CFMeV_kT[ikT]->GetXaxis()->SetLabelFont(42);
    CFMeV_kT[ikT]->GetYaxis()->SetLabelSize(0.045);
    CFMeV_kT[ikT]->GetYaxis()->SetTitleSize(0.05);
    CFMeV_kT[ikT]->GetYaxis()->SetLabelOffset(0.01);
    CFMeV_kT[ikT]->GetYaxis()->SetTitleOffset(1.);
    CFMeV_kT[ikT]->SetMarkerSize(1.5);
    CFMeV_kT[ikT]->SetLineWidth(2);
    CFMeV_kT[ikT]->GetXaxis()->SetRangeUser(0, 200);
    if (ikT == 0) {
      CFMeV_kT[ikT]->GetYaxis()->SetRangeUser(0.9, 2.8);
    } else if (ikT == 1) {
      CFMeV_kT[ikT]->GetYaxis()->SetRangeUser(0.8, 2.8);
    } else if (ikT == 2) {
      CFMeV_kT[ikT]->GetYaxis()->SetRangeUser(0.3, 3.3);
    } else if (ikT == 3) {
      CFMeV_kT[ikT]->GetYaxis()->SetRangeUser(0.4, 4.2);
    }
    CFMeV_kT[ikT]->SetStats(0);
    CFMeV_kT[ikT]->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
    CFMeV_kT[ikT]->Draw();
    FitResult_pp[ikT].Draw("L3 same");
    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize() * 0.85);
    BeamText.SetNDC(kTRUE);
    BeamText.DrawLatex(
        0.3,
        0.875,
        Form("Radius for ikT = %i : ( %.3f #pm %.3f (stat.)) fm", ikT,
             fTest1[ikT]->GetParameter("pp", DLM_Fitter1::p_sor0),
             fTest1[ikT]->GetParError("pp", DLM_Fitter1::p_sor0)));
    cankTBins[ikT]->SaveAs(Form("%sppFit_ikT_%i.png", OutputDir, ikT));
    cankTBins[ikT]->SaveAs(Form("%sppFit_ikT_%i.pdf", OutputDir, ikT));
    cankTBins[ikT]->Write();
    double dummy = 0;
    double theMeankT = 0;
    double theMeankTErr = 0;
    meankT->GetPoint(ikT, dummy, theMeankT);
    theMeankTErr = meankT->GetErrorY(ikT);
    RadiusResults->SetPoint(
        ikT, theMeankT, fTest1[ikT]->GetParameter("pp", DLM_Fitter1::p_sor0));
    RadiusResults->SetPointError(
        ikT, theMeankTErr, fTest1[ikT]->GetParError("pp", DLM_Fitter1::p_sor0));
  }
  canRadius->cd();
  RadiusResults->SetMarkerStyle(fMarkers[2]);
  RadiusResults->SetMarkerSize(1.1);
  RadiusResults->SetMarkerColor(fColors[2]);
  RadiusResults->GetXaxis()->SetTitle("<k_{T}> (GeV/c)");
//  RadiusResults->GetXaxis()->SetRangeUser(-0.5, 3.5);
  RadiusResults->GetYaxis()->SetTitle("Radius (fm)");
  RadiusResults->GetXaxis()->SetLabelSize(0.045);
  RadiusResults->GetXaxis()->SetTitleSize(0.05);
  RadiusResults->GetXaxis()->SetLabelOffset(0.01);
  RadiusResults->GetXaxis()->SetTitleOffset(1.);
  RadiusResults->GetXaxis()->SetLabelFont(42);
  RadiusResults->GetYaxis()->SetLabelSize(0.045);
  RadiusResults->GetYaxis()->SetTitleSize(0.05);
  RadiusResults->GetYaxis()->SetLabelOffset(0.01);
  RadiusResults->GetYaxis()->SetTitleOffset(1.3);
  RadiusResults->Draw("AP");
  canRadius->SaveAs(Form("%skTRadius.png", OutputDir));
  canRadius->SaveAs(Form("%skTRadius.pdf", OutputDir));
  canRadius->Write();

}
