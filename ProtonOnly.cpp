

#include "ProtonOnly.h"
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

void PROTON_ONLY(TString InputDir,TString FileName,TString OutputDir) {

  std::vector<int> fFillColors = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-3, kYellow-7};
  std::vector<int> fColors     = {kBlack, kRed+1 , kBlue+2, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2};
  std::vector<int> fMarkers    = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};

  //ALICE_pp_13TeV
  //ALICE_pPb_5TeV
  //    const TString DataSet = "ALICE_pp_13TeV";
  const TString DataSet = "ALICE_pPb_5TeV";

  const int EPOS_MAX_PAIRS = 128e6;
  const int EPOS_MIX_DEP = 16;
  const bool EPOS_THETA_DEP = false;
  //const bool RUN_ON_NX1 = false;

  const unsigned NumMomBins_pp = 50;
  const double kMin_pp = 0;
  const double kMax_pp = kMin_pp+4*NumMomBins_pp;

  const unsigned NumMomBins_pL = 13;
  const double kMin_pL = 0;
  const double kMax_pL = kMin_pL+20*NumMomBins_pL;

  //the first one is the main
  unsigned NumMultBins = 4;
  double* NumPairs;
  //    if(DataSet=="ALICE_pp_13TeV"){
  //        NumMultBins = 4;
  //        NumPairs = new double[NumMultBins];
  //        NumPairs[0] = 28270+96648+239974;
  //        NumPairs[1] = 28270;
  //        NumPairs[2] = 96648;
  //        NumPairs[3] = 239974;
  //    }
  //    else if(DataSet=="ALICE_pPb_5TeV"){
  //        NumMultBins = 4;
  //        NumPairs = new double[NumMultBins];
  //        NumPairs[0] = 42116+156822+126164;
  //        NumPairs[1] = 42116;
  //        NumPairs[2] = 156822;
  //        NumPairs[3] = 126164;
  //    }


  TString OliFileName_pp[NumMultBins];
  TString MultHistoName[NumMultBins];
  TString SEHistoName[NumMultBins];
  double UnitsTransform[NumMultBins];


  TString ResMatrixFileName;
  TString ResMatrixHistName;
  TString SigmaMatrixFileName;
  TString SigmaMatrixHistName;

  //1/FractionOfBins th number of original bins of the correction matrix are taken into account
  //originally: 1000 MeV => 1/2 should be good most of the time
  int Fraction_Res;
  int Fraction_Sig;
  double UnitConv_Res;
  double UnitConv_Sig;


  OliFileName_pp[0] = TString::Format("%s/CFOutput_pp_Rebin_1.root",InputDir.Data());
  OliFileName_pp[2] = TString::Format("%s/%s",InputDir.Data(),FileName.Data());
  OliFileName_pp[1] = OliFileName_pp[2];
  OliFileName_pp[3] = OliFileName_pp[2];

  MultHistoName[0] = "hCkTotNormWeight";
  MultHistoName[1] = "CF_Mult_0";
  MultHistoName[2] = "CF_Mult_1";
  MultHistoName[3] = "CF_Mult_2";

  SEHistoName[0] = "ppSE";
  SEHistoName[1] = "SE_Mult_0";
  SEHistoName[2] = "SE_Mult_1";
  SEHistoName[3] = "SE_Mult_2";

  UnitsTransform[0] = 1;
  UnitsTransform[1] = 1000;
  UnitsTransform[2] = 1000;
  UnitsTransform[3] = 1000;
//
//
//  OliFileName_pp[0] = TString::Format("%s/CFOutput_pp_Rebin_1.root",InputDir.Data());
//  OliFileName_pp[1] = TString::Format("%s/%s",InputDir.Data(),FileName.Data());
//  OliFileName_pp[2] = OliFileName_pp[1];
//  OliFileName_pp[3] = OliFileName_pp[1];
//
//  MultHistoName[0] = "hCkTotNormWeight";
//  MultHistoName[1] = "CF_Mult_0";
//  MultHistoName[2] = "CF_Mult_1";
//  MultHistoName[3] = "CF_Mult_2";
//
//  SEHistoName[0] = "ppSE";
//  SEHistoName[1] = "SE_Mult_0";
//  SEHistoName[2] = "SE_Mult_1";
//  SEHistoName[3] = "SE_Mult_2";
//
//  UnitsTransform[0] = 1;
//  UnitsTransform[1] = 1000;
//  UnitsTransform[2] = 1000;
//  UnitsTransform[3] = 1000;

  ResMatrixFileName = TString::Format("%s/run2_decay_matrices_old.root",InputDir.Data());
  ResMatrixHistName = "hRes_pp_pL";
  SigmaMatrixFileName = TString::Format("%s/Sample3_MeV_compact.root",InputDir.Data());
  SigmaMatrixHistName = "hSigmaMeV_Proton_Proton";

  Fraction_Res = 1;
  Fraction_Sig = 1;
  UnitConv_Res = 1;
  UnitConv_Sig = 1;


  TH2F* hRes_pp_pL;
  TH2F* hSigma_pp;

  TFile* FileRes = new TFile(ResMatrixFileName, "read");
  TFile* FileSigma = new TFile(SigmaMatrixFileName, "read");

  //for the Xi we will make the assumption that all residuals are flat since we do not know better
  FileRes->cd();
  hRes_pp_pL = (TH2F*)FileRes->Get(ResMatrixHistName);

  FileSigma->cd();
  hSigma_pp = (TH2F*)FileSigma->Get(SigmaMatrixHistName);

  printf("hSigma_pp=%p\n",hSigma_pp);

  TH2F* hRes_pp_pL_MeV = new TH2F("hRes_pp_pL_MeV", "hRes_pp_pL_MeV",
                                  hRes_pp_pL->GetNbinsX()/Fraction_Res, hRes_pp_pL->GetXaxis()->GetBinLowEdge(1)*UnitConv_Res,hRes_pp_pL->GetXaxis()->GetBinUpEdge(hRes_pp_pL->GetNbinsX()/Fraction_Res)*UnitConv_Res,
                                  hRes_pp_pL->GetNbinsY()/Fraction_Res, hRes_pp_pL->GetYaxis()->GetBinLowEdge(1)*UnitConv_Res,hRes_pp_pL->GetXaxis()->GetBinUpEdge(hRes_pp_pL->GetNbinsY()/Fraction_Res)*UnitConv_Res);
  TH2F* hSigma_pp_MeV = new TH2F("hSigma_pp_MeV", "hSigma_pp_MeV",
                                 hSigma_pp->GetNbinsX()/Fraction_Sig, hSigma_pp->GetXaxis()->GetBinLowEdge(1)*UnitConv_Sig,hSigma_pp->GetXaxis()->GetBinUpEdge(hSigma_pp->GetNbinsX()/Fraction_Sig)*UnitConv_Sig,
                                 hSigma_pp->GetNbinsY()/Fraction_Sig, hSigma_pp->GetYaxis()->GetBinLowEdge(1)*UnitConv_Sig,hSigma_pp->GetXaxis()->GetBinUpEdge(hSigma_pp->GetNbinsY()/Fraction_Sig)*UnitConv_Sig);

  for(int iBinX=1; iBinX<=hRes_pp_pL->GetNbinsX()/Fraction_Res; iBinX++){
    for(int iBinY=1; iBinY<=hRes_pp_pL->GetNbinsY()/Fraction_Res; iBinY++){
      hRes_pp_pL_MeV->SetBinContent(iBinX, iBinY, hRes_pp_pL->GetBinContent(iBinX, iBinY));
    }
  }
  for(int iBinX=1; iBinX<=hSigma_pp->GetNbinsX()/Fraction_Sig; iBinX++){
    for(int iBinY=1; iBinY<=hSigma_pp->GetNbinsY()/Fraction_Sig; iBinY++){
      hSigma_pp_MeV->SetBinContent(iBinX, iBinY, hSigma_pp->GetBinContent(iBinX, iBinY));
    }
  }



  TFile** OliFile_pp = new TFile* [NumMultBins];
  TH1F** OliHisto_pp = new TH1F* [NumMultBins];
  TH1F** OliHisto_pp_SE = new TH1F* [NumMultBins];
  TH1F** OliHisto_pp_MeV = new TH1F* [NumMultBins];
  for(unsigned uMult=0; uMult<NumMultBins; uMult++){
    OliFile_pp[uMult] = new TFile(OliFileName_pp[uMult], "read");
    printf("OliFile_pp[%u] = %s\n", uMult, OliFileName_pp[uMult].Data());

    OliHisto_pp[uMult] = (TH1F*)OliFile_pp[uMult]->Get(MultHistoName[uMult]);

    printf("OliHisto_pp[%u] = %p\n", uMult, OliHisto_pp[uMult]);
    OliHisto_pp_MeV[uMult] = new TH1F(TString::Format("OliHisto_pp_MeV_%u",uMult), TString::Format("OliHisto_pp_MeV_%u",uMult),
                                      OliHisto_pp[uMult]->GetNbinsX(),OliHisto_pp[uMult]->GetBinLowEdge(1)*UnitsTransform[uMult],
                                      OliHisto_pp[uMult]->GetXaxis()->GetBinUpEdge(OliHisto_pp[uMult]->GetNbinsX())*UnitsTransform[uMult]);
    for(unsigned uBin=0; uBin<unsigned(OliHisto_pp[uMult]->GetNbinsX()); uBin++){
      OliHisto_pp_MeV[uMult]->SetBinContent(uBin+1, OliHisto_pp[uMult]->GetBinContent(uBin+1));
      OliHisto_pp_MeV[uMult]->SetBinError(uBin+1, OliHisto_pp[uMult]->GetBinError(uBin+1));
    }
    OliHisto_pp_MeV[uMult]->SetLineColor(kBlue+2);
    OliHisto_pp_MeV[uMult]->SetLineWidth(6);
    OliHisto_pp_SE[uMult] = (TH1F*)OliFile_pp[uMult]->Get(SEHistoName[uMult].Data());
  }


  NumPairs = new double[NumMultBins];
  NumPairs[0] = OliHisto_pp_SE[0]->Integral(OliHisto_pp_SE[0]->FindBin(0),OliHisto_pp_SE[0]->FindBin(200));
  NumPairs[1] = OliHisto_pp_SE[1]->Integral(OliHisto_pp_SE[1]->FindBin(0),OliHisto_pp_SE[1]->FindBin(0.200));
  NumPairs[2] = OliHisto_pp_SE[2]->Integral(OliHisto_pp_SE[2]->FindBin(0),OliHisto_pp_SE[2]->FindBin(0.200));
  NumPairs[3] = OliHisto_pp_SE[3]->Integral(OliHisto_pp_SE[3]->FindBin(0),OliHisto_pp_SE[3]->FindBin(0.200));
  std::cout << NumPairs[0] << '\t' << NumPairs[1] << '\t' << NumPairs[2] << '\t' << NumPairs[3] << '\n';
  //    DLM_DtColor DlmCol;

  TFile* GraphFile = new TFile(TString::Format("%sGraphFileMULT.root",OutputDir.Data()), "recreate");

  TGraph* FitResult_pp = new TGraph [NumMultBins];

  TGraph GraphWeightedSum_pp;
  GraphWeightedSum_pp.Set(NumMomBins_pp);
  GraphWeightedSum_pp.SetName("GraphWeightedSum_pp");
  GraphWeightedSum_pp.SetMarkerColor(fColors[NumMultBins]);
  GraphWeightedSum_pp.SetMarkerStyle(20);
  GraphWeightedSum_pp.SetMarkerSize(0);
  GraphWeightedSum_pp.SetLineColor(fColors[NumMultBins]);
  GraphWeightedSum_pp.SetLineWidth(6);
  GraphWeightedSum_pp.SetFillColor(kWhite);

  CATShisto<double> hWeightedSum_pp(NumMomBins_pp,kMin_pp,kMax_pp);
  TGraph GraphRatio_pp;
  GraphRatio_pp.Set(NumMomBins_pp);
  GraphRatio_pp.SetName("GraphRatio_pp");
  GraphRatio_pp.SetMarkerColor(kBlack);
  GraphRatio_pp.SetMarkerStyle(20);
  GraphRatio_pp.SetMarkerSize(0);
  GraphRatio_pp.SetLineColor(kBlack);
  GraphRatio_pp.SetLineWidth(6);
  GraphRatio_pp.SetFillColor(kWhite);

  TGraph GraphWeightedSumTheory_pp;
  GraphWeightedSumTheory_pp.Set(NumMomBins_pp);
  GraphWeightedSumTheory_pp.SetName("GraphWeightedSumTheory_pp");
  GraphWeightedSumTheory_pp.SetMarkerColor(NumMultBins+3);
  GraphWeightedSumTheory_pp.SetMarkerStyle(20);
  GraphWeightedSumTheory_pp.SetMarkerSize(0);
  GraphWeightedSumTheory_pp.SetLineColor(NumMultBins+3);
  GraphWeightedSumTheory_pp.SetLineWidth(8);
  GraphWeightedSumTheory_pp.SetFillColor(kWhite);

  CATShisto<double> hWeightedSumTheory_pp(NumMomBins_pp,kMin_pp,kMax_pp);
  TGraph GraphCkTheory_pp;
  GraphCkTheory_pp.Set(NumMomBins_pp);
  GraphCkTheory_pp.SetName("GraphCkTheory_pp");
  GraphCkTheory_pp.SetMarkerColor(fColors[NumMultBins+2]);
  GraphCkTheory_pp.SetMarkerStyle(20);
  GraphCkTheory_pp.SetMarkerSize(0);
  GraphCkTheory_pp.SetLineColor(fColors[NumMultBins+2]);
  GraphCkTheory_pp.SetLineWidth(6);
  GraphCkTheory_pp.SetFillColor(kWhite);

  //CATShisto<double> hCkTheory(NumMomBins_pp,kMin_pp,kMax_pp);
  TGraph GraphRatioTheory_pp;
  GraphRatioTheory_pp.Set(NumMomBins_pp);
  GraphRatioTheory_pp.SetName("GraphRatioTheory_pp");
  GraphRatioTheory_pp.SetMarkerColor(kBlack);
  GraphRatioTheory_pp.SetMarkerStyle(20);
  GraphRatioTheory_pp.SetMarkerSize(0);
  GraphRatioTheory_pp.SetLineColor(kBlack);
  GraphRatioTheory_pp.SetLineWidth(6);
  GraphRatioTheory_pp.SetFillColor(kWhite);
  //CATShisto<double> hRatioTheory(NumMomBins_pp,kMin_pp,kMax_pp);




  TGraph GraphWeightedSumTheory_pL;
  GraphWeightedSumTheory_pL.Set(NumMomBins_pL);
  GraphWeightedSumTheory_pL.SetName("GraphWeightedSumTheory_pL");
  GraphWeightedSumTheory_pL.SetMarkerColor(NumMultBins+3);
  GraphWeightedSumTheory_pL.SetMarkerStyle(20);
  GraphWeightedSumTheory_pL.SetMarkerSize(0);
  GraphWeightedSumTheory_pL.SetLineColor(NumMultBins+3);
  GraphWeightedSumTheory_pL.SetLineWidth(8);
  GraphWeightedSumTheory_pL.SetFillColor(kWhite);

  CATShisto<double> hWeightedSumTheory_pL(NumMomBins_pL,kMin_pL,kMax_pL);
  TGraph GraphCkTheory_pL;
  GraphCkTheory_pL.Set(NumMomBins_pL);
  GraphCkTheory_pL.SetName("GraphCkTheory_pL");
  GraphCkTheory_pL.SetMarkerColor(fColors[NumMultBins+2]);
  GraphCkTheory_pL.SetMarkerStyle(20);
  GraphCkTheory_pL.SetMarkerSize(0);
  GraphCkTheory_pL.SetLineColor(fColors[NumMultBins+2]);
  GraphCkTheory_pL.SetLineWidth(6);
  GraphCkTheory_pL.SetFillColor(kWhite);

  //CATShisto<double> hCkTheory(NumMomBins_pL,kMin_pL,kMax_pL);
  TGraph GraphRatioTheory_pL;
  GraphRatioTheory_pL.Set(NumMomBins_pL);
  GraphRatioTheory_pL.SetName("GraphRatioTheory_pL");
  GraphRatioTheory_pL.SetMarkerColor(kBlack);
  GraphRatioTheory_pL.SetMarkerStyle(20);
  GraphRatioTheory_pL.SetMarkerSize(0);
  GraphRatioTheory_pL.SetLineColor(kBlack);
  GraphRatioTheory_pL.SetLineWidth(6);
  GraphRatioTheory_pL.SetFillColor(kWhite);
  //CATShisto<double> hRatioTheory(NumMomBins_pL,kMin_pL,kMax_pL);


  const unsigned PixelsX = 1920;
  const unsigned PixelsY = 1080;
  const double MarginL = 0.12;
  const double MarginR = 0.05;
  const double MarginB = 0.18;
  const double MarginT = 0.05;

  TCanvas* canWeightedSum_pp = new TCanvas("canWeightedSum_pp", "canWeightedSum_pp", 1);
  canWeightedSum_pp->cd(0); canWeightedSum_pp->SetCanvasSize(PixelsX, PixelsY); canWeightedSum_pp->SetMargin(MarginL,MarginR,MarginB,MarginT);//lrbt
  TLegend* legWeightedSum_pp = new TLegend(0.55,0.75,0.95,0.95);//lbrt
  legWeightedSum_pp->SetName(TString::Format("legWeightedSum_pp"));
  legWeightedSum_pp->SetTextSize(0.05);

  TCanvas* canRatio_pp = new TCanvas("canRatio_pp", "canRatio_pp", 1);
  canRatio_pp->cd(0); canRatio_pp->SetCanvasSize(PixelsX, PixelsY); canRatio_pp->SetMargin(MarginL,MarginR,MarginB,MarginT);//lrbt

  TCanvas** canMultBins_pp = new TCanvas* [NumMultBins];
  TLegend** legMultBins_pp = new TLegend* [NumMultBins];
  for(unsigned uMult=0; uMult<NumMultBins; uMult++){
    canMultBins_pp[uMult] = new TCanvas(TString::Format("canMultBins_pp_%u",uMult), TString::Format("canMultBins_pp_%u",uMult), 1);
    canMultBins_pp[uMult]->cd(0); canMultBins_pp[uMult]->SetCanvasSize(PixelsX, PixelsY); canMultBins_pp[uMult]->SetMargin(MarginL,MarginR,MarginB,MarginT);//lrbt

    legMultBins_pp[uMult] = new TLegend(0.5,0.75,0.95,0.95);//lbrt
    legMultBins_pp[uMult]->SetName(TString::Format("legMultBins_pp_%u",uMult));
    legMultBins_pp[uMult]->SetTextSize(0.045);
  }

  TCanvas* canWeightedSumTheory_pp = new TCanvas("canWeightedSumTheory_pp", "canWeightedSumTheory_pp", 1);
  canWeightedSumTheory_pp->cd(0); canWeightedSumTheory_pp->SetCanvasSize(PixelsX, PixelsY); canWeightedSumTheory_pp->SetMargin(MarginL,MarginR,MarginB,MarginT);//lrbt
  TLegend* legWeightedSumTheory_pp = new TLegend(0.55,0.75,0.95,0.95);//lbrt
  legWeightedSumTheory_pp->SetName(TString::Format("legWeightedSumTheory_pp"));
  legWeightedSumTheory_pp->SetTextSize(0.05);

  TCanvas* canRatioTheory_pp = new TCanvas("canRatioTheory_pp", "canRatioTheory_pp", 1);
  canRatioTheory_pp->cd(0); canRatioTheory_pp->SetCanvasSize(PixelsX, PixelsY); canRatioTheory_pp->SetMargin(MarginL,MarginR,MarginB,MarginT);//lrbt

  TH1F* hAxisCk_pp = new TH1F("hAxisCk_pp", "hAxisCk_pp", NumMomBins_pp, kMin_pp, kMax_pp);
  hAxisCk_pp->SetStats(false);
  hAxisCk_pp->SetTitle("");
  hAxisCk_pp->GetXaxis()->SetLabelSize(0.065);
  hAxisCk_pp->GetXaxis()->SetTitle("k (MeV)");
  hAxisCk_pp->GetXaxis()->CenterTitle();
  hAxisCk_pp->GetXaxis()->SetTitleOffset(1.15);
  hAxisCk_pp->GetXaxis()->SetLabelOffset(0.02);
  hAxisCk_pp->GetXaxis()->SetTitleSize(0.075);
  hAxisCk_pp->GetYaxis()->SetLabelSize(0.065);
  hAxisCk_pp->GetYaxis()->SetTitle("C(k)");
  hAxisCk_pp->GetYaxis()->CenterTitle();
  hAxisCk_pp->GetYaxis()->SetTitleOffset(0.75);
  hAxisCk_pp->GetYaxis()->SetTitleSize(0.075);

  if(DataSet=="ALICE_pp_13TeV"){hAxisCk_pp->GetYaxis()->SetRangeUser(0.5, 4.5);}
  else if(DataSet=="ALICE_pPb_5TeV"){hAxisCk_pp->GetYaxis()->SetRangeUser(0.5, 3.5);}

  //hAxisCk_pp->GetYaxis()->SetLimits(0.5, 4.5);

  TH1F* hAxisRatio_pp = new TH1F("hAxisRatio_pp", "hAxisRatio_pp", NumMomBins_pp,kMin_pp, kMax_pp);
  hAxisRatio_pp->SetStats(false);
  hAxisRatio_pp->SetTitle("");
  hAxisRatio_pp->GetXaxis()->SetLabelSize(0.065);
  hAxisRatio_pp->GetXaxis()->SetTitle("k (MeV)");
  hAxisRatio_pp->GetXaxis()->CenterTitle();
  hAxisRatio_pp->GetXaxis()->SetTitleOffset(1.15);
  hAxisRatio_pp->GetXaxis()->SetLabelOffset(0.02);
  hAxisRatio_pp->GetXaxis()->SetTitleSize(0.075);
  hAxisRatio_pp->GetYaxis()->SetLabelSize(0.065);
  hAxisRatio_pp->GetYaxis()->SetTitle("Ratio");
  hAxisRatio_pp->GetYaxis()->CenterTitle();
  hAxisRatio_pp->GetYaxis()->SetTitleOffset(0.75);
  hAxisRatio_pp->GetYaxis()->SetTitleSize(0.075);

  hAxisRatio_pp->GetYaxis()->SetRangeUser(0.9, 1.1);
  //hAxisRatio_pp->GetYaxis()->SetLimits(0.9, 1.1);


  TCanvas* canWeightedSumTheory_pL = new TCanvas("canWeightedSumTheory_pL", "canWeightedSumTheory_pL", 1);
  canWeightedSumTheory_pL->cd(0); canWeightedSumTheory_pL->SetCanvasSize(PixelsX, PixelsY); canWeightedSumTheory_pL->SetMargin(MarginL,MarginR,MarginB,MarginT);//lrbt
  TLegend* legWeightedSumTheory_pL = new TLegend(0.55,0.75,0.95,0.95);//lbrt
  legWeightedSumTheory_pL->SetName(TString::Format("legWeightedSumTheory_pL"));
  legWeightedSumTheory_pL->SetTextSize(0.05);

  TCanvas* canRatioTheory_pL = new TCanvas("canRatioTheory_pL", "canRatioTheory_pL", 1);
  canRatioTheory_pL->cd(0); canRatioTheory_pL->SetCanvasSize(PixelsX, PixelsY); canRatioTheory_pL->SetMargin(MarginL,MarginR,MarginB,MarginT);//lrbt

  TH1F* hAxisCk_pL = new TH1F("hAxisCk_pL", "hAxisCk_pL", NumMomBins_pL, kMin_pL, kMax_pL);
  hAxisCk_pL->SetStats(false);
  hAxisCk_pL->SetTitle("");
  hAxisCk_pL->GetXaxis()->SetLabelSize(0.065);
  hAxisCk_pL->GetXaxis()->SetTitle("k (MeV)");
  hAxisCk_pL->GetXaxis()->CenterTitle();
  hAxisCk_pL->GetXaxis()->SetTitleOffset(1.15);
  hAxisCk_pL->GetXaxis()->SetLabelOffset(0.02);
  hAxisCk_pL->GetXaxis()->SetTitleSize(0.075);
  hAxisCk_pL->GetYaxis()->SetLabelSize(0.065);
  hAxisCk_pL->GetYaxis()->SetTitle("C(k)");
  hAxisCk_pL->GetYaxis()->CenterTitle();
  hAxisCk_pL->GetYaxis()->SetTitleOffset(0.75);
  hAxisCk_pL->GetYaxis()->SetTitleSize(0.075);

  if(DataSet=="ALICE_pp_13TeV"){hAxisCk_pL->GetYaxis()->SetRangeUser(0.5, 4.5);}
  else if(DataSet=="ALICE_pPb_5TeV"){hAxisCk_pL->GetYaxis()->SetRangeUser(0.5, 3.5);}
  //hAxisCk_pL->GetYaxis()->SetRangeUser(0.5, 4.5);
  //hAxisCk_pL->GetYaxis()->SetLimits(0.5, 4.5);

  TH1F* hAxisRatio_pL = new TH1F("hAxisRatio_pL", "hAxisRatio_pL", NumMomBins_pL,kMin_pL, kMax_pL);
  hAxisRatio_pL->SetStats(false);
  hAxisRatio_pL->SetTitle("");
  hAxisRatio_pL->GetXaxis()->SetLabelSize(0.065);
  hAxisRatio_pL->GetXaxis()->SetTitle("k (MeV)");
  hAxisRatio_pL->GetXaxis()->CenterTitle();
  hAxisRatio_pL->GetXaxis()->SetTitleOffset(1.15);
  hAxisRatio_pL->GetXaxis()->SetLabelOffset(0.02);
  hAxisRatio_pL->GetXaxis()->SetTitleSize(0.075);
  hAxisRatio_pL->GetYaxis()->SetLabelSize(0.065);
  hAxisRatio_pL->GetYaxis()->SetTitle("Ratio");
  hAxisRatio_pL->GetYaxis()->CenterTitle();
  hAxisRatio_pL->GetYaxis()->SetTitleOffset(0.75);
  hAxisRatio_pL->GetYaxis()->SetTitleSize(0.075);

  hAxisRatio_pL->GetYaxis()->SetRangeUser(0.9, 1.1);
  hAxisRatio_pL->GetYaxis()->SetLimits(0.9, 1.1);


  int vSource;
  int vFit;

  TRandom3 rangen(1);

  vSource=0;//0 Gauss, 1 Cauchy, 2 DoubleGauss, 3 EPOS
  vFit=1;//1 = with BL


  double FemtoRegion_pp[2];
  FemtoRegion_pp[0]=0; FemtoRegion_pp[1]=200;

  double BlRegion[2];
  BlRegion[0]=300; BlRegion[1]=500;

  const double GaussSourceSize = 1.2;
  const double Mass_p=938.272; const double Mass_L=1115.683;

  //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
  double PotPars1S0[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
  double PotPars3P0[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
  double PotPars3P1[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
  double PotPars3P2[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};

  double pLamPotPars1S0[10]={0,0,pL_UsmaniOli,0,0,0,0,0,0,0};
  double pLamPotPars3S1[10]={0,0,pL_UsmaniOli,0,0,0,0,1,0,1};

  const double Weight1S0 = 3./12.;
  const double Weight3P0 = 1./12.;
  const double Weight3P1 = 3./12.;
  const double Weight3P2 = 5./12.;

  CATS AB_pp[NumMultBins];
  double Pars_pp[NumMultBins][6];


  CATS AB_pL[NumMultBins];
  double Pars_pL[NumMultBins][6];

  for(unsigned uMult=0; uMult<NumMultBins; uMult++){
    Pars_pp[uMult][0]=0; Pars_pp[uMult][1]=0; Pars_pp[uMult][2]=0;
    Pars_pp[uMult][3]=GaussSourceSize*1.2; Pars_pp[uMult][4]=GaussSourceSize/1.2; Pars_pp[uMult][5]=0.5;

    Pars_pL[uMult][0]=0; Pars_pL[uMult][1]=0; Pars_pL[uMult][2]=0;
    Pars_pL[uMult][3]=GaussSourceSize*1.2; Pars_pp[uMult][4]=GaussSourceSize/1.2; Pars_pp[uMult][5]=0.5;

    if(vSource==0){
      AB_pp[uMult].SetAnaSource(GaussSource, Pars_pp[uMult]);
      AB_pp[uMult].SetUseAnalyticSource(true);
      AB_pp[uMult].SetThetaDependentSource(false);
    }
    else if(vSource==1){
      AB_pp[uMult].SetAnaSource(CauchySource, Pars_pp[uMult]);
      AB_pp[uMult].SetUseAnalyticSource(true);
      AB_pp[uMult].SetThetaDependentSource(false);
    }
    else if(vSource==2){
      AB_pp[uMult].SetAnaSource(DoubleGaussSource, Pars_pp[uMult]);
      AB_pp[uMult].SetUseAnalyticSource(true);
      AB_pp[uMult].SetThetaDependentSource(false);
    }
    else if(vSource==3){
      AB_pp[uMult].SetUseAnalyticSource(false);
      AB_pp[uMult].SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19");
      if(EPOS_MAX_PAIRS) AB_pp[uMult].SetMaxPairsToRead(EPOS_MAX_PAIRS);
      AB_pp[uMult].SetMaxPairsPerBin(64000);
      AB_pp[uMult].SetMixingDepth(EPOS_MIX_DEP);
      AB_pp[uMult].SetThetaDependentSource(EPOS_THETA_DEP);
      AB_pp[uMult].SetTransportRenorm(1);
      if(EPOS_THETA_DEP){
        AB_pp[uMult].SetGridEpsilon(1./8192);
        AB_pp[uMult].SetGridMaxDepth(8);
      }
    }

    AB_pp[uMult].SetExcludeFailedBins(false);
    AB_pp[uMult].SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);

    AB_pp[uMult].SetNumChannels(4);
    AB_pp[uMult].SetNumPW(0,2);
    AB_pp[uMult].SetNumPW(1,2);
    AB_pp[uMult].SetNumPW(2,2);
    AB_pp[uMult].SetNumPW(3,2);
    AB_pp[uMult].SetSpin(0,0);
    AB_pp[uMult].SetSpin(1,1);
    AB_pp[uMult].SetSpin(2,1);
    AB_pp[uMult].SetSpin(3,1);
    AB_pp[uMult].SetChannelWeight(0, Weight1S0);
    AB_pp[uMult].SetChannelWeight(1, Weight3P0);
    AB_pp[uMult].SetChannelWeight(2, Weight3P1);
    AB_pp[uMult].SetChannelWeight(3, Weight3P2);

    AB_pp[uMult].SetQ1Q2(1);
    AB_pp[uMult].SetPdgId(2212, 2212);
    AB_pp[uMult].SetRedMass( 0.5*Mass_p );

    AB_pp[uMult].SetShortRangePotential(0,0,fDlmPot,PotPars1S0);
    AB_pp[uMult].SetShortRangePotential(1,1,fDlmPot,PotPars3P0);
    AB_pp[uMult].SetShortRangePotential(2,1,fDlmPot,PotPars3P1);
    AB_pp[uMult].SetShortRangePotential(3,1,fDlmPot,PotPars3P2);

    AB_pp[uMult].KillTheCat();

    if(vSource==0){
      AB_pL[uMult].SetAnaSource(GaussSource, Pars_pL[uMult]);
      AB_pL[uMult].SetUseAnalyticSource(true);
      AB_pL[uMult].SetThetaDependentSource(false);
    }
    else if(vSource==1){
      AB_pL[uMult].SetAnaSource(CauchySource, Pars_pL[uMult]);
      AB_pL[uMult].SetUseAnalyticSource(true);
      AB_pL[uMult].SetThetaDependentSource(false);
    }
    else if(vSource==2){
      AB_pL[uMult].SetAnaSource(DoubleGaussSource, Pars_pL[uMult]);
      AB_pL[uMult].SetUseAnalyticSource(true);
      AB_pL[uMult].SetThetaDependentSource(false);
    }
    else if(vSource==3){
      //            AB_pL[uMult].SetUseAnalyticSource(false);
      //            AB_pL[uMult].SetInputFileName(COMPUTER_ID==cNX1?"/home/gu47det/EPOS/EPOS_FILE_READER/OutputFiles/f19/Scratch9_OSCAR1997_AllLrzLambdaFiles.f19":
      //                                   "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzLambdaFiles.f19");
      //            if(EPOS_MAX_PAIRS) AB_pL[uMult].SetMaxPairsToRead(EPOS_MAX_PAIRS);
      //            AB_pL[uMult].SetMaxPairsPerBin(64000);
      //            AB_pL[uMult].SetMixingDepth(EPOS_MIX_DEP);
      //            AB_pL[uMult].SetThetaDependentSource(EPOS_THETA_DEP);
      //            AB_pL[uMult].SetTransportRenorm(1);
      //            if(EPOS_THETA_DEP){
      //                AB_pL[uMult].SetGridEpsilon(1./8192);
      //                AB_pL[uMult].SetGridMaxDepth(8);
      //            }
    }
    AB_pL[uMult].SetMomBins(NumMomBins_pL,kMin_pL,kMax_pL);

    AB_pL[uMult].SetNumChannels(2);
    AB_pL[uMult].SetNumPW(0,1);
    AB_pL[uMult].SetNumPW(1,1);
    AB_pL[uMult].SetSpin(0,0);
    AB_pL[uMult].SetSpin(1,1);
    AB_pL[uMult].SetChannelWeight(0, 0.25);
    AB_pL[uMult].SetChannelWeight(1, 0.75);

    AB_pL[uMult].SetQ1Q2(0);
    AB_pL[uMult].SetPdgId(2212, 3122);
    AB_pL[uMult].SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

    AB_pL[uMult].SetShortRangePotential(0,0,fDlmPot,pLamPotPars1S0);
    AB_pL[uMult].SetShortRangePotential(1,0,fDlmPot,pLamPotPars3S1);

    AB_pL[uMult].KillTheCat();
  }

  //!IMPORTANT SETTINGS
  const unsigned NumSourcePars =  vSource==0?1:
      vSource==1?1:
          vSource==2?3:
              1;

  DLM_Ck** Ck_pp = new DLM_Ck* [NumMultBins];
  DLM_Ck** Ck_pL = new DLM_Ck* [NumMultBins];

  DLM_CkDecomposition** CkDec_pp = new DLM_CkDecomposition*[NumMultBins];
  DLM_CkDecomposition** CkDec_pL = new DLM_CkDecomposition*[NumMultBins];

  double lam_pp;
  double lam_pp_pL;
  double lam_pp_fake;

  double lam_pL;
  double lam_pL_Feed;
  double lam_pL_fake;

  if(DataSet=="ALICE_pp_13TeV"){
    lam_pp = 0.752;
    lam_pp_pL = 0.151;
    lam_pp_fake = 1.-lam_pp-lam_pp_pL;

    lam_pL = 0.519;
    lam_pL_fake = 0.043;
    lam_pL_Feed = 1.-lam_pL-lam_pL_fake;
  }
  else if(DataSet=="ALICE_pPb_5TeV"){
    lam_pp = 0.720;
    lam_pp_pL = 0.161;
    lam_pp_fake = 1.-lam_pp-lam_pp_pL;

    lam_pL = 0.589;
    lam_pL_fake = 0.083;
    lam_pL_Feed = 1.-lam_pL-lam_pL_fake;
  }


  DLM_Fitter1** fTest1 = new DLM_Fitter1*[NumMultBins];

  double CkVal;
  double Momentum;
  double ExpRad=0;

  for(unsigned uMult=0; uMult<NumMultBins; uMult++){
    Ck_pp[uMult] = new DLM_Ck(NumSourcePars,0,AB_pp[uMult]);
    Ck_pL[uMult] = new DLM_Ck(NumSourcePars,0,AB_pL[uMult]);

    Ck_pp[uMult]->Update();
    Ck_pL[uMult]->Update();

    CkDec_pp[uMult] = new DLM_CkDecomposition("pp",3,*Ck_pp[uMult],hSigma_pp_MeV);
    CkDec_pL[uMult] = new DLM_CkDecomposition("pLambda",2,*Ck_pL[uMult],NULL);

    CkDec_pp[uMult]->AddContribution(0,lam_pp_pL,DLM_CkDecomposition::cFeedDown,CkDec_pL[uMult],hRes_pp_pL_MeV);
    CkDec_pp[uMult]->AddContribution(1,1.-lam_pp-lam_pp_pL-lam_pp_fake,DLM_CkDecomposition::cFeedDown);
    CkDec_pp[uMult]->AddContribution(2,lam_pp_fake,DLM_CkDecomposition::cFake);//0.02

    CkDec_pL[uMult]->AddContribution(0,lam_pL_Feed,DLM_CkDecomposition::cFeedDown);
    CkDec_pL[uMult]->AddContribution(1,lam_pL_fake,DLM_CkDecomposition::cFake);//0.03

    fTest1[uMult] = new DLM_Fitter1(1);
    fTest1[uMult]->SetOutputDir(OutputDir.Data());

    if(vFit==0){
      fTest1[uMult]->SetSystem(0,*OliHisto_pp_MeV[uMult],1,*CkDec_pp[uMult],
                               FemtoRegion_pp[0],FemtoRegion_pp[1],
                               FemtoRegion_pp[1],FemtoRegion_pp[1]);
    }
    else{
      fTest1[uMult]->SetSystem(0,*OliHisto_pp_MeV[uMult],1,*CkDec_pp[uMult],
                               FemtoRegion_pp[0],FemtoRegion_pp[1],
                               BlRegion[0],BlRegion[1]);
    }

    fTest1[uMult]->SetSeparateBL(0,vFit!=0);

    fTest1[uMult]->SetParameter("pp",DLM_Fitter1::p_a,1.0,0.7,1.3);
    if(vFit<1) fTest1[uMult]->FixParameter("pp",DLM_Fitter1::p_b,0);
    else fTest1[uMult]->SetParameter("pp",DLM_Fitter1::p_b,1e-4,0,2e-3);
    if(vFit<2) fTest1[uMult]->FixParameter("pp",DLM_Fitter1::p_c,0);
    else fTest1[uMult]->SetParameter("pp",DLM_Fitter1::p_c,0,-2e-4,2e-4);
    if(vSource==0){
      fTest1[uMult]->SetParameter("pp",DLM_Fitter1::p_sor0,GaussSourceSize,0.,2.0);
      //fTest1->FixParameter("pp",DLM_Fitter1::p_sor0,1.2);
    }
    else if(vSource==1){
      fTest1[uMult]->SetParameter("pp",DLM_Fitter1::p_sor0,GaussSourceSize/1.4,0.5,1.6);
    }
    else if(vSource==2){
      fTest1[uMult]->SetParameter("pp",DLM_Fitter1::p_sor0,0.8,0.4,3.6);
      fTest1[uMult]->SetParameter("pp",DLM_Fitter1::p_sor1,2.4,0.4,3.6);
      fTest1[uMult]->SetParameter("pp",DLM_Fitter1::p_sor2,0.5,0,1);
    }
    else if(vSource==3){
      fTest1[uMult]->SetParameter("pp",DLM_Fitter1::p_sor0,1.5,1.35,1.65);
      //fTest1->FixParameter("pp",DLM_Fitter1::p_sor0,1.5);
    }

    fTest1[uMult]->SetParameter("pp",DLM_Fitter1::p_Cl,-0.9,-1.2,-0.8);
    fTest1[uMult]->SetSeparateBL(0,false);

    CkDec_pp[uMult]->Update();
    CkDec_pL[uMult]->Update();
    std::cout << "I am in Iteration Number: " << uMult <<std::endl;
    fTest1[uMult]->GoBabyGo();


    AB_pL[uMult].SetAnaSource(0, fTest1[uMult]->GetParameter("pp",DLM_Fitter1::p_sor0));
    AB_pL[uMult].KillTheCat();

    printf("For uMult=%u I got r = %.3f fm (w=%.3f)\n", uMult, fTest1[uMult]->GetParameter("pp",DLM_Fitter1::p_sor0),NumPairs[uMult]/NumPairs[0]);
    if(uMult!=0) ExpRad += fTest1[uMult]->GetParameter("pp",DLM_Fitter1::p_sor0)*NumPairs[uMult]/NumPairs[0];

    FitResult_pp[uMult].SetName(TString::Format("FitResult_pp_%u",uMult));
    fTest1[uMult]->GetFitGraph(0, FitResult_pp[uMult]);
    FitResult_pp[uMult].SetMarkerColor(fColors[uMult]);
    FitResult_pp[uMult].SetMarkerStyle(20);
    FitResult_pp[uMult].SetMarkerSize(0);
    FitResult_pp[uMult].SetLineColor(fColors[uMult]);
    FitResult_pp[uMult].SetLineWidth(6);
    FitResult_pp[uMult].SetFillColor(kWhite);

    GraphFile->cd();
    FitResult_pp[uMult].Write();

    TString FitString = TString::Format("Fit (r=%.2f#pm%.2f fm",fTest1[uMult]->GetParameter("pp",DLM_Fitter1::p_sor0),fTest1[uMult]->GetParError("pp",DLM_Fitter1::p_sor0));
    if(uMult) FitString += TString::Format(", w=%.0f%%)",100.*NumPairs[uMult]/NumPairs[0]);
    else FitString += ")";
    canMultBins_pp[uMult]->cd(0);
    TString legName =uMult?TString::Format("Data for MB %u", uMult):"Total data";
    legMultBins_pp[uMult]->AddEntry(OliHisto_pp_MeV[uMult], legName.Data());
    legMultBins_pp[uMult]->AddEntry(&FitResult_pp[uMult], FitString);
    hAxisCk_pp->Draw("axis");
    OliHisto_pp_MeV[uMult]->Draw("same");
    FitResult_pp[uMult].Draw("C,same");
    legMultBins_pp[uMult]->Draw("same");


    for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
      if(uMult!=0){
        FitResult_pp[uMult].GetPoint(uBin,Momentum,CkVal);
        hWeightedSum_pp.Add(uBin,CkVal*NumPairs[uMult]/NumPairs[0]);
        hWeightedSumTheory_pp.Add(uBin,CkDec_pp[uMult]->EvalCk(Momentum)*NumPairs[uMult]/NumPairs[0]);
        //GraphWeightedSumTheory.SetPoint(uBin,Momentum,CkDec_pp[uMult]->EvalCk(Momentum)*NumPairs[uMult]/NumPairs[0]);
      }
      else{
        Momentum = AB_pp[uMult].GetMomentum(uBin);
        GraphCkTheory_pp.SetPoint(uBin,Momentum,CkDec_pp[uMult]->EvalCk(Momentum));

      }
    }

    for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){
      if(uMult!=0){
        Momentum = Ck_pL[uMult]->GetBinCenter(uBin);
        printf("pL: %.2f MeV: %.3f\n", Momentum, AB_pL[uMult].GetCorrFun(uBin));
        hWeightedSumTheory_pL.Add(uBin,AB_pL[uMult].GetCorrFun(uBin)*NumPairs[uMult]/NumPairs[0]);
      }
      else{
        Momentum = Ck_pL[uMult]->GetBinCenter(uBin);
        GraphCkTheory_pL.SetPoint(uBin,Momentum,AB_pL[uMult].GetCorrFun(uBin));
      }
    }


  }//uMult

  printf("Expected global radius (weighted)  = %.3f\n", ExpRad);

  for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
    FitResult_pp[0].GetPoint(uBin,Momentum,CkVal);
    GraphWeightedSum_pp.SetPoint(uBin,hWeightedSum_pp.GetBinCenter(uBin),hWeightedSum_pp.GetBinContent(uBin));
    GraphRatio_pp.SetPoint(uBin,hWeightedSum_pp.GetBinCenter(uBin),hWeightedSum_pp.GetBinContent(uBin)/CkVal);

    GraphCkTheory_pp.GetPoint(uBin,Momentum,CkVal);
    GraphWeightedSumTheory_pp.SetPoint(uBin,hWeightedSumTheory_pp.GetBinCenter(uBin),hWeightedSumTheory_pp.GetBinContent(uBin));
    GraphRatioTheory_pp.SetPoint(uBin,hWeightedSumTheory_pp.GetBinCenter(uBin),hWeightedSumTheory_pp.GetBinContent(uBin)/CkVal);
  }

  canRatio_pp->cd(0);
  hAxisRatio_pp->Draw("axis");
  GraphRatio_pp.Draw("C,same");
  canRatio_pp->SaveAs(TString::Format("%s/canRatio_pp.png",OutputDir.Data()));

  canRatioTheory_pp->cd(0);
  hAxisRatio_pp->Draw("axis");
  GraphRatioTheory_pp.Draw("C,same");
  canRatioTheory_pp->SaveAs(TString::Format("%s/canRatioTheory_pp.png",OutputDir.Data()));

  canWeightedSum_pp->cd(0);
  legWeightedSum_pp->AddEntry(OliHisto_pp_MeV[0],"Total data");
  legWeightedSum_pp->AddEntry(&GraphWeightedSum_pp,"Weighted C(k)");
  legWeightedSum_pp->AddEntry(&FitResult_pp[0],"C(k) (fit)");
  hAxisCk_pp->Draw("axis");
  OliHisto_pp_MeV[0]->Draw("same");
  GraphWeightedSum_pp.Draw("C,same");
  FitResult_pp[0].Draw("C,same");
  legWeightedSum_pp->Draw("same");
  canWeightedSum_pp->SaveAs(TString::Format("%s/canWeightedSum_pp.png",OutputDir.Data()));

  canWeightedSumTheory_pp->cd(0);
  legWeightedSumTheory_pp->AddEntry(&GraphWeightedSumTheory_pp,"Weighted C_{th}(k)");
  legWeightedSumTheory_pp->AddEntry(&GraphCkTheory_pp,"C_{th}(k)");
  hAxisCk_pp->Draw("axis");
  //OliHisto_pp_MeV[0]->Draw("same");
  GraphWeightedSumTheory_pp.Draw("C,same");
  GraphCkTheory_pp.Draw("C,same");
  legWeightedSumTheory_pp->Draw("same");
  canWeightedSumTheory_pp->SaveAs(TString::Format("%s/canWeightedSumTheory_pp.png",OutputDir.Data()));

  for(unsigned uMult=0; uMult<NumMultBins; uMult++){
    canMultBins_pp[uMult]->SaveAs(TString::Format("%s/canMultBins_pp_%u.png",OutputDir.Data(),uMult));
  }

  GraphFile->cd();
  GraphWeightedSum_pp.Write();
  GraphRatio_pp.Write();

  GraphCkTheory_pp.Write();
  GraphWeightedSumTheory_pp.Write();
  GraphRatioTheory_pp.Write();






  for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){
    GraphCkTheory_pL.GetPoint(uBin,Momentum,CkVal);
    GraphWeightedSumTheory_pL.SetPoint(uBin,hWeightedSumTheory_pL.GetBinCenter(uBin),hWeightedSumTheory_pL.GetBinContent(uBin));
    GraphRatioTheory_pL.SetPoint(uBin,hWeightedSumTheory_pL.GetBinCenter(uBin),hWeightedSumTheory_pL.GetBinContent(uBin)/CkVal);
  }

  canRatioTheory_pL->cd(0);
  hAxisRatio_pL->Draw("axis");
  GraphRatioTheory_pL.Draw("C,same");
  canRatioTheory_pL->SaveAs(TString::Format("%s/canRatioTheory_pL.png",OutputDir.Data()));


  canWeightedSumTheory_pL->cd(0);
  legWeightedSumTheory_pL->AddEntry(&GraphWeightedSumTheory_pL,"Weighted C_{th}(k)");
  legWeightedSumTheory_pL->AddEntry(&GraphCkTheory_pL,"C_{th}(k)");
  hAxisCk_pL->Draw("axis");
  //OliHisto_pL_MeV[0]->Draw("same");
  GraphWeightedSumTheory_pL.Draw("C,same");
  GraphCkTheory_pL.Draw("C,same");
  legWeightedSumTheory_pL->Draw("same");
  canWeightedSumTheory_pL->SaveAs(TString::Format("%s/canWeightedSumTheory_pL.png",OutputDir.Data()));


  GraphCkTheory_pL.Write();
  GraphWeightedSumTheory_pL.Write();
  GraphRatioTheory_pL.Write();




  delete hRes_pp_pL_MeV;
  delete hSigma_pp_MeV;

  delete hAxisCk_pp;
  delete hAxisRatio_pp;

  for(unsigned uMult=0; uMult<NumMultBins; uMult++){
    delete OliHisto_pp_MeV[uMult];
    delete Ck_pp[uMult];
    delete Ck_pL[uMult];
    delete CkDec_pp[uMult];
    delete CkDec_pL[uMult];
    delete fTest1[uMult];
  }
  delete [] OliHisto_pp_MeV;
  delete [] Ck_pp;
  delete [] Ck_pL;
  delete [] CkDec_pp;
  delete [] CkDec_pL;
  delete [] fTest1;

  //delete [] FitResult_pp;

  delete [] NumPairs;

  //delete [] FitResult_pL;

  delete hAxisCk_pL;
  delete hAxisRatio_pL;



  delete canRatio_pp;
  delete canRatioTheory_pp;
  delete canWeightedSum_pp;
  delete canWeightedSumTheory_pp;
  for(unsigned uMult=0; uMult<NumMultBins; uMult++) delete canMultBins_pp[uMult];
  delete [] canMultBins_pp;

  delete canRatioTheory_pL;
  delete canWeightedSumTheory_pL;
  return;
  delete GraphFile;

  delete FileSigma;
  delete FileRes;

  for(unsigned uMult=0; uMult<NumMultBins; uMult++){
    delete OliFile_pp[uMult];

  }
  delete [] OliFile_pp;

}

void EXECUTE_PROTON_ONLY(int argc, char *argv[]) {
  if (atoi(argv[1]) == 1) {
    std::cout << "We Are Running Multiplicity Binning \n";
    PROTON_ONLY("~/Analysis/CATS_Input/CentMultPPCF","MultBinnedProtonProtonCF.root",
                "~/Analysis/CATS_OUTPUT/CentMultPPCF/MultOutput");
  } else if (atoi(argv[1]) == 2) {
    std::cout << "We Are Running Centrality Binning \n";
    PROTON_ONLY("~/Analysis/CATS_Input/CentMultPPCF","CentBinnedProtonProtonCF.root",
                "~/Analysis/CATS_OUTPUT/CentMultPPCF/CentOutput");
  } else {
    std::cout << "Oh buoy this is not a viable option to choose input files. Rethink your life! \n ";
  }
    return;
  }
