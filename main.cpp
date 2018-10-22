#include <iostream>

#include "CATStools.h"
#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"

#include "TGraph.h"
#include "TFile.h"
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
#include "ForBernie.h"
#include "ProtonOnly.h"
#include "ProtonkT.h"

using namespace std;
void RUN2_SYSTEMATICS_MEDIAN(TString InputFolder, TString OutDirName);


int main(int argc, char *argv[])
{
  //FitkTBins(argv[1], argv[2], atoi(argv[3]));
  //EXECUTE_PROTON_ONLY(argc, argv);
  CALL_BERNIE_AND_VALE(argc,argv);
	//RUN2_SYSTEMATICS_MEDIAN(argv[1],argv[1]);
    return 0;
}



//computes the median and so on.
void RUN2_SYSTEMATICS_MEDIAN(TString InputFolder, TString OutDirName){
  const TString& SystematicsType="CutVarAdd";
  const TString& Scenario="Global_Radius";
  //const TString Scenario = "Global_Radius";
  //const TString Scenario = "pp_Radius";
  //const TString Scenario = "pL_Radius";

  //perform the systematics by using the different cut variations as further permutations
  //TString SystematicsType = "CutVarIterate";
  //perform the systematics by adding the systematics errors of the bin quadratically to the data points
  //TString SystematicsType = "CutVarAdd";

  //const double DefaultRad = SystematicsType=="CutVarIterate"?1.153:1.150;
  double DefaultRad; //= SystematicsType=="CutVarIterate"?1.196:1.185;//for the Gauss+LO run2
  const bool ShowEPOS = false;

  //const TString InputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputFilesCATS2.0/Systematics_February2018/Iter16384_FINAL2/";
//  const TString InputFolder = "";
  //const TString InputRootName = SystematicsType+"_Iter6144.root";
  //const unsigned NumIterResults = ShowEPOS?359:511;
  const unsigned NumIterResults = 2048;
  //const unsigned NumIterResults = 1;
  const TString InputRootName = SystematicsType+TString::Format("_Iter%u.root",NumIterResults);
  const TString InputGraphBase = "GraphFile_"+SystematicsType+"_Iter2048_uIter";

  //const TString InputRootName = "CutVarIterate_Iter6144.root";
  //const TString InputGraphBase = "GraphFile_CutVarIterate_Iter16384_uIter";


  const unsigned MinIter = 0;
  const unsigned MaxIter = (NumIterResults)-1+1;

  //const unsigned MaxIter = 5120-1;
  const unsigned NumIter = MaxIter-MinIter+1;
//  const TString OutDirName = ShowEPOS?"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputRun2/RUN2_SYSTEMATICS_Express2/pPb5TeV/SystFarm020718/":
//      "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputRun2/RUN2_SYSTEMATICS_Express2/pPb5TeV/SystFarm020718/";
  //const TString DataDescr = "NEW";
  //const TString DataDescr = "QuadrSyst";
  //const TString DataDescr = "Bernie";
  const TString OutFileName = OutDirName+"SYSTEMATICS_"+SystematicsType+"_"+Scenario+"_"+TString::Format(ShowEPOS?"EPOS":"Normal")+".root";

  TFile* InRootFile = new TFile(InputFolder+InputRootName,"read");
  TNtuple* ntResult = (TNtuple*)InRootFile->Get("ntResult");

  bool* CutUsed = new bool[NumIterResults];
  for(unsigned uBin=0; uBin<NumIterResults; uBin++){
    CutUsed[uBin] = false;
  }

  const TString BranchName[29] = {
      "IterID",//0
      "vCutID",//1
      "vSource",//2
      "vFemReg_pp",//3
      "vFemReg_pL",//4
      "vFemReg_LL",//5
      "vFemReg_pXim",//6
      "vBlReg",//7
      "vMod_pL",//8
      "vFrac_pp_pL",//9
      "vFrac_pL_pSigma0",//10
      "vFrac_pL_pXim",//11
      "vFit",//12
      "vStartPar_LL",//13
      "vSameRad",//14
      "Radius_pp",//15
      "RadiusErr_pp",//16
      "Radius_pL",//17
      "RadiusErr_pL",//18
      "Radius_LL",//19
      "RadiusErr_LL",//20
      "Radius_pXim",//21
      "RadiusErr_pXim",//22
      "a0",//23
      "a0Err",//24
      "rEff",//25
      "rEffErr",//26
      "Chi2Ndf",//27
      "pval"//28
  };
  Float_t ntBuffer[29];
  for(unsigned uBuff=0; uBuff<29; uBuff++){
    ntResult->SetBranchAddress(BranchName[uBuff],&ntBuffer[uBuff]);
  }

  float* Radius = new float [NumIter];
  unsigned* IterID = new unsigned [NumIter];
  //unsigned* EntryID = new unsigned [NumIter];

  const unsigned NumRadPts = 50;
  //const double Rmin = 0.975;
  //const double Rmax = 1.475;
  //const double Rmin = ShowEPOS?1.43:1.15;
  //const double Rmax = ShowEPOS?1.50:1.25;
  //const double Rmin = ShowEPOS?1.44:1.15;//13 TeV new
  //const double Rmax = ShowEPOS?1.52:1.25;
  const double Rmin = ShowEPOS?1.43:1.2;//pPb 5 TeV
  const double Rmax = ShowEPOS?1.50:1.5;

  TH1F* hRadius = new TH1F("hRadius", "hRadius", NumRadPts, Rmin, Rmax);
  unsigned AcceptedNumIter=0;
  for(unsigned uIter=MinIter; uIter<MaxIter; uIter++){

    ntResult->GetEntry(uIter);
    //if(ntBuffer[2]==3)
    //!fit conditions
    //if(ntBuffer[12]!=0) continue; //only norm
    if(ntBuffer[12]!=1) continue; //only 1st order
    //if(ntBuffer[12]!=2) continue; //only 2nd order

    //if(ntBuffer[12]==0) continue; //everything BUT 0th order
    //if(ntBuffer[12]==1) continue; //everything BUT 1th order
    //if(ntBuffer[12]==2) continue; //everything BUT 2nd order

    //if(ntBuffer[8]!=0) continue; //exact NLO
    //if(ntBuffer[8]!=2) continue; //Ledni LO

    //if(ntBuffer[3]==0) continue; //exclude the small femto range in pp
    //if(ntBuffer[4]==0) continue; //exclude the small femto range in pL
    //if(ntBuffer[5]==0) continue; //exclude the small femto range in LL
    //if(ntBuffer[6]==0) continue; //exclude the small femto range in pXim

    //if(ntBuffer[8]!=2) continue;//LO

    if(ntBuffer[2]!=0 && ShowEPOS==false) continue;//Gauss
    if(ntBuffer[2]!=3 && ShowEPOS==true) continue;//EPOS

    //if(ntBuffer[0]>=2000) continue;

    //if(ntBuffer[27]<0.05) continue;

    //if(CutUsed[int(round(ntBuffer[1]))]==true) continue;

    //!CONDITIONS
    if(Scenario=="Global_Radius"){
      if( ntBuffer[14]!=1 ) continue;
      IterID[AcceptedNumIter] = ntBuffer[0];
      Radius[AcceptedNumIter] = ntBuffer[15];
      hRadius->Fill(Radius[AcceptedNumIter]);
      AcceptedNumIter++;
      CutUsed[int(round(ntBuffer[1]))]=true;
    }
    else if(Scenario=="pp_Radius"){
      if( ntBuffer[14]!=0 ) continue;
      IterID[AcceptedNumIter] = ntBuffer[0];
      Radius[AcceptedNumIter] = ntBuffer[15];
      hRadius->Fill(Radius[AcceptedNumIter]);
      AcceptedNumIter++;
      CutUsed[int(round(ntBuffer[1]))]=true;
    }
    else if(Scenario=="pL_Radius"){
      if( ntBuffer[14]!=0 ) continue;
      IterID[AcceptedNumIter] = ntBuffer[0];
      Radius[AcceptedNumIter] = ntBuffer[17];
      hRadius->Fill(Radius[AcceptedNumIter]);
      AcceptedNumIter++;
      CutUsed[int(round(ntBuffer[1]))]=true;
    }
    else{
      continue;
    }






  }
  printf("AcceptedNumIter=%u\n",AcceptedNumIter);

  //sort the radius
  DLM_MergeSort < float, unsigned > SortTool;
  SortTool.SetData(Radius,AcceptedNumIter);
  SortTool.MergeSort();
  SortTool.GetSortedData(Radius,Radius);

  //resort the IterID
  unsigned* Temp_IterID = new unsigned[AcceptedNumIter];
  for(unsigned uEl=0; uEl<AcceptedNumIter; uEl++){
    //printf("uEl=%u --> r=%f\n",uEl,Radius[uEl]);
    Temp_IterID[uEl] = IterID[SortTool.GetKey()[uEl]];
  }
  for(unsigned uEl=0; uEl<AcceptedNumIter; uEl++){
    IterID[uEl] = Temp_IterID[uEl];
  }
  delete [] Temp_IterID;

  //if we have an even number of iterations, there will be two candidates for the median
  unsigned MedianLow = AcceptedNumIter%2==0?AcceptedNumIter/2-1:AcceptedNumIter/2;
  unsigned MedianUp = AcceptedNumIter/2;

  //number of entries to remove on each side to form the final result
  unsigned Central68=floor(double(AcceptedNumIter)*0.16);
  unsigned Central95=floor(double(AcceptedNumIter)*0.025);

  TFile* OutputFile = new TFile(OutFileName, "recreate");

  TFile* FileDefault = new TFile(InputFolder+InputGraphBase+TString::Format("%u.root",0),"read");
  TF1* FitDefault = FileDefault?(TF1*)FileDefault->Get(TString::Format("GlobalFit_%u",0)):nullptr;
  if(FitDefault){
    TGraph* ppGraphDefault = (TGraph*)FileDefault->Get(TString::Format("FitResult_pp_%u",0));
    TGraph* pLamGraphDefault = (TGraph*)FileDefault->Get(TString::Format("FitResult_pL_%u",0));
    TGraph* pLamGraphDefault_LO = (TGraph*)FileDefault->Get(TString::Format("FitResult_pL_LO_%u",0));
    TGraph* pLamGraphDefault_NLO = (TGraph*)FileDefault->Get(TString::Format("FitResult_pL_NLO_%u",0));
    TGraph* LamLamGraphDefault = (TGraph*)FileDefault->Get(TString::Format("FitResult_LL_%u",0));
    TGraph* LamLamGraphDefault_STAR = (TGraph*)FileDefault->Get(TString::Format("FitResult_LL_STAR_%u",0));
    TGraph* pXimGraphDefault = (TGraph*)FileDefault->Get(TString::Format("FitResult_pXim_%u",0));
    TGraph* pXimGraphDefault_COULOMB = (TGraph*)FileDefault->Get(TString::Format("FitResult_pXim_COULOMB_%u",0));
    TF1 Copy_FitDefault(*FitDefault);
    Copy_FitDefault.SetName("FitDefault");
    TGraph Copy_ppGraphDefault(*ppGraphDefault);
    Copy_ppGraphDefault.SetName("ppGraphDefault");
    TGraph Copy_pLamGraphDefault(*pLamGraphDefault);
    Copy_pLamGraphDefault.SetName("pLamGraphDefault");
    TGraph Copy_pLamGraphDefault_LO(*pLamGraphDefault_LO);
    Copy_pLamGraphDefault_LO.SetName("pLamGraphDefault_LO");
    TGraph Copy_pLamGraphDefault_NLO(*pLamGraphDefault_NLO);
    Copy_pLamGraphDefault_NLO.SetName("pLamGraphDefault_NLO");
    TGraph Copy_LamLamGraphDefault(*LamLamGraphDefault);
    Copy_LamLamGraphDefault.SetName("LamLamGraphDefault");
    TGraph Copy_LamLamGraphDefault_STAR(*LamLamGraphDefault_STAR);
    Copy_LamLamGraphDefault_STAR.SetName("LamLamGraphDefault_STAR");
    TGraph Copy_pXimGraphDefault(*pXimGraphDefault);
    Copy_pXimGraphDefault.SetName("pXimGraphDefault");
    TGraph Copy_pXimGraphDefault_COULOMB(*pXimGraphDefault_COULOMB);
    Copy_pXimGraphDefault_COULOMB.SetName("pXimGraphDefault_COULOMB");
    OutputFile->cd();
    Copy_FitDefault.Write();
    Copy_ppGraphDefault.Write();
    Copy_pLamGraphDefault.Write();
    Copy_pLamGraphDefault_LO.Write();
    Copy_pLamGraphDefault_NLO.Write();
    Copy_LamLamGraphDefault.Write();
    Copy_LamLamGraphDefault_STAR.Write();
    Copy_pXimGraphDefault.Write();
    Copy_pXimGraphDefault_COULOMB.Write();
  }

  TFile* FileMedianLow = new TFile(InputFolder+InputGraphBase+TString::Format("%u.root",IterID[MedianLow]),"read");
  TF1* FitMedianLow = (TF1*)FileMedianLow->Get(TString::Format("GlobalFit_%u",IterID[MedianLow]));
  TGraph* ppGraphMedianLow = (TGraph*)FileMedianLow->Get(TString::Format("FitResult_pp_%u",IterID[MedianLow]));
  TGraph* pLamGraphMedianLow = (TGraph*)FileMedianLow->Get(TString::Format("FitResult_pL_%u",IterID[MedianLow]));
  TGraph* pLamGraphMedianLow_LO = (TGraph*)FileMedianLow->Get(TString::Format("FitResult_pL_LO_%u",IterID[MedianLow]));
  TGraph* pLamGraphMedianLow_NLO = (TGraph*)FileMedianLow->Get(TString::Format("FitResult_pL_NLO_%u",IterID[MedianLow]));
  TGraph* LamLamGraphMedianLow = (TGraph*)FileMedianLow->Get(TString::Format("FitResult_LL_%u",IterID[MedianLow]));
  TGraph* LamLamGraphMedianLow_STAR = (TGraph*)FileMedianLow->Get(TString::Format("FitResult_LL_STAR_%u",IterID[MedianLow]));
  TGraph* pXimGraphMedianLow = (TGraph*)FileMedianLow->Get(TString::Format("FitResult_pXim_%u",IterID[MedianLow]));
  TGraph* pXimGraphMedianLow_COULOMB = (TGraph*)FileMedianLow->Get(TString::Format("FitResult_pXim_COULOMB_%u",IterID[MedianLow]));
  TF1 Copy_FitMedianLow(*FitMedianLow);
  Copy_FitMedianLow.SetName("FitMedianLow");
  TGraph Copy_ppGraphMedianLow(*ppGraphMedianLow);
  Copy_ppGraphMedianLow.SetName("ppGraphMedianLow");
  TGraph Copy_pLamGraphMedianLow(*pLamGraphMedianLow);
  Copy_pLamGraphMedianLow.SetName("pLamGraphMedianLow");
  TGraph Copy_pLamGraphMedianLow_LO(*pLamGraphMedianLow_LO);
  Copy_pLamGraphMedianLow_LO.SetName("pLamGraphMedianLow_LO");
  TGraph Copy_pLamGraphMedianLow_NLO(*pLamGraphMedianLow_NLO);
  Copy_pLamGraphMedianLow_NLO.SetName("pLamGraphMedianLow_NLO");
  TGraph Copy_LamLamGraphMedianLow(*LamLamGraphMedianLow);
  Copy_LamLamGraphMedianLow.SetName("LamLamGraphMedianLow");
  TGraph Copy_LamLamGraphMedianLow_STAR(*LamLamGraphMedianLow_STAR);
  Copy_LamLamGraphMedianLow_STAR.SetName("LamLamGraphMedianLow_STAR");
  TGraph Copy_pXimGraphMedianLow(*pXimGraphMedianLow);
  Copy_pXimGraphMedianLow.SetName("pXimGraphMedianLow");
  TGraph Copy_pXimGraphMedianLow_COULOMB(*pXimGraphMedianLow_COULOMB);
  Copy_pXimGraphMedianLow_COULOMB.SetName("pXimGraphMedianLow_COULOMB");
  OutputFile->cd();
  Copy_FitMedianLow.Write();
  Copy_ppGraphMedianLow.Write();
  Copy_pLamGraphMedianLow.Write();
  Copy_pLamGraphMedianLow_LO.Write();
  Copy_pLamGraphMedianLow_NLO.Write();
  Copy_LamLamGraphMedianLow.Write();
  Copy_LamLamGraphMedianLow_STAR.Write();
  Copy_pXimGraphMedianLow.Write();
  Copy_pXimGraphMedianLow_COULOMB.Write();

  TFile* FileMedianUp = new TFile(InputFolder+InputGraphBase+TString::Format("%u.root",IterID[MedianUp]),"read");
  TF1* FitMedianUp = (TF1*)FileMedianUp->Get(TString::Format("GlobalFit_%u",IterID[MedianUp]));
  TGraph* ppGraphMedianUp = (TGraph*)FileMedianUp->Get(TString::Format("FitResult_pp_%u",IterID[MedianUp]));
  TGraph* pLamGraphMedianUp = (TGraph*)FileMedianUp->Get(TString::Format("FitResult_pL_%u",IterID[MedianUp]));
  TGraph* pLamGraphMedianUp_LO = (TGraph*)FileMedianUp->Get(TString::Format("FitResult_pL_LO_%u",IterID[MedianUp]));
  TGraph* pLamGraphMedianUp_NLO = (TGraph*)FileMedianUp->Get(TString::Format("FitResult_pL_NLO_%u",IterID[MedianUp]));
  TGraph* LamLamGraphMedianUp = (TGraph*)FileMedianUp->Get(TString::Format("FitResult_LL_%u",IterID[MedianUp]));
  TGraph* LamLamGraphMedianUp_STAR = (TGraph*)FileMedianUp->Get(TString::Format("FitResult_LL_STAR_%u",IterID[MedianUp]));
  TGraph* pXimGraphMedianUp = (TGraph*)FileMedianUp->Get(TString::Format("FitResult_pXim_%u",IterID[MedianUp]));
  TGraph* pXimGraphMedianUp_COULOMB = (TGraph*)FileMedianUp->Get(TString::Format("FitResult_pXim_COULOMB_%u",IterID[MedianUp]));
  TF1 Copy_FitMedianUp(*FitMedianUp);
  Copy_FitMedianUp.SetName("FitMedianUp");
  TGraph Copy_ppGraphMedianUp(*ppGraphMedianUp);
  Copy_ppGraphMedianUp.SetName("ppGraphMedianUp");
  TGraph Copy_pLamGraphMedianUp(*pLamGraphMedianUp);
  Copy_pLamGraphMedianUp.SetName("pLamGraphMedianUp");
  TGraph Copy_pLamGraphMedianUp_LO(*pLamGraphMedianUp_LO);
  Copy_pLamGraphMedianUp_LO.SetName("pLamGraphMedianUp_LO");
  TGraph Copy_pLamGraphMedianUp_NLO(*pLamGraphMedianUp_NLO);
  Copy_pLamGraphMedianUp_NLO.SetName("pLamGraphMedianUp_NLO");
  TGraph Copy_LamLamGraphMedianUp(*LamLamGraphMedianUp);
  Copy_LamLamGraphMedianUp.SetName("LamLamGraphMedianUp");
  TGraph Copy_LamLamGraphMedianUp_STAR(*LamLamGraphMedianUp_STAR);
  Copy_LamLamGraphMedianUp_STAR.SetName("LamLamGraphMedianUp_STAR");
  TGraph Copy_pXimGraphMedianUp(*pXimGraphMedianUp);
  Copy_pXimGraphMedianUp.SetName("pXimGraphMedianUp");
  TGraph Copy_pXimGraphMedianUp_COULOMB(*pXimGraphMedianUp_COULOMB);
  Copy_pXimGraphMedianUp_COULOMB.SetName("pXimGraphMedianUp_COULOMB");
  OutputFile->cd();
  Copy_FitMedianUp.Write();
  Copy_ppGraphMedianUp.Write();
  Copy_pLamGraphMedianUp.Write();
  Copy_pLamGraphMedianUp_LO.Write();
  Copy_pLamGraphMedianUp_NLO.Write();
  Copy_LamLamGraphMedianUp.Write();
  Copy_LamLamGraphMedianUp_STAR.Write();
  Copy_pXimGraphMedianUp.Write();
  Copy_pXimGraphMedianUp_COULOMB.Write();


  TFile* FileLowerLim = new TFile(InputFolder+InputGraphBase+TString::Format("%u.root",IterID[Central68]),"read");
  TF1* FitLowerLim = (TF1*)FileLowerLim->Get(TString::Format("GlobalFit_%u",IterID[Central68]));
  TGraph* ppGraphLowerLim = (TGraph*)FileLowerLim->Get(TString::Format("FitResult_pp_%u",IterID[Central68]));
  TGraph* pLamGraphLowerLim = (TGraph*)FileLowerLim->Get(TString::Format("FitResult_pL_%u",IterID[Central68]));
  TGraph* pLamGraphLowerLim_LO = (TGraph*)FileLowerLim->Get(TString::Format("FitResult_pL_LO_%u",IterID[Central68]));
  TGraph* pLamGraphLowerLim_NLO = (TGraph*)FileLowerLim->Get(TString::Format("FitResult_pL_NLO_%u",IterID[Central68]));
  TGraph* LamLamGraphLowerLim = (TGraph*)FileLowerLim->Get(TString::Format("FitResult_LL_%u",IterID[Central68]));
  TGraph* LamLamGraphLowerLim_STAR = (TGraph*)FileLowerLim->Get(TString::Format("FitResult_LL_STAR_%u",IterID[Central68]));
  TGraph* pXimGraphLowerLim = (TGraph*)FileLowerLim->Get(TString::Format("FitResult_pXim_%u",IterID[Central68]));
  TGraph* pXimGraphLowerLim_COULOMB = (TGraph*)FileLowerLim->Get(TString::Format("FitResult_pXim_COULOMB_%u",IterID[Central68]));
  TF1 Copy_FitLowerLim(*FitLowerLim);
  Copy_FitLowerLim.SetName("FitLowerLim");
  TGraph Copy_ppGraphLowerLim(*ppGraphLowerLim);
  Copy_ppGraphLowerLim.SetName("ppGraphLowerLim");
  TGraph Copy_pLamGraphLowerLim(*pLamGraphLowerLim);
  Copy_pLamGraphLowerLim.SetName("pLamGraphLowerLim");
  TGraph Copy_pLamGraphLowerLim_LO(*pLamGraphLowerLim_LO);
  Copy_pLamGraphLowerLim_LO.SetName("Copy_pLamGraphLowerLim_LO");
  TGraph Copy_pLamGraphLowerLim_NLO(*pLamGraphLowerLim_NLO);
  Copy_pLamGraphLowerLim_NLO.SetName("Copy_pLamGraphLowerLim_NLO");
  TGraph Copy_LamLamGraphLowerLim(*LamLamGraphLowerLim);
  Copy_LamLamGraphLowerLim.SetName("LamLamGraphLowerLim");
  TGraph Copy_LamLamGraphLowerLim_STAR(*LamLamGraphLowerLim_STAR);
  Copy_LamLamGraphLowerLim_STAR.SetName("LamLamGraphLowerLim_STAR");
  TGraph Copy_pXimGraphLowerLim(*pXimGraphLowerLim);
  Copy_pXimGraphLowerLim.SetName("pXimGraphLowerLim");
  TGraph Copy_pXimGraphLowerLim_COULOMB(*pXimGraphLowerLim_COULOMB);
  Copy_pXimGraphLowerLim_COULOMB.SetName("pXimGraphLowerLim_COULOMB");
  OutputFile->cd();
  Copy_FitLowerLim.Write();
  Copy_ppGraphLowerLim.Write();
  Copy_pLamGraphLowerLim.Write();
  Copy_pLamGraphLowerLim_LO.Write();
  Copy_pLamGraphLowerLim_NLO.Write();
  Copy_LamLamGraphLowerLim.Write();
  Copy_LamLamGraphLowerLim_STAR.Write();
  Copy_pXimGraphLowerLim.Write();
  Copy_pXimGraphLowerLim_COULOMB.Write();
  std::cout << "Ahhhhh:\t" << InputFolder+InputGraphBase+TString::Format("%u.root",IterID[AcceptedNumIter-1-Central68]) << std::endl;
  TFile* FileUpperLim = new TFile(InputFolder+InputGraphBase+TString::Format("%u.root",IterID[AcceptedNumIter-1-Central68]),"read");
  TF1* FitUpperLim = (TF1*)FileUpperLim->Get(TString::Format("GlobalFit_%u",IterID[AcceptedNumIter-1-Central68]));
  TGraph* ppGraphUpperLim = (TGraph*)FileUpperLim->Get(TString::Format("FitResult_pp_%u",IterID[AcceptedNumIter-1-Central68]));
  TGraph* pLamGraphUpperLim = (TGraph*)FileUpperLim->Get(TString::Format("FitResult_pL_%u",IterID[AcceptedNumIter-1-Central68]));
  TGraph* pLamGraphUpperLim_LO = (TGraph*)FileUpperLim->Get(TString::Format("FitResult_pL_LO_%u",IterID[AcceptedNumIter-1-Central68]));
  TGraph* pLamGraphUpperLim_NLO = (TGraph*)FileUpperLim->Get(TString::Format("FitResult_pL_NLO_%u",IterID[AcceptedNumIter-1-Central68]));
  TGraph* LamLamGraphUpperLim = (TGraph*)FileUpperLim->Get(TString::Format("FitResult_LL_%u",IterID[AcceptedNumIter-1-Central68]));
  TGraph* LamLamGraphUpperLim_STAR = (TGraph*)FileUpperLim->Get(TString::Format("FitResult_LL_STAR_%u",IterID[AcceptedNumIter-1-Central68]));
  TGraph* pXimGraphUpperLim = (TGraph*)FileUpperLim->Get(TString::Format("FitResult_pXim_%u",IterID[AcceptedNumIter-1-Central68]));
  TGraph* pXimGraphUpperLim_COULOMB = (TGraph*)FileUpperLim->Get(TString::Format("FitResult_pXim_COULOMB_%u",IterID[AcceptedNumIter-1-Central68]));
  TF1 Copy_FitUpperLim(*FitUpperLim);
  Copy_FitUpperLim.SetName("FitUpperLim");
  TGraph Copy_ppGraphUpperLim(*ppGraphUpperLim);
  Copy_ppGraphUpperLim.SetName("ppGraphUpperLim");
  TGraph Copy_pLamGraphUpperLim(*pLamGraphUpperLim);
  Copy_pLamGraphUpperLim.SetName("pLamGraphUpperLim");
  TGraph Copy_pLamGraphUpperLim_LO(*pLamGraphUpperLim_LO);
  Copy_pLamGraphUpperLim_LO.SetName("pLamGraphUpperLim_LO");
  TGraph Copy_pLamGraphUpperLim_NLO(*pLamGraphUpperLim_NLO);
  Copy_pLamGraphUpperLim_NLO.SetName("pLamGraphUpperLim_NLO");
  TGraph Copy_LamLamGraphUpperLim(*LamLamGraphUpperLim);
  Copy_LamLamGraphUpperLim.SetName("LamLamGraphUpperLim");
  TGraph Copy_LamLamGraphUpperLim_STAR(*LamLamGraphUpperLim_STAR);
  Copy_LamLamGraphUpperLim_STAR.SetName("LamLamGraphUpperLim_STAR");
  TGraph Copy_pXimGraphUpperLim(*pXimGraphUpperLim);
  Copy_pXimGraphUpperLim.SetName("pXimGraphUpperLim");
  TGraph Copy_pXimGraphUpperLim_COULOMB(*pXimGraphUpperLim_COULOMB);
  Copy_pXimGraphUpperLim_COULOMB.SetName("pXimGraphUpperLim_COULOMB");
  OutputFile->cd();
  Copy_FitUpperLim.Write();
  Copy_ppGraphUpperLim.Write();
  Copy_pLamGraphUpperLim.Write();
  Copy_pLamGraphUpperLim_LO.Write();
  Copy_pLamGraphUpperLim_NLO.Write();
  Copy_LamLamGraphUpperLim.Write();
  Copy_LamLamGraphUpperLim_STAR.Write();
  Copy_pXimGraphUpperLim.Write();
  Copy_pXimGraphUpperLim_COULOMB.Write();

  //TH1F* hLimits = new TH1F("hLimits", "hLimits", 1, Radius[Central68], Radius[AcceptedNumIter-1-Central68]);
  const unsigned LimitLowBin = hRadius->FindBin(Radius[Central68]);
  const unsigned LimitLowUp = hRadius->FindBin(Radius[AcceptedNumIter-1-Central68]);
  TH1F* hOuterLow = new TH1F("hOuterLow", "hOuterLow", LimitLowBin-1, Rmin, hRadius->GetBinLowEdge(LimitLowBin));
  TH1F* hInner = new TH1F("hInner", "hInner", LimitLowUp-LimitLowBin+1, hRadius->GetBinLowEdge(LimitLowBin), hRadius->GetXaxis()->GetBinUpEdge(LimitLowUp));
  TH1F* hOuterUp = new TH1F("hOuterUp", "hOuterUp", NumRadPts-LimitLowUp, hRadius->GetXaxis()->GetBinUpEdge(LimitLowUp), Rmax);
  unsigned CurrentBin=0;
  for(unsigned uBin=0; uBin<LimitLowBin-1; uBin++){
    hOuterLow->SetBinContent(uBin+1, hRadius->GetBinContent(CurrentBin));
    CurrentBin++;
  }
  //printf("hi!\n");
  for(unsigned uBin=0; uBin<LimitLowUp-LimitLowBin+1; uBin++){
    hInner->SetBinContent(uBin+1, hRadius->GetBinContent(CurrentBin));
    CurrentBin++;
  }
  //printf("hi!\n");
  for(unsigned uBin=0; uBin<NumRadPts-LimitLowUp; uBin++){
    hOuterUp->SetBinContent(uBin+1, hRadius->GetBinContent(CurrentBin));
    CurrentBin++;
  }
  printf("NumRadPts=%u =? CurrentBin=%u\n", NumRadPts, CurrentBin);

  hRadius->SetStats(false);
  hRadius->SetTitle("");
  hRadius->GetXaxis()->SetLabelSize(0.065);
  hRadius->GetXaxis()->SetTitle(ShowEPOS?"N_{R}":"r (fm)");
  hRadius->GetXaxis()->CenterTitle();
  hRadius->GetXaxis()->SetTitleOffset(1.15);
  hRadius->GetXaxis()->SetLabelOffset(0.02);
  hRadius->GetXaxis()->SetTitleSize(0.075);
  hRadius->GetYaxis()->SetLabelSize(0.065);
  hRadius->GetYaxis()->SetTitle("Counts");
  hRadius->GetYaxis()->CenterTitle();
  hRadius->GetYaxis()->SetTitleOffset(0.7);
  hRadius->GetYaxis()->SetTitleSize(0.075);
  hRadius->GetYaxis()->SetLimits(0, hRadius->GetBinContent(hRadius->GetMaximumBin())*1.2);
  hRadius->GetYaxis()->SetRangeUser(0, hRadius->GetBinContent(hRadius->GetMaximumBin())*1.2);

  hOuterLow->SetLineColor(kGray);
  hOuterLow->SetFillColor(kGray);
  hInner->SetLineColor(kAzure+1);
  hInner->SetFillColor(kAzure+1);
  hOuterUp->SetLineColor(kGray);
  hOuterUp->SetFillColor(kGray);

  const unsigned NumFitPars = 13;
  printf("Extended info for the default:\n");
  ntResult->GetEntry(0);
  for(unsigned uBuff=0; uBuff<29; uBuff++){
    if(uBuff<14) printf(" %s = %.0f\n", BranchName[uBuff].Data(), ntBuffer[uBuff]);
    else printf(" %s = %.3f\n", BranchName[uBuff].Data(), ntBuffer[uBuff]);
  }
  if(FitDefault){
    printf(" BL_a(pp) = %.3f +/- %.3f\n", FitDefault->GetParameter(DLM_Fitter1::p_a), FitDefault->GetParError(DLM_Fitter1::p_a));
    printf(" BL_b(pp) = %.3f +/- %.3f\n", FitDefault->GetParameter(DLM_Fitter1::p_b)*1000, FitDefault->GetParError(DLM_Fitter1::p_b)*1000);
    printf(" Cl(pp) = %.3f +/- %.3f\n", FitDefault->GetParameter(DLM_Fitter1::p_Cl), FitDefault->GetParError(DLM_Fitter1::p_Cl));
    printf(" BL_a(pL) = %.3f +/- %.3f\n", FitDefault->GetParameter(DLM_Fitter1::p_a+NumFitPars), FitDefault->GetParError(DLM_Fitter1::p_a+NumFitPars));
    printf(" BL_b(pL) = %.3f +/- %.3f\n", FitDefault->GetParameter(DLM_Fitter1::p_b+NumFitPars)*1000, FitDefault->GetParError(DLM_Fitter1::p_b+NumFitPars)*1000);
    printf(" Cl(pL) = %.3f +/- %.3f\n", FitDefault->GetParameter(DLM_Fitter1::p_Cl+NumFitPars), FitDefault->GetParError(DLM_Fitter1::p_Cl+NumFitPars));
    printf(" BL_a(LL) = %.3f +/- %.3f\n", FitDefault->GetParameter(DLM_Fitter1::p_a+2*NumFitPars), FitDefault->GetParError(DLM_Fitter1::p_a+2*NumFitPars));
    printf(" BL_b(LL) = %.3f +/- %.3f\n", FitDefault->GetParameter(DLM_Fitter1::p_b+2*NumFitPars)*1000, FitDefault->GetParError(DLM_Fitter1::p_b+2*NumFitPars)*1000);
    printf(" Cl(LL) = %.3f +/- %.3f\n", FitDefault->GetParameter(DLM_Fitter1::p_Cl+2*NumFitPars), FitDefault->GetParError(DLM_Fitter1::p_Cl+2*NumFitPars));
    printf(" BL_a(pXi) = %.3f +/- %.3f\n", FitDefault->GetParameter(DLM_Fitter1::p_a+3*NumFitPars), FitDefault->GetParError(DLM_Fitter1::p_a+3*NumFitPars));
    printf(" BL_b(pXi) = %.3f +/- %.3f\n", FitDefault->GetParameter(DLM_Fitter1::p_b+3*NumFitPars)*1000, FitDefault->GetParError(DLM_Fitter1::p_b+3*NumFitPars)*1000);
    printf(" Cl(pXi) = %.3f +/- %.3f\n", FitDefault->GetParameter(DLM_Fitter1::p_Cl+3*NumFitPars), FitDefault->GetParError(DLM_Fitter1::p_Cl+3*NumFitPars));

    DefaultRad = FitDefault->GetParameter(DLM_Fitter1::p_sor0);
  }
  //printf("DefaultRad=%f\n",DefaultRad);return;
  //TLine PublishedRadius(1.159,0,1.159,hRadius->GetBinContent(hRadius->GetMaximumBin())*1.2);
  TLine PublishedRadius(DefaultRad,0,DefaultRad,hRadius->GetBinContent(hRadius->GetMaximumBin())*1.2);
  PublishedRadius.SetLineColor(kRed);
  PublishedRadius.SetLineWidth(8);

  TPaveText* myPT = new TPaveText(0.7,0.7,0.975,0.975, "blNDC");//lbrt
  myPT->SetName("myPT");
  myPT->SetBorderSize(1);
  myPT->SetTextSize(0.045);
  myPT->SetFillColor(kWhite);
  myPT->SetTextFont(22);
  if(FitDefault) myPT->AddText(TString::Format("#color[%u]{Default     R = %.3f%s}",kRed,DefaultRad,ShowEPOS?"":" fm"));
  else myPT->AddText(TString::Format(" "));

  //myPT->AddText(TString::Format("--------------"));
  if(ShowEPOS){
    myPT->AddText(TString::Format("#color[%u]{Minimum N_{R} = %.3f%s}", kAzure+1, Radius[Central68],ShowEPOS?"":" fm"));
    myPT->AddText(TString::Format("#color[%u]{Median     N_{R} = %.3f%s}", kAzure+1, Radius[MedianLow],ShowEPOS?"":" fm"));
    myPT->AddText(TString::Format("#color[%u]{Maximum N_{R} = %.3f%s}", kAzure+1, Radius[AcceptedNumIter-1-Central68],ShowEPOS?"":" fm"));
  }
  else{
    myPT->AddText(TString::Format("#color[%u]{Minimum R = %.3f%s}", kAzure+1, Radius[Central68],ShowEPOS?"":" fm"));
    myPT->AddText(TString::Format("#color[%u]{Median     R = %.3f%s}", kAzure+1, Radius[MedianLow],ShowEPOS?"":" fm"));
    myPT->AddText(TString::Format("#color[%u]{Maximum R = %.3f%s}", kAzure+1, Radius[AcceptedNumIter-1-Central68],ShowEPOS?"":" fm"));
  }


  TCanvas* cRadius = new TCanvas("cRadius", "cRadius", 1);
  cRadius->cd(0); cRadius->SetCanvasSize(1920, 1080); cRadius->SetMargin(0.12,0.05,0.18,0.05);//lrbt
  hRadius->Draw("axis");
  //hRadius->Draw("");
  hInner->Draw("same");
  hOuterLow->Draw("same");
  hOuterUp->Draw("same");
  PublishedRadius.Draw("same");
  myPT->Draw("same");
  cRadius->SaveAs(TString::Format("%scRadius_%s_%s_%s.png",OutDirName.Data(),SystematicsType.Data(),Scenario.Data(),ShowEPOS?"EPOS":"Normal"));
  PublishedRadius.SetLineWidth(PublishedRadius.GetLineWidth()/2.5);
  cRadius->SaveAs(TString::Format("%scRadius_%s_%s_%s.pdf",OutDirName.Data(),SystematicsType.Data(),Scenario.Data(),ShowEPOS?"EPOS":"Normal"));
  PublishedRadius.SetLineWidth(PublishedRadius.GetLineWidth()*2.5);

  printf("-----------------------\n");
  printf("Extended info for the first median:\n");
  ntResult->GetEntry(IterID[MedianLow]);
  for(unsigned uBuff=0; uBuff<29; uBuff++){
    if(uBuff<14) printf(" %s = %.0f\n", BranchName[uBuff].Data(), ntBuffer[uBuff]);
    else printf(" %s = %.3f\n", BranchName[uBuff].Data(), ntBuffer[uBuff]);
  }

  printf(" BL_a(pp) = %.3f +/- %.3f\n", FitMedianLow->GetParameter(DLM_Fitter1::p_a), FitMedianLow->GetParError(DLM_Fitter1::p_a));
  printf(" BL_b(pp) = %.3f +/- %.3f\n", FitMedianLow->GetParameter(DLM_Fitter1::p_b)*1000, FitMedianLow->GetParError(DLM_Fitter1::p_b)*1000);
  printf(" Cl(pp) = %.3f +/- %.3f\n", FitMedianLow->GetParameter(DLM_Fitter1::p_Cl), FitMedianLow->GetParError(DLM_Fitter1::p_Cl));
  printf(" BL_a(pL) = %.3f +/- %.3f\n", FitMedianLow->GetParameter(DLM_Fitter1::p_a+NumFitPars), FitMedianLow->GetParError(DLM_Fitter1::p_a+NumFitPars));
  printf(" BL_b(pL) = %.3f +/- %.3f\n", FitMedianLow->GetParameter(DLM_Fitter1::p_b+NumFitPars)*1000, FitMedianLow->GetParError(DLM_Fitter1::p_b+NumFitPars)*1000);
  printf(" Cl(pL) = %.3f +/- %.3f\n", FitMedianLow->GetParameter(DLM_Fitter1::p_Cl+NumFitPars), FitMedianLow->GetParError(DLM_Fitter1::p_Cl+NumFitPars));
  printf(" BL_a(LL) = %.3f +/- %.3f\n", FitMedianLow->GetParameter(DLM_Fitter1::p_a+2*NumFitPars), FitMedianLow->GetParError(DLM_Fitter1::p_a+2*NumFitPars));
  printf(" BL_b(LL) = %.3f +/- %.3f\n", FitMedianLow->GetParameter(DLM_Fitter1::p_b+2*NumFitPars)*1000, FitMedianLow->GetParError(DLM_Fitter1::p_b+2*NumFitPars)*1000);
  printf(" Cl(LL) = %.3f +/- %.3f\n", FitMedianLow->GetParameter(DLM_Fitter1::p_Cl+2*NumFitPars), FitMedianLow->GetParError(DLM_Fitter1::p_Cl+2*NumFitPars));
  printf(" BL_a(pXi) = %.3f +/- %.3f\n", FitMedianLow->GetParameter(DLM_Fitter1::p_a+3*NumFitPars), FitMedianLow->GetParError(DLM_Fitter1::p_a+3*NumFitPars));
  printf(" BL_b(pXi) = %.3f +/- %.3f\n", FitMedianLow->GetParameter(DLM_Fitter1::p_b+3*NumFitPars)*1000, FitMedianLow->GetParError(DLM_Fitter1::p_b+3*NumFitPars)*1000);
  printf(" Cl(pXi) = %.3f +/- %.3f\n", FitMedianLow->GetParameter(DLM_Fitter1::p_Cl+3*NumFitPars), FitMedianLow->GetParError(DLM_Fitter1::p_Cl+3*NumFitPars));

  printf("-----------------------\n");
  printf("Extended info for the second median:\n");
  ntResult->GetEntry(IterID[MedianUp]);
  for(unsigned uBuff=0; uBuff<29; uBuff++){
    if(uBuff<14) printf(" %s = %.0f\n", BranchName[uBuff].Data(), ntBuffer[uBuff]);
    else printf(" %s = %.3f\n", BranchName[uBuff].Data(), ntBuffer[uBuff]);
  }
  printf("Extended info for the 1σ lower limit:\n");
  ntResult->GetEntry(IterID[Central68]);
  for(unsigned uBuff=0; uBuff<29; uBuff++){
    if(uBuff<14) printf(" %s = %.0f\n", BranchName[uBuff].Data(), ntBuffer[uBuff]);
    else printf(" %s = %.3f\n", BranchName[uBuff].Data(), ntBuffer[uBuff]);
  }
  printf(" BL_a(pp) = %.3f +/- %.3f\n", FitMedianUp->GetParameter(DLM_Fitter1::p_a), FitMedianUp->GetParError(DLM_Fitter1::p_a));
  printf(" BL_b(pp) = %.3f +/- %.3f\n", FitMedianUp->GetParameter(DLM_Fitter1::p_b)*1000, FitMedianUp->GetParError(DLM_Fitter1::p_b)*1000);
  printf(" Cl(pp) = %.3f +/- %.3f\n", FitMedianUp->GetParameter(DLM_Fitter1::p_Cl), FitMedianUp->GetParError(DLM_Fitter1::p_Cl));
  printf(" BL_a(pL) = %.3f +/- %.3f\n", FitMedianUp->GetParameter(DLM_Fitter1::p_a+NumFitPars), FitMedianUp->GetParError(DLM_Fitter1::p_a+NumFitPars));
  printf(" BL_b(pL) = %.3f +/- %.3f\n", FitMedianUp->GetParameter(DLM_Fitter1::p_b+NumFitPars)*1000, FitMedianUp->GetParError(DLM_Fitter1::p_b+NumFitPars)*1000);
  printf(" Cl(pL) = %.3f +/- %.3f\n", FitMedianUp->GetParameter(DLM_Fitter1::p_Cl+NumFitPars), FitMedianUp->GetParError(DLM_Fitter1::p_Cl+NumFitPars));
  printf(" BL_a(LL) = %.3f +/- %.3f\n", FitMedianUp->GetParameter(DLM_Fitter1::p_a+2*NumFitPars), FitMedianUp->GetParError(DLM_Fitter1::p_a+2*NumFitPars));
  printf(" BL_b(LL) = %.3f +/- %.3f\n", FitMedianUp->GetParameter(DLM_Fitter1::p_b+2*NumFitPars)*1000, FitMedianUp->GetParError(DLM_Fitter1::p_b+2*NumFitPars)*1000);
  printf(" Cl(LL) = %.3f +/- %.3f\n", FitMedianUp->GetParameter(DLM_Fitter1::p_Cl+2*NumFitPars), FitMedianUp->GetParError(DLM_Fitter1::p_Cl+2*NumFitPars));
  printf(" BL_a(pXi) = %.3f +/- %.3f\n", FitMedianUp->GetParameter(DLM_Fitter1::p_a+3*NumFitPars), FitMedianUp->GetParError(DLM_Fitter1::p_a+3*NumFitPars));
  printf(" BL_b(pXi) = %.3f +/- %.3f\n", FitMedianUp->GetParameter(DLM_Fitter1::p_b+3*NumFitPars)*1000, FitMedianUp->GetParError(DLM_Fitter1::p_b+3*NumFitPars)*1000);
  printf(" Cl(pXi) = %.3f +/- %.3f\n", FitMedianUp->GetParameter(DLM_Fitter1::p_Cl+3*NumFitPars), FitMedianUp->GetParError(DLM_Fitter1::p_Cl+3*NumFitPars));

  printf("Extended info for the 1σ upper limit:\n");
  ntResult->GetEntry(IterID[AcceptedNumIter-1-Central68]);
  for(unsigned uBuff=0; uBuff<29; uBuff++){
    if(uBuff<14) printf(" %s = %.0f\n", BranchName[uBuff].Data(), ntBuffer[uBuff]);
    else printf(" %s = %.3f\n", BranchName[uBuff].Data(), ntBuffer[uBuff]);
  }
  printf(" BL_a(pp) = %.3f +/- %.3f\n", FitUpperLim->GetParameter(DLM_Fitter1::p_a), FitUpperLim->GetParError(DLM_Fitter1::p_a));
  printf(" BL_b(pp) = %.3f +/- %.3f\n", FitUpperLim->GetParameter(DLM_Fitter1::p_b)*1000, FitUpperLim->GetParError(DLM_Fitter1::p_b)*1000);
  printf(" Cl(pp) = %.3f +/- %.3f\n", FitUpperLim->GetParameter(DLM_Fitter1::p_Cl), FitUpperLim->GetParError(DLM_Fitter1::p_Cl));
  printf(" BL_a(pL) = %.3f +/- %.3f\n", FitUpperLim->GetParameter(DLM_Fitter1::p_a+NumFitPars), FitUpperLim->GetParError(DLM_Fitter1::p_a+NumFitPars));
  printf(" BL_b(pL) = %.3f +/- %.3f\n", FitUpperLim->GetParameter(DLM_Fitter1::p_b+NumFitPars)*1000, FitUpperLim->GetParError(DLM_Fitter1::p_b+NumFitPars)*1000);
  printf(" Cl(pL) = %.3f +/- %.3f\n", FitUpperLim->GetParameter(DLM_Fitter1::p_Cl+NumFitPars), FitUpperLim->GetParError(DLM_Fitter1::p_Cl+NumFitPars));
  printf(" BL_a(LL) = %.3f +/- %.3f\n", FitUpperLim->GetParameter(DLM_Fitter1::p_a+2*NumFitPars), FitUpperLim->GetParError(DLM_Fitter1::p_a+2*NumFitPars));
  printf(" BL_b(LL) = %.3f +/- %.3f\n", FitUpperLim->GetParameter(DLM_Fitter1::p_b+2*NumFitPars)*1000, FitUpperLim->GetParError(DLM_Fitter1::p_b+2*NumFitPars)*1000);
  printf(" Cl(LL) = %.3f +/- %.3f\n", FitUpperLim->GetParameter(DLM_Fitter1::p_Cl+2*NumFitPars), FitUpperLim->GetParError(DLM_Fitter1::p_Cl+2*NumFitPars));
  printf(" BL_a(pXi) = %.3f +/- %.3f\n", FitUpperLim->GetParameter(DLM_Fitter1::p_a+3*NumFitPars), FitUpperLim->GetParError(DLM_Fitter1::p_a+3*NumFitPars));
  printf(" BL_b(pXi) = %.3f +/- %.3f\n", FitUpperLim->GetParameter(DLM_Fitter1::p_b+3*NumFitPars)*1000, FitUpperLim->GetParError(DLM_Fitter1::p_b+3*NumFitPars)*1000);
  printf(" Cl(pXi) = %.3f +/- %.3f\n", FitUpperLim->GetParameter(DLM_Fitter1::p_Cl+3*NumFitPars), FitUpperLim->GetParError(DLM_Fitter1::p_Cl+3*NumFitPars));

  printf("Extended info for the 1σ lower limit:\n");
  ntResult->GetEntry(IterID[Central68]);
  for(unsigned uBuff=0; uBuff<29; uBuff++){
    if(uBuff<14) printf(" %s = %.0f\n", BranchName[uBuff].Data(), ntBuffer[uBuff]);
    else printf(" %s = %.3f\n", BranchName[uBuff].Data(), ntBuffer[uBuff]);
  }
  printf(" BL_a(pp) = %.3f +/- %.3f\n", FitLowerLim->GetParameter(DLM_Fitter1::p_a), FitLowerLim->GetParError(DLM_Fitter1::p_a));
  printf(" BL_b(pp) = %.3f +/- %.3f\n", FitLowerLim->GetParameter(DLM_Fitter1::p_b)*1000, FitLowerLim->GetParError(DLM_Fitter1::p_b)*1000);
  printf(" Cl(pp) = %.3f +/- %.3f\n", FitLowerLim->GetParameter(DLM_Fitter1::p_Cl), FitLowerLim->GetParError(DLM_Fitter1::p_Cl));
  printf(" BL_a(pL) = %.3f +/- %.3f\n", FitLowerLim->GetParameter(DLM_Fitter1::p_a+NumFitPars), FitLowerLim->GetParError(DLM_Fitter1::p_a+NumFitPars));
  printf(" BL_b(pL) = %.3f +/- %.3f\n", FitLowerLim->GetParameter(DLM_Fitter1::p_b+NumFitPars)*1000, FitLowerLim->GetParError(DLM_Fitter1::p_b+NumFitPars)*1000);
  printf(" Cl(pL) = %.3f +/- %.3f\n", FitLowerLim->GetParameter(DLM_Fitter1::p_Cl+NumFitPars), FitLowerLim->GetParError(DLM_Fitter1::p_Cl+NumFitPars));
  printf(" BL_a(LL) = %.3f +/- %.3f\n", FitLowerLim->GetParameter(DLM_Fitter1::p_a+2*NumFitPars), FitLowerLim->GetParError(DLM_Fitter1::p_a+2*NumFitPars));
  printf(" BL_b(LL) = %.3f +/- %.3f\n", FitLowerLim->GetParameter(DLM_Fitter1::p_b+2*NumFitPars)*1000, FitLowerLim->GetParError(DLM_Fitter1::p_b+2*NumFitPars)*1000);
  printf(" Cl(LL) = %.3f +/- %.3f\n", FitLowerLim->GetParameter(DLM_Fitter1::p_Cl+2*NumFitPars), FitLowerLim->GetParError(DLM_Fitter1::p_Cl+2*NumFitPars));
  printf(" BL_a(pXi) = %.3f +/- %.3f\n", FitLowerLim->GetParameter(DLM_Fitter1::p_a+3*NumFitPars), FitLowerLim->GetParError(DLM_Fitter1::p_a+3*NumFitPars));
  printf(" BL_b(pXi) = %.3f +/- %.3f\n", FitLowerLim->GetParameter(DLM_Fitter1::p_b+3*NumFitPars)*1000, FitLowerLim->GetParError(DLM_Fitter1::p_b+3*NumFitPars)*1000);
  printf(" Cl(pXi) = %.3f +/- %.3f\n", FitLowerLim->GetParameter(DLM_Fitter1::p_Cl+3*NumFitPars), FitLowerLim->GetParError(DLM_Fitter1::p_Cl+3*NumFitPars));


  printf("-----------------------\n");
  printf("Accepted number of iterations: %u\n",AcceptedNumIter);
  printf("Median-Low: %u\n",MedianLow);
  printf("Median-up: %u\n",MedianUp);
  printf("Lower limit: %u --> %u\n",Central68,Central95);
  printf("Upper limit: %u --> %u\n",AcceptedNumIter-1-Central68,AcceptedNumIter-1-Central95);
  printf("-----------------------\n");
  printf("The median is at:\n");
  printf(" IterID = %u or %u\n", IterID[MedianLow], IterID[MedianUp]);
  printf(" Radius = %f fm or %f fm\n", Radius[MedianLow], Radius[MedianUp]);
  printf("The lower limit is at:\n");
  printf(" IterID = %u(1σ), %u(2σ)\n", IterID[Central68], IterID[Central95]);
  printf(" Radius = %f fm(1σ), %f fm(2σ)\n", Radius[Central68], Radius[Central95]);
  printf("The upper limit is at:\n");
  printf(" IterID = %u(1σ), %u(2σ)\n", IterID[AcceptedNumIter-1-Central68], IterID[AcceptedNumIter-1-Central95]);
  printf(" Radius = %f fm(1σ), %f fm(2σ)\n", Radius[AcceptedNumIter-1-Central68], Radius[AcceptedNumIter-1-Central95]);
  printf("-----------------------\n");
  printf("The 1σ result to quote is:\n");
  printf(" r = %.3f +(%.3f) -(%.3f)\n", Radius[MedianLow],
         Radius[AcceptedNumIter-1-Central68]-Radius[MedianLow],
         Radius[MedianLow]-Radius[Central68]);
  printf(" or\n");
  printf(" r = %.3f +(%.3f) -(%.3f)\n", Radius[MedianUp],
         Radius[AcceptedNumIter-1-Central68]-Radius[MedianUp],
         Radius[MedianUp]-Radius[Central68]);
  printf("The 2σ result to quote is:\n");
  printf(" r = %.3f +(%.3f) -(%.3f)\n", Radius[MedianLow],
         Radius[AcceptedNumIter-1-Central95]-Radius[MedianLow],
         Radius[MedianLow]-Radius[Central95]);
  printf(" or\n");
  printf(" r = %.3f +(%.3f) -(%.3f)\n", Radius[MedianUp],
         Radius[AcceptedNumIter-1-Central95]-Radius[MedianUp],
         Radius[MedianUp]-Radius[Central95]);



  delete hRadius;
  delete hOuterLow;
  delete hInner;
  delete hOuterUp;
  delete myPT;

  delete InRootFile;
  delete [] Radius;
  delete [] IterID;

  if(FileDefault) delete FileDefault;
  delete FileMedianLow;
  delete FileMedianUp;
  delete FileLowerLim;
  delete FileUpperLim;

  delete OutputFile;

  delete [] CutUsed;
}
