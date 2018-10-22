#include "ForBernie.h"
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
#include <iostream>
#include "stdlib.h"
void FitPPVariations(const unsigned& NumIter, const unsigned& NumJobs,
		const unsigned& JobID, int system, TString InputDir, TString OutputDir);

int main(int argc, char *argv[]) {
	FitPPVariations(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
			argv[5], argv[6]);
	return 0;
}

void FitPPVariations(const unsigned& NumIter, const unsigned& NumJobs,
		const unsigned& JobID, int system, TString InputDir,
		TString OutputDir) {

	const bool FAST_PLOT = true;
	TString HistName = "hCk_ReweightedMeV_";
	TString HistppName = HistName + "0";
	TString InputFilePrefix = "CFOutput_";
	const unsigned NumMomBins_pp = 50;
	const double kMin_pp = 0;
	const double kMax_pp = kMin_pp + 4 * NumMomBins_pp;  //(4 is the bin width)

	double FemtoRegion_pp[3][2];
	FemtoRegion_pp[0][0] = kMin_pp;
	FemtoRegion_pp[0][1] = 120;
	FemtoRegion_pp[1][0] = kMin_pp;
	FemtoRegion_pp[1][1] = 160;
	FemtoRegion_pp[2][0] = kMin_pp;
	FemtoRegion_pp[2][1] = 200;
	//!The baseline region (the same for all systems)
	double BlRegion[3][2];
	BlRegion[0][0] = 320;
	BlRegion[0][1] = 480;
	BlRegion[1][0] = 300;
	BlRegion[1][1] = 500;
	BlRegion[2][0] = 300;
	BlRegion[2][1] = 540;
	double PurityProton;
	double pp_f0;
	double pp_f1;
	PurityProton = 0.984266;  //pPb 5 TeV
	pp_f0 = 0.862814;
	pp_f1 = 0.09603;

	double ProtonPrim = pp_f0;
	double arrayPercLamProton[3] = { pp_f1 / (1. - pp_f0) * 0.8, pp_f1
			/ (1. - pp_f0), pp_f1 / (1. - pp_f0) * 1.2 }; //+/- 20%

	//following my lambda pars with the 3 possible modifications
	//for the proton:
	//0 = primary
	//1 = from Lambda
	//2 = other feeddown (flat)
	//3 = missidentified
	const unsigned NumChannels_p = 4;
	double** Purities_p = new double*[3];
	double** Fraction_p = new double*[3];
	for (unsigned uVar = 0; uVar < 3; uVar++) {
		Purities_p[uVar] = new double[NumChannels_p];
		Fraction_p[uVar] = new double[NumChannels_p];

		Purities_p[uVar][0] = PurityProton;
		Purities_p[uVar][1] = PurityProton;
		Purities_p[uVar][2] = PurityProton;
		Purities_p[uVar][3] = 1. - PurityProton;

		Fraction_p[uVar][0] = ProtonPrim;
		Fraction_p[uVar][1] = (1. - ProtonPrim) * (arrayPercLamProton[uVar]);
		Fraction_p[uVar][2] = (1. - ProtonPrim)
				* (1. - arrayPercLamProton[uVar]);
		Fraction_p[uVar][3] = 1.;
	}
	//starting value, do not worry about it too much
	const double GaussSourceSize = 1.2;

	//ADVANCED***
	//!DO NOT TOUCH UNLESS YOU CHANGE THE ABOVE FILES, IN WHICH CASE ASK DIMI FOR HELP
	//1/FractionOfBins th number of original bins of the correction matrix are taken into account
	//originally: 1000 MeV => 1/2 should be good most of the time
	const int Fraction_Res = 2;
	const int Fraction_Sig = 1;
	const double UnitConv_Res = 1;
	const double UnitConv_Sig = 1;
	TString CalibBaseDir = "";
	std::cout << "SYSTEM: " << system << std::endl;
	if (system == 0) { // pPb MB
		CalibBaseDir +=
				"/home/hohlweger/cernbox/SystematicsAndCalib/pPbRun2_MB";
	} else if (system == 1) { // pp MB
		CalibBaseDir += "/home/hohlweger/cernbox/SystematicsAndCalib/ppRun2_MB";
	} else if (system == 2) { // pp HM
		CalibBaseDir += "/home/hohlweger/cernbox/SystematicsAndCalib/ppRun2_HM";
	} else if (system == 11) { //for use on the cluster
		CalibBaseDir +=
				"/home/hohlweger/cernbox/SystematicsAndCalib/pPbRun2_MB";
	}

	const TString ResMatrixFileName = Form("%s/run2_decay_matrices_old.root",
			CalibBaseDir.Data());
	const TString SigmaMatrixFileName = Form("%s/Sample3_MeV_compact.root",
			CalibBaseDir.Data());
	TH2F* hRes_pp_pL;
	TFile* FileRes = new TFile(ResMatrixFileName, "read");
	FileRes->cd();
	hRes_pp_pL = (TH2F*) FileRes->Get("hRes_pp_pL");
	TH2F* hRes_pp_pL_MeV = new TH2F("hRes_pp_pL_MeV", "hRes_pp_pL_MeV",
			hRes_pp_pL->GetNbinsX() / Fraction_Res,
			hRes_pp_pL->GetXaxis()->GetBinLowEdge(1) * UnitConv_Res,
			hRes_pp_pL->GetXaxis()->GetBinUpEdge(
					hRes_pp_pL->GetNbinsX() / Fraction_Res) * UnitConv_Res,
			hRes_pp_pL->GetNbinsY() / Fraction_Res,
			hRes_pp_pL->GetYaxis()->GetBinLowEdge(1) * UnitConv_Res,
			hRes_pp_pL->GetXaxis()->GetBinUpEdge(
					hRes_pp_pL->GetNbinsY() / Fraction_Res) * UnitConv_Res);

	for (int iBinX = 1; iBinX <= hRes_pp_pL->GetNbinsX() / Fraction_Res;
			iBinX++) {
		for (int iBinY = 1; iBinY <= hRes_pp_pL->GetNbinsY() / Fraction_Res;
				iBinY++) {
			hRes_pp_pL_MeV->SetBinContent(iBinX, iBinY,
					hRes_pp_pL->GetBinContent(iBinX, iBinY));
		}
	}
	TH2F* hSigma_pp;
	TFile* FileSigma = new TFile(SigmaMatrixFileName, "read");
	FileSigma->cd();
	hSigma_pp = (TH2F*) FileSigma->Get("hSigmaMeV_Proton_Proton");

	TH2F* hSigma_pp_MeV = new TH2F("hSigma_pp_MeV", "hSigma_pp_MeV",
			hSigma_pp->GetNbinsX() / Fraction_Sig,
			hSigma_pp->GetXaxis()->GetBinLowEdge(1) * UnitConv_Sig,
			hSigma_pp->GetXaxis()->GetBinUpEdge(
					hSigma_pp->GetNbinsX() / Fraction_Sig) * UnitConv_Sig,
			hSigma_pp->GetNbinsY() / Fraction_Sig,
			hSigma_pp->GetYaxis()->GetBinLowEdge(1) * UnitConv_Sig,
			hSigma_pp->GetXaxis()->GetBinUpEdge(
					hSigma_pp->GetNbinsY() / Fraction_Sig) * UnitConv_Sig);
	for (int iBinX = 1; iBinX <= hSigma_pp->GetNbinsX() / Fraction_Sig;
			iBinX++) {
		for (int iBinY = 1; iBinY <= hSigma_pp->GetNbinsY() / Fraction_Sig;
				iBinY++) {
			hSigma_pp_MeV->SetBinContent(iBinX, iBinY,
					hSigma_pp->GetBinContent(iBinX, iBinY));
		}
	}
	int vFemReg_pp;  //which femto region we use for pp (1 = default)
	int vBlReg;  //which baseline region to use (1 = default)
	int vFrac_pp_pL; //fraction of protons coming from Lambda variation (1 = default)

	//each JOB produces a separate output file
	TFile* OutFile = new TFile(
			TString::Format("%sOutFile_%s_Iter%u_JOBS%u_ID%u.root",
					OutputDir.Data(), "CutVarAdd", NumIter, NumJobs, JobID),
			"recreate");
	//you save a lot of stuff in an NTuple
	TNtuple* ntResult =
			new TNtuple("ntResult", "ntResult",
					"IterID:vFemReg_pp:vBlReg:vFrac_pp_pL:Radius_pp:RadiusErr_pp:p_a:p_b:p_Cl:Chi2Ndf:pval");
	unsigned WhichPart = JobID + 1;
	unsigned NumJobsPart = NumIter;

	unsigned iSplitInto = NumJobs;
	unsigned FirstIter = 0;
	unsigned LastIter = 0;
	Printf("JobID = %i", JobID);
	while (iSplitInto > 0 && WhichPart > 0) {
		FirstIter = LastIter + bool(LastIter);
		LastIter += NumJobsPart / iSplitInto - 1 + bool(LastIter);
		NumJobsPart = NumIter - LastIter - 1;
		iSplitInto--;
		WhichPart--;
	}
	if (NumIter == NumJobs) {
		FirstIter = JobID;
		LastIter = JobID;
	}

	Printf(" FirstIter = %i", FirstIter);
	Printf(" LastIter = %i", LastIter);
	Printf("  Nruns = %i", LastIter - FirstIter + 1);
	//***
	TRandom3 rangen(1 + JobID);
	//the 0 iter is always the default!
	for (unsigned uIter = FirstIter; uIter <= LastIter; uIter++) {
		vFemReg_pp = rangen.Integer(3);
		vBlReg = rangen.Integer(3);
		vFrac_pp_pL = rangen.Integer(3);

		//The defaults
		if (uIter == 0) {
			vFemReg_pp = 1;
			vBlReg = 1;
			vFrac_pp_pL = 1;
		}

		//true => do the BL separately (as an RUN1), false => fit femto and BL region together (RUN2)
		bool SEPARATE_BL = false;
		//if true renormalization is NOT allowed
		bool FIX_CL = false;

		//ADVANCED***
		//#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
		double PotPars1S0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0,
				0 };
		double PotPars3P0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1,
				0 };
		double PotPars3P1[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1,
				1 };
		double PotPars3P2[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1,
				2 };

		const double Weight1S0 = 3. / 12.;
		const double Weight3P0 = 1. / 12.;
		const double Weight3P1 = 3. / 12.;
		const double Weight3P2 = 5. / 12.;
		const double Mass_p = 938.272;
		const double Mass_L = 1115.683;

		CATS AB_pp;
		double Pars_pp[6] = { 0, 0, 0, GaussSourceSize * 1.2, GaussSourceSize
				/ 1.2, 0.5 };
		AB_pp.SetAnaSource(GaussSource, Pars_pp);
		AB_pp.SetUseAnalyticSource(true);
		AB_pp.SetThetaDependentSource(false);

		AB_pp.SetExcludeFailedBins(false);
		AB_pp.SetMomBins(NumMomBins_pp, kMin_pp, kMax_pp);

		AB_pp.SetNumChannels(4);
		AB_pp.SetNumPW(0, 2);
		AB_pp.SetNumPW(1, 2);
		AB_pp.SetNumPW(2, 2);
		AB_pp.SetNumPW(3, 2);
		AB_pp.SetSpin(0, 0);
		AB_pp.SetSpin(1, 1);
		AB_pp.SetSpin(2, 1);
		AB_pp.SetSpin(3, 1);
		AB_pp.SetChannelWeight(0, Weight1S0);
		AB_pp.SetChannelWeight(1, Weight3P0);
		AB_pp.SetChannelWeight(2, Weight3P1);
		AB_pp.SetChannelWeight(3, Weight3P2);

		AB_pp.SetQ1Q2(1);
		AB_pp.SetPdgId(2212, 2212);
		AB_pp.SetRedMass(0.5 * Mass_p);

		AB_pp.SetShortRangePotential(0, 0, fDlmPot, PotPars1S0);
		AB_pp.SetShortRangePotential(1, 1, fDlmPot, PotPars3P0);
		AB_pp.SetShortRangePotential(2, 1, fDlmPot, PotPars3P1);
		AB_pp.SetShortRangePotential(3, 1, fDlmPot, PotPars3P2);

		AB_pp.KillTheCat();

		CATS AB_pL;
		double Pars_pL[6] = { 0, 0, 0, GaussSourceSize * 1.2, GaussSourceSize
				/ 1.2, 0.5 };
		AB_pL.SetAnaSource(GaussSource, Pars_pL);
		AB_pL.SetUseAnalyticSource(true);
		AB_pL.SetThetaDependentSource(false);
		AB_pL.SetMomBins(NumMomBins_pp, kMin_pp, kMax_pp);
		AB_pL.SetExcludeFailedBins(false);

		AB_pL.SetNumChannels(2);
		AB_pL.SetNumPW(0, 1);
		AB_pL.SetNumPW(1, 1);
		AB_pL.SetSpin(0, 0);
		AB_pL.SetSpin(1, 1);
		AB_pL.SetChannelWeight(0, 0.25);
		AB_pL.SetChannelWeight(1, 0.75);

		AB_pL.SetQ1Q2(0);
		AB_pL.SetPdgId(2212, 3122);
		AB_pL.SetRedMass((Mass_p * Mass_L) / (Mass_p + Mass_L));

		//!TEMPORARY: WE USE USMANI, ASK DIMI ABOUT FURTHER CHANGES (LIKE USE THE NLO)
		//#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
		double pLamPotPars1S0[10] = { 0, 0, pL_UsmaniOli, 0, 0, 0, 0, 0, 0, 0 };
		double pLamPotPars3S1[10] = { 0, 0, pL_UsmaniOli, 0, 0, 0, 0, 1, 0, 1 };
		AB_pL.SetShortRangePotential(0, 0, fDlmPot, pLamPotPars1S0);
		AB_pL.SetShortRangePotential(1, 0, fDlmPot, pLamPotPars3S1);

		AB_pL.KillTheCat();

		std::cout << "Reading Data \n";
		//! DATA FILE
		//    TString OliFileName_pp =
		TString OliFileName_pp = TString::Format("%s%spp.root", InputDir.Data(),
				InputFilePrefix.Data());
		TFile* OliFile_pp =
				OliFileName_pp != "" ? new TFile(OliFileName_pp, "read") : NULL;
		TH1F* OliHisto_pp =
				OliFile_pp ? (TH1F*) OliFile_pp->Get(HistppName.Data()) : NULL;
		if (OliHisto_pp)
			std::cout << OliHisto_pp->GetName() << std::endl;
		else {
			std::cout << OliFileName_pp.Data() << "--" << HistppName.Data()
					<< " Missing" << std::endl;
		}
		//!CHANGE PATH HERE
		TString SystErrFileName_pp = TString::Format("%s/C2totalsysPP.root",
				CalibBaseDir.Data());
		TFile* SystErrFile_pp =
				SystErrFileName_pp != "" ?
						new TFile(SystErrFileName_pp, "read") : NULL;
		TH1F* outputParamPP = (TH1F*) SystErrFile_pp->Get("SysParamPP");
		std::cout << "PP" << std::endl;
		std::cout << outputParamPP->GetBinContent(1) << std::endl;
		std::cout << outputParamPP->GetBinContent(2) << std::endl;
		std::cout << outputParamPP->GetBinContent(3) << std::endl;
		TF1 *RelSystPP = new TF1("sysPP", "pol2", 0, 3);
		RelSystPP->SetParameter(0, outputParamPP->GetBinContent(1));
		RelSystPP->SetParameter(1, outputParamPP->GetBinContent(2));
		RelSystPP->SetParameter(2, outputParamPP->GetBinContent(3));

		int NumSEB_pp = RelSystPP == NULL ? 0 : OliHisto_pp->FindBin(500);

		for (int iBin = 0; iBin < NumSEB_pp; iBin++) {
			const float x = OliHisto_pp->GetBinCenter(iBin + 1);
			const float y = OliHisto_pp->GetBinContent(iBin + 1);
			OliHisto_pp->SetBinError(iBin + 1,
					sqrt(
							pow(OliHisto_pp->GetBinError(iBin + 1), 2.)
									+ pow(y * RelSystPP->Eval(x / 1000.), 2.)));
		}
		const unsigned NumSourcePars = 1;
		DLM_Ck* Ck_pp = new DLM_Ck(NumSourcePars, 0, AB_pp);
		DLM_Ck* Ck_pL = new DLM_Ck(NumSourcePars, 0, AB_pL);
		Ck_pp->Update();
		Ck_pL->Update();
		DLM_CkDecomposition CkDec_pp("pp", 2, *Ck_pp, hSigma_pp_MeV);
		DLM_CkDecomposition CkDec_pL("pLambda", 0, *Ck_pL, NULL);
		double lam_pp = Purities_p[vFrac_pp_pL][0] * Fraction_p[vFrac_pp_pL][0]
				* Purities_p[vFrac_pp_pL][0] * Fraction_p[vFrac_pp_pL][0];
		double lam_pp_pL = Purities_p[vFrac_pp_pL][0]
				* Fraction_p[vFrac_pp_pL][0] * Purities_p[vFrac_pp_pL][1]
				* Fraction_p[vFrac_pp_pL][1] * 2;
		double lam_pp_fake = Purities_p[vFrac_pp_pL][3]
				* Purities_p[vFrac_pp_pL][0]
				+ Purities_p[vFrac_pp_pL][0] * Purities_p[vFrac_pp_pL][3]
				+ Purities_p[vFrac_pp_pL][3] * Purities_p[vFrac_pp_pL][3];
		printf("lam_pp = %.3f\n", lam_pp);
//		CkDec_pp.AddContribution(0, lam_pp_pL, DLM_CkDecomposition::cFeedDown,
//				&CkDec_pL, hRes_pp_pL_MeV);
//		CkDec_pp.AddContribution(0, 1. - lam_pp - lam_pp_pL - lam_pp_fake,
//				DLM_CkDecomposition::cFeedDown);
		CkDec_pp.AddContribution(0, 1. - lam_pp - lam_pp_fake,
				DLM_CkDecomposition::cFeedDown);
		CkDec_pp.AddContribution(1, lam_pp_fake, DLM_CkDecomposition::cFake); //0.02
		DLM_Fitter1* fitter = new DLM_Fitter1(1);
		fitter->SetOutputDir(OutputDir.Data());
		std::cout << "pp" << std::endl;
		fitter->SetSystem(0, *OliHisto_pp, 1, CkDec_pp,
				FemtoRegion_pp[vFemReg_pp][0], FemtoRegion_pp[vFemReg_pp][1],
				BlRegion[vBlReg][0], BlRegion[vBlReg][1]);
		fitter->SetSeparateBL(0, false);
//		fitter->AddSameSource("pLambda", "pp", 1);
//		fitter->SetParameter("pp", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
//		fitter->SetParameter("pp", DLM_Fitter1::p_b, 1e-4, 0, 2e-3);
//		fitter->FixParameter("pp", DLM_Fitter1::p_c, 0);
//		fitter->SetParameter("pp", DLM_Fitter1::p_Cl, -0.9, -1.2, -0.8);
//		fitter->SetParameter("pp", DLM_Fitter1::p_sor0, 1.2, 0.8, 1.8);

		fitter->FixParameter("pp", DLM_Fitter1::p_a, 0.990232);
		fitter->FixParameter("pp", DLM_Fitter1::p_b, 0.000190822);
		fitter->FixParameter("pp", DLM_Fitter1::p_c, 0);
		fitter->FixParameter("pp", DLM_Fitter1::p_Cl, -0.947533);
		fitter->FixParameter("pp", DLM_Fitter1::p_sor0, 1.359);

		std::cout << "updating \n";
		CkDec_pp.Update();
		//CkDec_pL.Update();
		std::cout << "Fitting \n";
		fitter->GoBabyGo();
		std::cout << "Done Fitting \n";
		float ntBuffer[11];
		ntBuffer[0] = uIter;
		ntBuffer[1] = vFemReg_pp;
		ntBuffer[2] = vBlReg;
		ntBuffer[3] = vFrac_pp_pL;
		ntBuffer[4] = fitter->GetParameter("pp", DLM_Fitter1::p_sor0);
		ntBuffer[5] = fitter->GetParError("pp", DLM_Fitter1::p_sor0);
		ntBuffer[6] = fitter->GetParameter("pp", DLM_Fitter1::p_a);
		ntBuffer[7] = fitter->GetParameter("pp", DLM_Fitter1::p_b);
		ntBuffer[8] = fitter->GetParameter("pp", DLM_Fitter1::p_Cl);
		ntBuffer[9] = fitter->GetChi2Ndf();
		ntBuffer[10] = fitter->GetPval();
		ntResult->Fill(ntBuffer);
		std::cout << "fast plot \n";

		if (FAST_PLOT) {
			TGraph FitResult_pp;
			FitResult_pp.SetName(TString::Format("FitResult_pp_%u", uIter));
			fitter->GetFitGraph(0, FitResult_pp);

			double Chi2_pp = 0;
			unsigned EffNumBins_pp = 0;
			for (unsigned uBin = 0; uBin < NumMomBins_pp; uBin++) {

				double mom = AB_pp.GetMomentum(uBin);
				double dataY;
				double dataErr;
				double theoryX;
				double theoryY;

				if (mom > FemtoRegion_pp[vFemReg_pp][1])
					continue;

				FitResult_pp.GetPoint(uBin, theoryX, theoryY);
				if (mom != theoryX) {
					std::cout << mom << '\t' << theoryX << std::endl;
					printf("  PROBLEM pp!\n");
				}
				dataY = OliHisto_pp->GetBinContent(uBin + 1);
				dataErr = OliHisto_pp->GetBinError(uBin + 1);

				Chi2_pp += (dataY - theoryY) * (dataY - theoryY)
						/ (dataErr * dataErr);
				EffNumBins_pp++;
			}

			TPaveText* info1 = new TPaveText(0.45, 0.65, 0.9, 0.95, "blNDC"); //lbrt
			info1->SetName("info1");
			info1->SetBorderSize(1);
			info1->SetTextSize(0.04);
			info1->SetFillColor(kWhite);
			info1->SetTextFont(22);
			TString SOURCE_NAME = "Gauss";

			info1->AddText(
					TString::Format("R(%s)=%.3f#pm%.3f", SOURCE_NAME.Data(),
							ntBuffer[4], ntBuffer[5]));
			info1->AddText(
					TString::Format("C_{l}=%.3f#pm%.3f",
							fitter->GetParameter("pp", DLM_Fitter1::p_Cl),
							fitter->GetParError("pp", DLM_Fitter1::p_Cl)));
			info1->AddText(
					TString::Format("Local #chi^{2}_{ndf}=%.2f, p_{val}=%.3f",
							Chi2_pp / double(EffNumBins_pp), ntBuffer[10]));
			double Yoffset = 1.2;
			TH1F* hAxis_pp = new TH1F("hAxis_pp", "hAxis_pp", 600, 0, 600);
			hAxis_pp->SetStats(false);
			hAxis_pp->SetTitle("");
			hAxis_pp->GetXaxis()->SetLabelSize(0.065);
			hAxis_pp->GetXaxis()->CenterTitle();
			hAxis_pp->GetXaxis()->SetTitleOffset(1.35);
			hAxis_pp->GetXaxis()->SetLabelOffset(0.02);
			hAxis_pp->GetXaxis()->SetTitleSize(0.075);
			hAxis_pp->GetYaxis()->SetLabelSize(0.065);
			hAxis_pp->GetYaxis()->CenterTitle();
			hAxis_pp->GetYaxis()->SetTitleOffset(Yoffset);
			hAxis_pp->GetYaxis()->SetTitleSize(0.075);
			hAxis_pp->GetXaxis()->SetRangeUser(0, kMax_pp);
			hAxis_pp->GetYaxis()->SetRangeUser(0.5, 3);    //pPb

			TCanvas* cfast = new TCanvas("cfast", "cfast", 1);
			cfast->cd(0);
			cfast->SetCanvasSize(1920, 1280);
			cfast->SetMargin(0.15, 0.05, 0.2, 0.05);    //lrbt

			OliHisto_pp->SetStats(false);
			OliHisto_pp->SetTitle("pp");
			OliHisto_pp->SetLineWidth(2);
			OliHisto_pp->SetLineColor(kBlack);
			FitResult_pp.SetLineWidth(2);
			FitResult_pp.SetLineColor(kRed);
			FitResult_pp.SetMarkerStyle(24);
			FitResult_pp.SetMarkerColor(kRed);
			FitResult_pp.SetMarkerSize(1);

			hAxis_pp->Draw("axis");
			OliHisto_pp->Draw("same");
			FitResult_pp.Draw("CP,same");
			info1->Draw("same");
			cfast->SaveAs(Form("%scfast_pponly.png", OutputDir.Data()));
		}
		OutFile->cd();
		ntResult->Write();
	}
	OutFile->Close();
}

