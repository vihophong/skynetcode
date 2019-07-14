/// \file r-process.cpp
/// \author jlippuner
/// \since Sep 3, 2014
///
/// \brief
///
///

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/SkyNetScreening.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Network/NSE.hpp"
#include "DensityProfiles/ExpTMinus3.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TF1.h"

using namespace std;
void plota(char* infilelist,char* innamelist=0)
{
    TCanvas *c1=new TCanvas("c1","aaa",900,700);
    c1->SetLogy();
    ifstream fdat;
    fdat.open("/home/phong/projects/SkyNet/runskynet/abundance_ref/r-abundance.txt");
    Double_t nmass[141];
    Double_t stan[141];
    Double_t min[141];
    Double_t max[141];
    Double_t temp2,fmin,fmax;
    for (Int_t i=0;i<141;i++){
      fdat>>temp2>>nmass[i]>>stan[i]>>fmin>>fmax;
      //cout<<i<<"\t"<<nmass[i]<<"\t"<<stan[i]<<endl;
      min[i]=stan[i]-fmin;
      max[i]=fmax-stan[i];
    }
    Double_t stan_153Eu=stan[84];

    TGraphAsymmErrors* gr1=new TGraphAsymmErrors(141,nmass,stan,0,0,min,max);
    gr1->SetMarkerStyle(20);
    gr1->SetFillColor(0);
    gr1->GetYaxis()->SetRangeUser(1e-3,1e3);
    gr1->Draw("AP");

    std::ifstream ifs(infilelist);
    string filelist[1000];
    Int_t nfiles=0;
    while (!ifs.eof()){
        ifs>>filelist[nfiles];
        cout<<filelist[nfiles]<<endl;
        nfiles++;
    }
    nfiles=nfiles-1;
    cout<<"There are "<<nfiles<<" files in total!"<<endl;

    string namelist[1000];
    Int_t nlines=0;
    if (innamelist){
        std::ifstream ifsname(innamelist);
        while (!ifsname.eof()){
            std::getline(ifsname,namelist[nlines]);
            nlines++;
        }
        nlines=nlines-1;
        cout<<"There are "<<nlines<<" lines in total!"<<endl;
    }

    TGraph *grs[nfiles];
    Int_t colorcode=2;
    for (Int_t i=0;i<nfiles;i++){
        if (colorcode==4) colorcode++;
        grs[i]=new TGraph(filelist[i].c_str(),"%lg %lg");
        grs[i]->SetLineColor(colorcode);
        grs[i]->SetFillColor(0);
        Double_t normfactor=1;
        for (Int_t j=0;j<grs[i]->GetN();j++){
            if (grs[i]->GetX()[j]==153) normfactor=grs[i]->GetY()[j];
        }
        normfactor=stan_153Eu/normfactor;

        for (Int_t j=0;j<grs[i]->GetN();j++){
            grs[i]->GetY()[j]=normfactor*grs[i]->GetY()[j];
        }
        grs[i]->Draw("SAME");
        colorcode++;
    }
    grs[0]->Draw("SAME");
    //grs[1]->Draw("SAME");

    TLegend* leg = new TLegend(0.1,0.65,0.38,0.9);
    leg->AddEntry(gr1,"SS r-abudances");
    for (Int_t i=0;i<nfiles;i++){
        if (innamelist) leg->AddEntry(grs[i], namelist[i].c_str());
        else leg->AddEntry(grs[i], filelist[i].c_str());
    }
    leg->SetBorderSize(0);
    leg->Draw();
}


int main(int argc, char**argv) {
    //! read inputs

    // default input: BH-BH merger
    double T0 = 6.0;
    double Ye = 0.01;
    double s = 10.0;
    double tau = 7.1;

    char outputname[500];
    sprintf(outputname,"SkyNet_r-process");

    char reaclibdatabasechar[500];
    std::string reaclibdatabase=SkyNetRoot + "data/reaclib";

    std::cout<<"\n\n ***************** \n\nUsaged: ./r-process initial_condition path_to_output_file path_to_reaclib \n\n ***************** \n\n"<<std::endl;

    if (argc>1&&argc<5){
        std::ifstream ifs(argv[1]);
        char temp[500];
        ifs>>temp>>T0;
        ifs>>temp>>Ye;
        ifs>>temp>>s;
        ifs>>temp>>tau;
        if (argc>2){
            sprintf(outputname,"%s",argv[2]);
        }
        if (argc==4){
            sprintf(reaclibdatabasechar,"%s",argv[3]);
            reaclibdatabase=std::string(reaclibdatabasechar);
        }
    }

    std::cout<<"R-process condition: T0 = "<<T0<<"\tYe = "<<Ye<<"\ts = "<<s<<"\ttau = "<<tau<<std::endl;
    std::cout<<"Path to output file: "<<outputname<<std::endl;
    std::cout<<"Path to reactlib: "<<reaclibdatabase<<std::endl<<std::endl;


    auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
                SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

    NetworkOptions opts;
    opts.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
    opts.MassDeviationThreshold = 1.0E-10;
    opts.IsSelfHeating = true;
    opts.EnableScreening = true;

    SkyNetScreening screen(nuclib);
    HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");

    REACLIBReactionLibrary strongReactionLibrary(reaclibdatabase,
                                                 ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
                                                 "Strong reactions", nuclib, opts, true);
    REACLIBReactionLibrary symmetricFission(SkyNetRoot +
                                            "/data/netsu_panov_symmetric_0neut", ReactionType::Strong, false,
                                            LeptonMode::TreatAllAsDecayExceptLabelEC,
                                            "Symmetric neutron induced fission with 0 free neutrons", nuclib, opts,
                                            false);
    REACLIBReactionLibrary spontaneousFission(SkyNetRoot +
                                              "/data/netsu_sfis_Roberts2010rates", ReactionType::Strong, false,
                                              LeptonMode::TreatAllAsDecayExceptLabelEC, "Spontaneous fission", nuclib,
                                              opts, false);

    // use only REACLIB weak rates
    REACLIBReactionLibrary weakReactionLibrary(reaclibdatabase,
                                               ReactionType::Weak, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
                                               "Weak reactions", nuclib, opts, true);

    // or use the following code to use FFN rates and weak REACLIB rates
    // pre-computed with the <SkyNetRoot>/examples/precompute_reaction_libs.py
    // script
    //  auto ffnMesaReactionLibrary = FFNReactionLibrary::ReadFromDisk(
    //      "ffnMesa_with_neutrino", opts);
    //  auto ffnReactionLibrary = FFNReactionLibrary::ReadFromDisk(
    //      "ffn_with_neutrino_ffnMesa", opts);
    //  auto weakReactionLibrary = REACLIBReactionLibrary::ReadFromDisk(
    //      "weak_REACLIB_with_neutrino_ffnMesa_ffn", opts);

    // add neutrino reactions, if desired (will need a neutrino distribution
    // function set with net.LoadNeutrinoHistory(...))
    //  NeutrinoReactionLibrary neutrinoLibrary(SkyNetRoot
    //     + "/data/neutrino_reactions.dat", "Neutrino interactions", nuclib, opts,
    //     1.e-2, false, true);

    ReactionLibs reactionLibraries { &strongReactionLibrary, &symmetricFission,
                &spontaneousFission, &weakReactionLibrary };

    //  ReactionLibs reactionLibraries { &strongReactionLibrary, &symmetricFission,
    //    &spontaneousFission, &weakReactionLibrary,  &neutrinoLibrary,
    //    &ffnMesaReactionLibrary, &ffnReactionLibrary };

    ReactionNetwork net(nuclib, reactionLibraries, &helm, &screen, opts);
    NSE nse(net.GetNuclideLibrary(), &helm, &screen);

    // run NSE with the temperature and entropy to find the initial density
    auto nseResult = nse.CalcFromTemperatureAndEntropy(T0, s, Ye);

    auto densityProfile = ExpTMinus3(nseResult.Rho(), tau / 1000.0);

    auto output = net.EvolveSelfHeatingWithInitialTemperature(nseResult.Y(), 0.0,
                                                              1.0E9, T0, &densityProfile, outputname);

    std::vector<double> finalYVsA = output.FinalYVsA();
    std::vector<double> finalY = output.FinalY();

    char outputname_txt[500];
    sprintf(outputname_txt,"%s.txt",outputname);
    char outputname_all_txt[500];
    sprintf(outputname_all_txt,"%s_all.txt",outputname);

    FILE * f = fopen(outputname_txt, "w");
    for (unsigned int A = 0; A < finalYVsA.size(); ++A)
        fprintf(f, "%6i  %30.20E\n", A, finalYVsA[A]);
    FILE * f2 = fopen(outputname_txt, "w");
    for (unsigned int A = 0; A < finalY.size(); ++A)
        fprintf(f2, "%6i  %30.20E\n", A, finalY[A]);

}



