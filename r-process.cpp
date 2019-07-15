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


void makeyvsz(char* infile, char* outfile)
{
    string x;
    ifstream fin(infile);
    string breakstring="# [6] = initial abundance (of Z)";
    while (!fin.eof()){
        std::getline(fin,x);
        if(!x.compare(breakstring)) break;
    }
    Double_t temp;
    Double_t z;
    Double_t y;

    ofstream fout(outfile);

    while (!fin.eof()){
        fin>>temp>>temp>>temp>>z>>y>>temp;
        if (z==0&&y==0) break;
        fout<<z<<"\t"<<y<<endl;
    }
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


    //! Outputs
    std::vector<double> finalYVsA = output.FinalYVsA();

    char outputname_yvsa_txt[500];
    sprintf(outputname_yvsa_txt,"%s_YvsA.txt",outputname);

    FILE * f = fopen(outputname_yvsa_txt, "w");
    for (unsigned int A = 0; A < finalYVsA.size(); ++A)
        fprintf(f, "%6i  %30.20E\n", A, finalYVsA[A]);


    char outputname_h5[500];
    sprintf(outputname_h5,"%s.h5",outputname);

    output.MakeDatFile(outputname_h5);

    char inputdatfile[500];
    sprintf(inputdatfile,"%s.dat",outputname_h5);

    char outputname_yvsz_txt[500];
    sprintf(outputname_yvsz_txt,"%s_YvsZ.txt",outputname);
    makeyvsz(inputdatfile,outputname_yvsz_txt);
}



