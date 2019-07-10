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

int main(int, char**) {
  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

  NetworkOptions opts;
  opts.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
  opts.MassDeviationThreshold = 1.0E-10;
  opts.IsSelfHeating = true;
  opts.EnableScreening = true;

  SkyNetScreening screen(nuclib);
  HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");

  REACLIBReactionLibrary strongReactionLibrary(SkyNetRoot + "/data/reaclib",
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
  REACLIBReactionLibrary weakReactionLibrary(SkyNetRoot + "/data/reaclib",
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

  double T0 = 10.0;
  double Ye = 0.17;
  double s = 30.0;
  double tau = 25.0;

  // run NSE with the temperature and entropy to find the initial density
  auto nseResult = nse.CalcFromTemperatureAndEntropy(T0, s, Ye);

  auto densityProfile = ExpTMinus3(nseResult.Rho(), tau / 1000.0);

  auto output = net.EvolveSelfHeatingWithInitialTemperature(nseResult.Y(), 0.0,
      1.0E9, T0, &densityProfile, "SkyNet_r-process");

  std::vector<double> finalYVsA = output.FinalYVsA();

  FILE * f = fopen("final_y_r-process", "w");
  for (unsigned int A = 0; A < finalYVsA.size(); ++A)
    fprintf(f, "%6i  %30.20E\n", A, finalYVsA[A]);
}

