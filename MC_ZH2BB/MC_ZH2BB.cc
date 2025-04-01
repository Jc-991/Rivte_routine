// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  class MC_ZH2BB : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_ZH2BB()
      : Analysis("MC_ZH2BB")
    {    }

    //@}


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const double jetR = getOption("antikt", 0.4);
      const double Zptcut = getOption("Zptcut", 400) * GeV;

      FinalState fs;
      declare(fs, "FinalState");
      declare(FastJets(fs, JetAlg::ANTIKT, jetR), "AntiKT");
      declare(MissingMomentum(FinalState()), "MissingET");  // Declare Missing Transverse Energy (MET)

      /// Book histograms
      book(_h_jet_bb_Delta_eta ,"jet_bb_Delta_eta", 50, 0, 4);
      book(_h_jet_bb_Delta_phi ,"jet_bb_Delta_phi", 50, 0, 1);
      book(_h_jet_bb_Delta_pT ,"jet_bb_Delta_pT", 50, 0, 500);
      book(_h_jet_bb_Delta_R ,"jet_bb_Delta_R", 50, 0, 5);
      book(_h_jet_b_jet_eta ,"jet_b_jet_eta", 50, -4, 4);
      book(_h_jet_b_jet_multiplicity ,"jet_b_jet_multiplicity", 11, -0.5, 10.5);
      book(_h_jet_b_jet_phi ,"jet_b_jet_phi", 50, 0, 1);
      book(_h_jet_b_jet_pT ,"jet_b_jet_pT", 50, 0, 3000);
      book(_h_jet_H_eta_using_bb ,"jet_H_eta_using_bb", 50, -4, 4);
      book(_h_jet_H_mass_using_bb ,"jet_H_mass_using_bb", 50, 0, 200);
      book(_h_jet_H_phi_using_bb ,"jet_H_phi_using_bb", 50, 0, 1);
      book(_h_jet_H_pT_using_bb ,"jet_H_pT_using_bb", 50, 0, 3000);
      book(_h_jet_eta ,"jet_eta", 50, -4, 4);
      book(_h_jet_multiplicity ,"jet_multiplicity", 20, -0.5, 19.5);
      book(_h_jet_phi ,"jet_phi", 50, 0, 1);
      book(_h_jet_pT ,"jet_pT", 50, 0, 3000);

      book(_h_jet_H_eta_using_b ,"jet_H_eta_using_b", 50, -4, 4);
      book(_h_jet_H_mass_using_b ,"jet_H_mass_using_b", 50, 0, 200);
      book(_h_jet_H_phi_using_b ,"jet_H_phi_using_b", 50, 0, 1);
      book(_h_jet_H_pT_using_b ,"jet_H_pT_using_b", 50, 0, 3000);
      book(_h_jet_H_pT_using_b_boosted ,"jet_H_pT_using_b_boosted", 50, 0, 3000);

      book(_h_jet_H_eta_using_single_jet ,"jet_H_eta_using_single_jet", 50, -4, 4);
      book(_h_jet_H_mass_using_single_jet ,"jet_H_mass_using_single_jet", 50, 0, 200);
      book(_h_jet_H_phi_using_single_jet ,"jet_H_phi_using_single_jet", 50, 0, 1);
      book(_h_jet_H_pT_using_single_jet ,"jet_H_pT_using_single_jet", 50, 0, 3000);
      book(_h_jet_H_pT_using_single_jet_boosted ,"jet_H_pT_using_single_jet_boosted", 50, 0, 3000);

      book(_h_jet_H_eta_using_double_b_jets ,"jet_H_eta_using_double_b_jets", 50, -4, 4);
      book(_h_jet_H_mass_using_double_b_jets ,"jet_H_mass_using_double_b_jets", 50, 0, 200);
      book(_h_jet_H_phi_using_double_b_jets ,"jet_H_phi_using_double_b_jets", 50, 0, 1);
      book(_h_jet_H_pT_using_double_b_jets ,"jet_H_pT_using_double_b_jets", 50, 0, 3000);
      book(_h_jet_H_pT_using_double_b_jets_boosted ,"jet_H_pT_using_double_b_jets_boosted", 50, 0, 3000);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double JETPTCUT = 20*GeV;
      //const double Zptcut = getOption("Zptcut", 400) * GeV; 
      
      // Apply jet algorithm 
      const Jets jets = apply<FastJets>(event, "AntiKT").jetsByPt(Cuts::pT > JETPTCUT);
      _h_jet_multiplicity->fill(jets.size());


      // Missing transverse momentum, hence the Z momentum
      //const MissingMomentum& metProj = apply<MissingMomentum>(event, "MissingET");
      const FourMomentum met = apply<MissingMomentum>(event, "MissingET").missingMomentum();
      
      double Z_pT = met.pT();
      if (Z_pT < Zptcut) vetoEvent;

      // Identify the b-jets
      Jets bjets;
      Jets double_bjets;
      for (const Jet& jet : jets) {
        const double jetEta = jet.eta();
        const double jetPhi = jet.phi();
        const double jetPt = jet.pT(); //Change the Z pT
        const double jetMass = jet.mass();
        

        _h_jet_eta->fill(jetEta);
        _h_jet_phi->fill(jetPhi/2/M_PI);
        _h_jet_pT->fill(jetPt/GeV);

        if (jet.bTagged() && jet.pT() > JETPTCUT) {
          bjets.push_back(jet);
          _h_jet_b_jet_eta->fill(jetEta);
          _h_jet_b_jet_phi->fill(jetPhi/2/M_PI);
          _h_jet_b_jet_pT->fill(jetPt);
        }

        // Identidify the jets with 2 b-hadrons
        Particles bTags = jet.bTags();

        if (bTags.size() == 2) {
          double_bjets.push_back(jet);
          _h_jet_H_eta_using_double_b_jets->fill(jetEta);
          _h_jet_H_mass_using_double_b_jets->fill(jetMass);
          _h_jet_H_phi_using_double_b_jets->fill(jetPhi/2/M_PI);
          _h_jet_H_pT_using_double_b_jets->fill(jetPt);
          
          // Select double tagged b-jets in the mass window
          if (jetMass > 115*GeV && jetMass < 140*GeV) {
            _h_jet_H_pT_using_double_b_jets_boosted->fill(jetPt);
          }
        }
      }
      _h_jet_b_jet_multiplicity->fill(bjets.size());
      
      // rest of analysis requires at least 1 b jets
      if(bjets.empty()) vetoEvent;

      // Construct Higgs candidates from single b-jet
      if (bjets.size() == 1) {

	      const Jet& jetSingle = bjets[0];
	      const FourMomentum phiggsSingle = jetSingle.momentum();

        // Construct Higgs pT for events with mass in the mass window
        if (phiggsSingle.mass() > 115*GeV && phiggsSingle.mass() < 140*GeV) {

          _h_jet_H_pT_using_b_boosted->fill(phiggsSingle.pT());
        }
	      _h_jet_H_eta_using_b->fill(phiggsSingle.eta());
        _h_jet_H_mass_using_b->fill(phiggsSingle.mass());
        _h_jet_H_phi_using_b->fill(phiggsSingle.phi()/2/M_PI);
        _h_jet_H_pT_using_b->fill(phiggsSingle.pT());


        // Construct Higgs candidates from single jet (also b-jet)
        if (jets.size() == 1){
          const Jet& bjetSingle = bjets[0];
          const FourMomentum phiggsSingle_b = bjetSingle.momentum();

          _h_jet_H_eta_using_single_jet->fill(phiggsSingle_b.eta());
          _h_jet_H_mass_using_single_jet->fill(phiggsSingle_b.mass());
          _h_jet_H_phi_using_single_jet->fill(phiggsSingle_b.phi()/2/M_PI);
          _h_jet_H_pT_using_single_jet->fill(phiggsSingle_b.pT());

          if (phiggsSingle_b.mass() > 115*GeV && phiggsSingle_b.mass() < 140*GeV) {
            _h_jet_H_pT_using_single_jet_boosted->fill(phiggsSingle_b.pT());
          }
        }
      }
      
      if (bjets.size() == 2) {

        // Bubble sort the bjets based on their momentum
        for (size_t i = 0; i < bjets.size(); ++i) {
          for (size_t j = 0; j < bjets.size() - i - 1; ++j) {
	          if (bjets[j].pT() < bjets[j + 1].pT()) {
	            // Swap element
	            const Jet& temp = bjets[j];
	            bjets[j] = bjets[j + 1];
	            bjets[j + 1] = temp;
            }
          }
        }
        // Construct Higgs candidates from pairs of b-jets
        const Jet& bjet1 = bjets[0];
        const Jet& bjet2 = bjets[1];
        
        const double deltaEtaJJ = fabs(bjet1.eta() - bjet2.eta());
        const double deltaPhiJJ = deltaPhi(bjet1.momentum(), bjet2.momentum());
        const double deltaRJJ = deltaR(bjet1.momentum(), bjet2.momentum());
        const double deltaPtJJ = fabs(bjet1.pT() - bjet2.pT());
        _h_jet_bb_Delta_eta->fill(deltaEtaJJ);
        _h_jet_bb_Delta_phi->fill(deltaPhiJJ/M_PI);
        _h_jet_bb_Delta_pT->fill(deltaPtJJ);
        _h_jet_bb_Delta_R->fill(deltaRJJ);

        // 2 jets event with 2 b-jets
        if (jets.size() == 2) { 
          const FourMomentum phiggs2 = bjet1.momentum() + bjet2.momentum();
          _h_jet_H_eta_using_bb->fill(phiggs2.eta());
          _h_jet_H_mass_using_bb->fill(phiggs2.mass());
          _h_jet_H_phi_using_bb->fill(phiggs2.phi()/2/M_PI);
          _h_jet_H_pT_using_bb->fill(phiggs2.pT());
        }
 
      }
    } 

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_jet_bb_Delta_eta, crossSection() / sumOfWeights());
      scale(_h_jet_bb_Delta_phi, crossSection() / sumOfWeights());
      scale(_h_jet_bb_Delta_pT, crossSection() / sumOfWeights());
      scale(_h_jet_bb_Delta_R, crossSection() / sumOfWeights());
      scale(_h_jet_b_jet_eta, crossSection() / sumOfWeights());
      scale(_h_jet_b_jet_multiplicity, crossSection() / sumOfWeights());
      scale(_h_jet_b_jet_phi, crossSection() / sumOfWeights());
      scale(_h_jet_b_jet_pT, crossSection() / sumOfWeights());
      scale(_h_jet_H_eta_using_bb, crossSection() / sumOfWeights());
      scale(_h_jet_H_mass_using_bb, crossSection() / sumOfWeights());
      scale(_h_jet_H_phi_using_bb, crossSection() / sumOfWeights());
      scale(_h_jet_H_pT_using_bb, crossSection() / sumOfWeights());
      scale(_h_jet_eta, crossSection() / sumOfWeights());
      scale(_h_jet_multiplicity, crossSection() / sumOfWeights());
      scale(_h_jet_phi, crossSection() / sumOfWeights());
      scale(_h_jet_pT, crossSection() / sumOfWeights());
      scale(_h_jet_H_eta_using_b, crossSection() / sumOfWeights());
      scale(_h_jet_H_mass_using_b, crossSection() / sumOfWeights());
      scale(_h_jet_H_phi_using_b, crossSection() / sumOfWeights());
      scale(_h_jet_H_pT_using_b, crossSection() / sumOfWeights());
      scale(_h_jet_H_pT_using_b_boosted, crossSection() / sumOfWeights());
      scale(_h_jet_H_eta_using_single_jet, crossSection() / sumOfWeights());
      scale(_h_jet_H_mass_using_single_jet, crossSection() / sumOfWeights());
      scale(_h_jet_H_phi_using_single_jet, crossSection() / sumOfWeights());
      scale(_h_jet_H_pT_using_single_jet, crossSection() / sumOfWeights());
      scale(_h_jet_H_pT_using_single_jet_boosted, crossSection() / sumOfWeights());
      scale(_h_jet_H_eta_using_double_b_jets, crossSection() / sumOfWeights());
      scale(_h_jet_H_mass_using_double_b_jets, crossSection() / sumOfWeights());
      scale(_h_jet_H_phi_using_double_b_jets, crossSection() / sumOfWeights());
      scale(_h_jet_H_pT_using_double_b_jets, crossSection() / sumOfWeights());
      scale(_h_jet_H_pT_using_double_b_jets_boosted, crossSection() / sumOfWeights());
      scale(_h_jet_pT, crossSection() / sumOfWeights());


      // Count the number of 1 b-jet events
      int one_bjet_events = _h_jet_H_pT_using_b->numEntries();

      // Count the number of 1 b-jet events in the mass window
      int one_bjet_events_in = _h_jet_H_pT_using_b_boosted->numEntries();

      // Calculate efficiency
      double efficiency1 = static_cast<double>(one_bjet_events_in) / one_bjet_events;

      // Count the number of 1 jet, 1 b-jet events
      int one_jet_bjet_events = _h_jet_H_pT_using_single_jet->numEntries();

      // Count the number of 1 jet, 1 b-jet events in the mass window
      int one_jet_bjet_events_in = _h_jet_H_pT_using_single_jet_boosted->numEntries();

      // Calculate efficiency
      double efficiency11 = static_cast<double>(one_jet_bjet_events_in) / one_jet_bjet_events;

      // Count the number of double tagged b-jets
      int two_b_hadron_jets = _h_jet_H_pT_using_double_b_jets->numEntries();

      // Count the number of double tagged b jets in the mass window 
      int two_b_hadron_jets_in = _h_jet_H_pT_using_double_b_jets_boosted->numEntries();

      // Count the total number of events 
      int total_events = _h_jet_multiplicity->numEntries();

      // Calculate efficiency
      double efficiency2 = static_cast<double>(two_b_hadron_jets_in) / two_b_hadron_jets;
      double efficiency3 = static_cast<double>(two_b_hadron_jets_in) / total_events;
      double xs_pb = crossSection() / picobarn;  // Convert to pb

      // print the efficiencies
      cout << "Number of events with 1 b-jets: " << one_bjet_events << endl;
      cout << "Number of events with 1 b-jets in the mass window: "<< one_bjet_events_in << endl;
      cout << "The efficiency is " << efficiency1 << endl;

      cout << "Number of events with 1 jet, 1 b-jets: " << one_jet_bjet_events << endl;
      cout << "Number of events with 1 jet, 1 b-jets in the mass window: "<< one_jet_bjet_events_in << endl;
      cout << "The efficiency is " << efficiency11 << endl;

      cout << "Number of double b-tagged jets: " << two_b_hadron_jets << endl;
      cout << "Number of double b-tagged jets in the mass window: "<< two_b_hadron_jets_in << endl;
      cout << "The relative efficiency is " << efficiency2 << endl; 
      cout << "The absolute efficiency is " << efficiency3 << endl; 

      cout << "Cross section: " << xs_pb << " pb" << endl;    
      cout << "Total sum of event weights: " << sumOfWeights() << endl;
      cout << "Integrated Weight (MadGraph): " << xs_pb / sumOfWeights() << " pb" << endl;
      cout << "Applied Zptcut: " << Zptcut / GeV << " GeV" << endl;

    }


  private:

    /// @name Histograms
    //@{

    Histo1DPtr _h_jet_bb_Delta_eta, _h_jet_bb_Delta_phi, _h_jet_bb_Delta_pT, _h_jet_bb_Delta_R;
    Histo1DPtr _h_jet_b_jet_eta, _h_jet_b_jet_multiplicity, _h_jet_b_jet_phi, _h_jet_b_jet_pT;
    Histo1DPtr _h_jet_H_eta_using_bb, _h_jet_H_mass_using_bb, _h_jet_H_phi_using_bb, _h_jet_H_pT_using_bb;
    Histo1DPtr _h_jet_eta, _h_jet_multiplicity, _h_jet_phi, _h_jet_pT;
    Histo1DPtr _h_jet_H_eta_using_b, _h_jet_H_mass_using_b, _h_jet_H_phi_using_b, _h_jet_H_pT_using_b, _h_jet_H_pT_using_b_boosted;
    Histo1DPtr _h_jet_H_eta_using_single_jet, _h_jet_H_mass_using_single_jet, _h_jet_H_phi_using_single_jet, _h_jet_H_pT_using_single_jet, _h_jet_H_pT_using_single_jet_boosted;
    Histo1DPtr _h_jet_H_eta_using_double_b_jets, _h_jet_H_mass_using_double_b_jets, _h_jet_H_phi_using_double_b_jets, _h_jet_H_pT_using_double_b_jets, _h_jet_H_pT_using_double_b_jets_boosted;

    //@}

  };


  // This global object acts as a hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_ZH2BB);

}
