#ifndef CALOVALID_CALOVALID_H
#define CALOVALID_CALOVALID_H

#include <fun4all/SubsysReco.h>
#include <bitset>
#include <string>
#include <vector>
#include <iostream>  // for std::cout, std::endl
#include <TH1F.h>    // for TH1F forward usage
#include <calotrigger/TriggerAnalyzer.h>  // or the correct path in your setup
#include <TH3F.h>


// Forward declarations
class PHCompositeNode;
class TH1;
class TH2;
class TH2F;
class TH3F;
class TProfile2D;
class TProfile;

class CaloValid : public SubsysReco
{
 public:
  //! constructor
  CaloValid(const std::string& name = "CaloValid");
  //! destructor
  ~CaloValid() override;

  //! full initialization
  int Init(PHCompositeNode*) override;

  //! event processing method
  int process_event(PHCompositeNode*) override;

  //! end of run method
  int End(PHCompositeNode*) override;

  // -----------------------------------------------------------------------
  // Optional methods. You can keep them empty or define them if needed.
  // -----------------------------------------------------------------------
  int process_g4hits(PHCompositeNode*);
  int process_g4cells(PHCompositeNode*);
  int process_towers(PHCompositeNode*);
  int process_clusters(PHCompositeNode*);

  // -----------------------------------------------------------------------
  // Simple parameter toggles
  // -----------------------------------------------------------------------
  void set_timing_cut_width(const int& t) { _range = t; }
  void set_debug(bool debug) { m_debug = debug; }

  // -----------------------------------------------------------------------
  // Utility method to produce log-scale Y-axis TH2 (for QA or debugging)
  // -----------------------------------------------------------------------
  TH2* LogYHist2D(const std::string& name,
                  const std::string& title,
                  int xbins, double xmin, double xmax,
                  int ybins, double ymin, double ymax);


 private:
  // internal method to set up histograms
  void createHistos();
  // utility for histogram reflection
  void MirrorHistogram(TH1* histogram);
  // prefix used by QA hist manager
  std::string getHistoPrefix() const;
  TriggerAnalyzer* trigAna{nullptr};

  // -------------------------------------------------------------------
  // Data members
  // -------------------------------------------------------------------

  // A 3D histogram for (triggerBit, flattened IB index, pi0Mass)
  TH3F* h_pi0_trigIB_mass = nullptr;

  // Additional triggers of interest
  std::vector<int> triggerIndices = {10, 16, 17, 18, 19, 24, 25, 26, 27};

  // Tower-level pedestal + energy histos
  TH1* h_cemc_channel_pedestal[128 * 192]{nullptr};
  TH1* h_ihcal_channel_pedestal[32 * 48]{nullptr};
  TH1* h_ohcal_channel_pedestal[32 * 48]{nullptr};
  TH1* h_cemc_channel_energy[128 * 192]{nullptr};
  TH1* h_ihcal_channel_energy[32 * 48]{nullptr};
  TH1* h_ohcal_channel_energy[32 * 48]{nullptr};

  // QAHists
  TH2* h_emcal_mbd_correlation{nullptr};
  TH1* h_mbd_hits{nullptr};
  TH2* h_ohcal_mbd_correlation{nullptr};
  TH2* h_ihcal_mbd_correlation{nullptr};
  TH2* h_emcal_hcal_correlation{nullptr};

  // EMCal, HCal, ...
  TH2* h_cemc_etaphi{nullptr};
  TH2* h_ihcal_etaphi{nullptr};
  TH2* h_ohcal_etaphi{nullptr};

  TH2* h_cemc_etaphi_wQA{nullptr};
  TH2* h_ihcal_etaphi_wQA{nullptr};
  TH2* h_ohcal_etaphi_wQA{nullptr};

  TH1* h_ihcal_status{nullptr};
  TH1* h_ohcal_status{nullptr};
  TH1* h_cemc_status{nullptr};

  TH2* h_cemc_e_chi2{nullptr};
  TH2* h_ihcal_e_chi2{nullptr};
  TH2* h_ohcal_e_chi2{nullptr};

  // time profiles
  TProfile2D* h_cemc_etaphi_time{nullptr};
  TProfile2D* h_ihcal_etaphi_time{nullptr};
  TProfile2D* h_ohcal_etaphi_time{nullptr};

  // fracHit
  TProfile2D* h_cemc_etaphi_fracHitADC{nullptr};
  TProfile2D* h_cemc_etaphi_fracHit{nullptr};
  TProfile2D* h_ihcal_etaphi_fracHitADC{nullptr};
  TProfile2D* h_ohcal_etaphi_fracHitADC{nullptr};

  // pedRMS
  TProfile2D* h_cemc_etaphi_pedRMS{nullptr};
  TProfile2D* h_ihcal_etaphi_pedRMS{nullptr};
  TProfile2D* h_ohcal_etaphi_pedRMS{nullptr};

  // ZSpedRMS
  TProfile2D* h_cemc_etaphi_ZSpedRMS{nullptr};
  TProfile2D* h_ohcal_etaphi_ZSpedRMS{nullptr};
  TProfile2D* h_ihcal_etaphi_ZSpedRMS{nullptr};

  // Chi2 profiles
  TProfile2D* h_cemc_etaphi_badChi2{nullptr};
  TProfile2D* h_ihcal_etaphi_badChi2{nullptr};
  TProfile2D* h_ohcal_etaphi_badChi2{nullptr};

  // Basic inv mass
  TH1* h_InvMass{nullptr};
  TH1* h_channel_pedestal_0{nullptr};

  // vertex, timing, tower energy, cluster, etc.
  TH1* h_vtx_z_raw{nullptr};
  TH1* h_vtx_z_cut{nullptr};
  TH1* h_emcaltime_cut{nullptr};
  TH1* h_ihcaltime_cut{nullptr};
  TH1* h_ohcaltime_cut{nullptr};
  TH1* h_emcaltime{nullptr};
  TH1* h_ihcaltime{nullptr};
  TH1* h_ohcaltime{nullptr};
  TH1* h_emcal_tower_e{nullptr};
  TH2* h_etaphi_clus{nullptr};
  TH1* h_clusE{nullptr};

  TProfile2D* h_cemc_etaphi_time_raw{nullptr};
  TProfile2D* h_ohcal_etaphi_time_raw{nullptr};
  TProfile2D* h_ihcal_etaphi_time_raw{nullptr};

  // Trigger histos
  TH1 *h_triggerVec{nullptr};
  TH2 *h_edist[64] = {nullptr};
  TH1 *h_ldClus_trig[64] = {nullptr};
  TProfile *pr_evtNum_ldClus_trig[64] = {nullptr};
  TProfile *pr_rejection[64] = {nullptr};
  TProfile *pr_livetime[64] = {nullptr};
  TProfile *pr_ldClus_trig{nullptr};
  std::vector<int> trigOfInterest = {3,10,11,21,22,23,25,26,27};

  // Additional example histos
  TH1* h_InvMass_trig10 {nullptr};
  TH2F* h_pi0_mass_vs_IB {nullptr};

  // internal counters
  int _eventcounter{0};
  int _range{1};
  bool m_debug{false};

  // Output management
  std::string m_outputFileName;
  std::string OutputFileName;
};

#endif  // CALOVALID_CALOVALID_H
