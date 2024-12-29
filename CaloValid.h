#ifndef CALOVALID_CALOVALID_H
#define CALOVALID_CALOVALID_H

#include <fun4all/SubsysReco.h>
#include <bitset>
#include <string>
#include <vector>
#include <iostream>  // for std::cout, std::endl
#include <TH1F.h>    // for TH1F forward usage

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

  // -------------------------------------------------------------------
  //   Inline helpers for extracting or checking trigger bits
  // -------------------------------------------------------------------
  inline std::vector<int> extractTriggerBits(uint64_t gl1_scaledvec, int entry)
  {
    std::vector<int> trig_bits;
    std::bitset<64> bits(gl1_scaledvec);

    if (m_debug)
    {
      std::cout << "[DEBUG] Processing entry " << entry
                << " bits: " << bits.to_string() << std::endl;
    }

    for (unsigned int bit = 0; bit < 64; bit++)
    {
      if (((gl1_scaledvec >> bit) & 0x1U) == 0x1U)
      {
        trig_bits.push_back(bit);
      }
    }
    return trig_bits;
  }

  inline bool checkTriggerCondition(const std::vector<int> &trig_bits, int inputBit)
  {
    for (const int &bit : trig_bits)
    {
      if (bit == inputBit)
      {
        if (m_debug)
        {
          std::cout << "[DEBUG] Trigger bit " << bit << " fired.\n";
        }
        return true;
      }
    }
    if (m_debug)
    {
      std::cout << "[DEBUG] No relevant trigger condition met.\n";
    }
    return false;
  }

  // -------------------------------------------------------------------
  //   Static mapping functions for sector + IB determination
  // -------------------------------------------------------------------
  static int custom_sector_mapping(unsigned int eta, unsigned int phi)
  {
    // If eta >= 48..96 => sector = (phi / 8) in [0..31]
    // If eta < 48      => sector = 32 + (phi / 8) in [32..63]
    // If sector < 0 => out-of-range => handle as needed
    int sector = -1;
    if (eta >= 48 && eta < 96)
    {
      sector = (phi / 8); // 0..31
    }
    else if (eta < 48)
    {
      sector = 32 + (phi / 8); // 32..63
    }
    return sector;
  }

  static int custom_ib_board(int eta, int phi)
  {
    // IB board determination (6 boards each region):
    //   If eta >= 48..56 => ib_board=0
    //   56..64 => 1, 64..72 => 2, 72..80 => 3, 80..88 => 4, 88..96 => 5
    //   (lower region) 40..48 =>0, 32..40 =>1, 24..32 =>2, etc.
    int ib_board;
    if      (eta >= 48 && eta < 56) { ib_board = 0; }
    else if (eta >= 56 && eta < 64) { ib_board = 1; }
    else if (eta >= 64 && eta < 72) { ib_board = 2; }
    else if (eta >= 72 && eta < 80) { ib_board = 3; }
    else if (eta >= 80 && eta < 88) { ib_board = 4; }
    else if (eta >= 88 && eta < 96) { ib_board = 5; }
    else if (eta >= 40 && eta < 48) { ib_board = 0; }
    else if (eta >= 32 && eta < 40) { ib_board = 1; }
    else if (eta >= 24 && eta < 32) { ib_board = 2; }
    else if (eta >= 16 && eta < 24) { ib_board = 3; }
    else if (eta >= 8  && eta < 16) { ib_board = 4; }
    else if (eta >= 0  && eta < 8)  { ib_board = 5; }
    else
    {
      ib_board = -1;  // out-of-range
    }
    return ib_board;
  }

 private:
  // internal method to set up histograms
  void createHistos();
  // utility for histogram reflection
  void MirrorHistogram(TH1* histogram);
  // prefix used by QA hist manager
  std::string getHistoPrefix() const;

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
