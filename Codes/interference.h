#ifndef INTERFERENCE_H
#define INTERFERENCE_H

#include "kloe_class.h"
#include "../const.h"
#include "TMath.h"
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include <string>
#include <iostream>
#include <TGraphErrors.h>
#include <TFile.h>
#include <vector>

namespace KLOE
{
  class interference : public pm00
  {
    public:

      std::vector<Double_t> time_diff[chann_num], time_diff_gen;
      std::vector<Double_t> time_diff_rand_mc[chann_num], time_diff_rand_data[chann_num], time_diff_gen_rand_mc, time_diff_gen_rand_data;

      //////////////////////////////////////////////////////////////////////////////////////////////////////

      interference(TString mode_init, UInt_t *evnum_init, UInt_t bin_num_init, Double_t x_min_init, Double_t x_max_init) : mode(mode_init)
      {

        bin_number = bin_num_init;
        x_min = x_min_init;
        x_max = x_max_init;

        evnum = new UInt_t[chann_num];

        for(Int_t i = 0; i < chann_num; i++) evnum[i] = evnum_init[i];

        //////////////////////////////////////////////////////////////////////////////////

        left_x_split = -30.0;
        center_x_split = 0.0;
        right_x_split = 30.0;

        //////////////////////////////////////////////////////////////////////////////////

        for (Int_t i = 0; i < chann_num; i++)
          frac[i] = new TH1F(("Fitted histo " + std::to_string(i)).c_str(), "", bin_number, x_min, x_max);

        //////////////////////////////////////////////////////////////////////////////////////////////////

        corr_vals = new Double_t[bin_number];

        TFile file("/internal/big_one/4/users/gamrat/scripts/Scripts/corr_file.root");

        corr_factor = (TGraphErrors *)file.Get("correction_factor");
        corr_vals = corr_factor->GetY();

        file.Close();

        ///////////////////////////////////////////////////////////////////////////////////
      };

      //////////////////////////////////////////////////////////////////////////////////////////////////////

      void bin_extraction(UInt_t channel, TH1 *histogram);

      void mc_randomize();

      //////////////////////////////////////////////////////////////////////////////////////////////////////

      Double_t interf_function(const Float_t x, Int_t check = 1, const Double_t *par = 0);

      Double_t interf_chi2_split(const Double_t *xx);
      Double_t interf_chi2_window(const Double_t *xx);
      Double_t interf_chi2_excluded(const Double_t *xx);
      Double_t interf_chi2_all(const Double_t *xx);

      Double_t interf_chi2_mc(const Double_t *xx);
      Double_t interf_chi2_bcg(const Double_t *xx);

      Double_t interf_chi2(const Double_t *xx);

      /////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////

    private:

      TString mode = "split"; // "split", "window", "excluded, "mc", "bcg", "all"
      TGraphErrors *corr_factor;

      /////////////////////////////////////////////////////////////////////////////////////////////////////

      std::vector<UInt_t> indices[chann_num];
      std::vector<Double_t> b[chann_num], e[chann_num];

      /////////////////////////////////////////////////////////////////////////////////////////////////////

      UInt_t *evnum;
      Double_t *exclusions, left_x_split, center_x_split, right_x_split;
      Double_t *corr_vals;

  };
} // namespace KLOE

#endif