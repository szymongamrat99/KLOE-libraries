#ifndef EFFICIENCY_PLOT_H
#define EFFICIENCY_PLOT_H

#include <iostream>
#include <vector>

#include <TMultiGraph.h>
#include <TGraph.h>
#include <TMath.h>

class EfficiencyPlots
{
  private:
    TMultiGraph *gathered_eff;
    TGraph *eff_graphs[chann_num];
    std::vector<Double_t> eff_vals_gr[chann_num];
    std::vector<Double_t> eff_vals_sm[chann_num];
    std::vector<Double_t> cut_var[chann_num];
    TString var_name;


  public:
    EfficiencyPlots(Int_t point_num)
    {
      for(Int_t i = 0; i < chann_num; i++)
      {
        eff_graphs[i] = new TGraph(point_num);
      }

      gathered_eff = new TMultiGraph("multi_eff", "");
    }

    void CutVarSet(Double_t cut_var_single, Int_t channel);
    void CutVarClean(Int_t channel);

    void EffCalc(TString cut_var_name);
};

#endif // !EFFICIENCY_PLOT_H