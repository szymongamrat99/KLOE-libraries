#include "efficiency_plot.h"
#include "../const.h"

void EfficiencyPlots::CutVarSet(Double_t cut_var_single, Int_t channel)
{
    cut_var[channel].push_back(cut_var_single);
}

void EfficiencyPlots::CutVarClean(Int_t channel)
{
    cut_var[channel].clear();
}

Bool_t CutVarCheck()
{
    return 
}

void EfficiencyPlots::EffCalc(TString cut_var_name)
{
    Int_t tot_el = 0;
    for(Int_t i = 0; i < channNum; i++)
    {
        tot_el = cut_var[i].size();

        std::count_if(cut_var[i].begin(), cut_var[i].end(), );
    }
}