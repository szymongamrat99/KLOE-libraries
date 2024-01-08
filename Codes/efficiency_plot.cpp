#include "efficiency_plot.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "../const.h"

TMultiGraph *efficiency_plot(Float_t *variable, Float_t *efficiency, TGraph *graphs[chann_num-2])
{
    TMultiGraph *merged;

    for(Int_t i = 0; i < chann_num - 2; i++)
        merged->Add(graphs[i]);
}