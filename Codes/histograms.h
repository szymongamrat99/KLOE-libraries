#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

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
  class histograms
  {
    
    private:
      std::vector<TH1*> hist; 

    public:
      histograms()
      {
        hist.push_back(new TH1D(""));
      };

  }
}

#endif