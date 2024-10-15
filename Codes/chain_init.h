#ifndef CHAIN_INIT_H
#define CHAIN_INIT_H

#include <TString.h>
#include <TChain.h>

#include <const.h>

inline void chain_init(TChain *chain_init, UInt_t first, UInt_t last)
{
    //char *env = "KLOEFILES";

    TString path(path_tmp);

    TString fullname = "",
    dirnamemc = "MONTE_CARLO", dirnamedata = "DATA",
    filenamemc = "mc_stream62_mccard2_", filenamedata = "data_stream42_",
    extension = ".root";

    for(Int_t i = first; i <= last; i++)
    {
    	fullname = path + "/" + dirnamedata + "/" + filenamedata + i + extension;
    	chain_init->Add(fullname);
    }
    for(Int_t i = first; i <= last; i++)
    {

        fullname = path + "/" + dirnamemc + "/" + filenamemc + i + extension;
    	chain_init->Add(fullname);
    }
};

#endif // !CHAIN_INIT_H