#include "chain_init.h"

void chain_init(TChain *chain_init, UInt_t first, UInt_t last)
{
    //char *env = "KLOEFILES";
    char *path_tmp = "/internal/big_one/4/users/gamrat/old_root_files";//std::getenv(env);

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