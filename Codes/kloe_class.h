#ifndef KLOE_CLASS_H
#define KLOE_CLASS_H

#include <TMath.h>
#include <TLorentzVector.h>
#include <TH1.h>

namespace KLOE
{
    class pm00
    {
        public:
            TLorentzVector phi_mom, phi_pos;
            TLorentzVector kaon_mom[2], kaon_pos[2];
            TLorentzVector pi_ch_mom[2], pi_ch_pos[2];
            TLorentzVector pi_ne_mom[2], pi_ne_pos[2];
            TLorentzVector photon_mom[4], photon_pos[4];

            UInt_t bin_number;
            Double_t x_min, x_max, y_min, y_max;

            Double_t inv_mass;
            TVector3 boost;
            Double_t angle_vec, transv, azim_angle;

            std::vector<TH1*> frac, frac_data;
            TH1 *data, *mc_sum;

            pm00(TLorentzVector *mom_list, TLorentzVector *pos_list);

            pm00();

            void inv_mass_calc(TLorentzVector four_mom);
            void angle(TLorentzVector vec1, TLorentzVector vec2);
            void cyl_comp(TLorentzVector vec);

            void boost_vector(TLorentzVector four_mom);
            void lorentz_transf(TLorentzVector four_mom);

    };
}

#endif