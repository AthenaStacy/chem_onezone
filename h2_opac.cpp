#include "h2_opac.dat"

void compute_h2_opacity(double temp, double N_H2_eff, double &H2_opacity_correction)
{

      double column_min, column_max, logN, logT, diff, dN, dT;
      double opac_tmp[2], opac;
      int ii, jj;

      column_min = pow(10,(h2_opac_column[1-1]));
      column_max = pow(10,(h2_opac_column[nh2op-1]));

      if (N_H2_eff <= column_min){
        opac = 1e0;
        H2_opacity_correction = opac;
        return;
        }
      else if(N_H2_eff >= column_max){
        opac = 0e0;
        H2_opacity_correction = opac;
        return;
        }
      else {
        logN = log10(N_H2_eff);
        jj   = 1 + int(10 * (logN - 17.0)) -1;
        diff = h2_opac_column[jj+1] - h2_opac_column[jj];
        dN   = (logN - h2_opac_column[jj]) / diff;

        logT = log10(temp);
        if (logT <= h2_opac_temp[1-1]) {
          ii = 1 -1;
          dT = 0e0;
          }
        else if (logT >= h2_opac_temp[nh2op-1]) {
          ii = nh2op -1;
          dT = 0e0;
          }
        else {
          ii   = 1 + int((logT - 1.5) / 0.03) -1;
          diff = h2_opac_temp[ii+1] - h2_opac_temp[ii];
          dT   = (logT - h2_opac_temp[ii]) / diff;
         }

        if (dT > 0e0) {
          opac_tmp[0] = h2_opac[ii][jj] + dT * (h2_opac[ii+1][jj]
                     - h2_opac[ii][jj]);
          opac_tmp[1] = h2_opac[ii][jj+1] + dT * (h2_opac[ii+1][jj+1]
                     - h2_opac[ii][jj+1]);
         }
        else {
          opac_tmp[0] = h2_opac[ii][jj];
          opac_tmp[1] = h2_opac[ii][jj+1];
         }

        opac = opac_tmp[0] + dN * (opac_tmp[1] - opac_tmp[0]);
        opac = pow(1e1, opac);
        }

       H2_opacity_correction = opac;

      return;

};

