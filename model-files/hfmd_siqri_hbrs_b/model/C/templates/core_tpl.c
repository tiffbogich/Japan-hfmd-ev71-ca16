/**************************************************************************
 *    This file is part of plom.
 *
 *    plom is free software: you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published
 *    by the Free Software Foundation, either version 3 of the
 *    License, or (at your option) any later version.
 *
 *    plom is distributed in the hope that it will be useful, but
 *    WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License along with plom.  If not, see
 *    <http://www.gnu.org/licenses/>.
 *************************************************************************/

#include "plom.h"

/* automatically generated code: order of the parameters */

#define ORDER_SS 0
#define ORDER_IS 1
#define ORDER_SI 2
#define ORDER_RS 3
#define ORDER_SR 4
#define ORDER_IR 5
#define ORDER_RI 6
#define ORDER_SQ 7
#define ORDER_QS 8
#define ORDER_QR 9
#define ORDER_RQ 10
#define ORDER_r0_1 11
#define ORDER_r0_2 12
#define ORDER_v 13
#define ORDER_e 14
#define ORDER_d 15
#define ORDER_sigma 16
#define ORDER_sto 17
#define ORDER_iota_1 18
#define ORDER_iota_2 19
#define ORDER_q 20
#define ORDER_z 21
#define ORDER_rep 22


#define ORDER_U 11
#define ORDER_RR 12




#define ORDER_N 0
#define ORDER_prop 1
#define ORDER_ptest 2
#define ORDER_hfmd 3
#define ORDER_mu_b 4
#define ORDER_mu_d 5

/**
 * Alloc memory for the psr implementation
 */
void build_psr(struct s_calc *p)
{
    unsigned int tab[13]; 

    /*automaticaly generated code: dimension of prob and inc*/
    
    tab[ORDER_SS] = 4;
    tab[ORDER_IS] = 3;
    tab[ORDER_SI] = 3;
    tab[ORDER_RS] = 5;
    tab[ORDER_SR] = 5;
    tab[ORDER_IR] = 3;
    tab[ORDER_RI] = 3;
    tab[ORDER_SQ] = 3;
    tab[ORDER_QS] = 3;
    tab[ORDER_QR] = 3;
    tab[ORDER_RQ] = 3;
    tab[ORDER_U] = 2;
    tab[ORDER_RR] = 5;

    p->prob = init2d_var_set0(13, tab);
    p->inc = init3u_varp2_set0(13, N_CAC, tab);

    //  p->gravity = init1d_set0(N_C);
}


void proj2obs(struct s_X *p_X, struct s_data *p_data)
{
    int o, ind_obs, ind_proj_inc;
    int n_ts_unique_o, n_stream_o_ts;

    struct s_obs2ts **obs2ts = p_data->obs2ts;

    

    ind_obs = 0;
    ind_proj_inc = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    /* extend incidence: duplicate ts with multiple data streams */
    for(o=0; o<N_OBS_INC; o++) {
        for(n_ts_unique_o=0; n_ts_unique_o< (obs2ts[o])->n_ts_unique; n_ts_unique_o++) {
            for(n_stream_o_ts=0; n_stream_o_ts< (obs2ts[o])->n_stream[n_ts_unique_o]; n_stream_o_ts++) {
                p_X->obs[ind_obs++] = p_X->proj[ind_proj_inc];
            }
            ind_proj_inc++;
        }
    }

    /* add prevalence: aggregate across c and ac to match ts and repeat to tacle multiple data streams */
    
}

//stepping functions for Poisson System with stochastic rates (psr)
void step_psr(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    /* t is the time in unit of the data */

    struct s_obs2ts **obs2ts = p_data->obs2ts;  /* syntaxic shortcut */
    struct s_router **routers = p_data->routers;   /* syntaxic shortcut */

    int c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;

    double sum, one_minus_exp_sum;

    double **par = p_par->natural;   

    double *X = p_X->proj;
    double dt = p_X->dt;


    /*automaticaly generated code:*/
    /*0-declaration of noise terms (if any)*/
    
    double white_noise__0;
    double white_noise__1;

    double _r[11];
    
    double _sf[3];

    


    for(c=0;c<N_C;c++) {
        for(ac=0;ac<N_AC;ac++) {
            cac = c*N_AC+ac;

            

            /*1-generate noise increments (if any) (automaticaly generated code)*/
            
            if(p_data->noises_off & PLOM_NO_ENV_STO){
                
                white_noise__0 = 1.0;
                white_noise__1 = 1.0;
            } else {
                
                white_noise__0 = gsl_ran_gamma(p_calc->randgsl, (dt)/ pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]], 2), pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]], 2))/dt;
                white_noise__1 = gsl_ran_gamma(p_calc->randgsl, (dt)/ pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]], 2), pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]], 2))/dt;
            }
            

            /*2-generate process increments (automaticaly generated code)*/
            
            _sf[0] = correct_rate(par[ORDER_q][routers[ORDER_q]->map[cac]],dt);
            _sf[1] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR));
            _sf[2] = correct_rate(par[ORDER_v][routers[ORDER_v]->map[cac]],dt);

            
            _r[0] = _sf[0];
            _r[1] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*white_noise__1*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[1]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
            _r[2] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*white_noise__0*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[1]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
            _r[3] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*white_noise__0*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[1]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
            _r[4] = _sf[2];
            _r[5] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*white_noise__1*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[1]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
            _r[6] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*white_noise__1*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[1]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
            _r[7] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*white_noise__0*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[1]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
            _r[8] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*white_noise__1*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[1]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
            _r[9] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
            _r[10] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*white_noise__0*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[1]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

            sum = _r[9]*dt+_r[2]*dt+_r[6]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_SS][0] = one_minus_exp_sum*((_r[9]*dt)/sum);
p_calc->prob[ORDER_SS][1] = one_minus_exp_sum*((_r[2]*dt)/sum);
p_calc->prob[ORDER_SS][2] = one_minus_exp_sum*((_r[6]*dt)/sum);
p_calc->prob[ORDER_SS][3] = 1.0 - p_calc->prob[ORDER_SS][0] - p_calc->prob[ORDER_SS][1] - p_calc->prob[ORDER_SS][2];
}
else{
p_calc->prob[ORDER_SS][0] = 0.0;
p_calc->prob[ORDER_SS][1] = 0.0;
p_calc->prob[ORDER_SS][2] = 0.0;
p_calc->prob[ORDER_SS][3] = 1.0;
}

sum = _r[9]*dt+_r[4]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_IS][0] = one_minus_exp_sum*((_r[9]*dt)/sum);
p_calc->prob[ORDER_IS][1] = one_minus_exp_sum*((_r[4]*dt)/sum);
p_calc->prob[ORDER_IS][2] = 1.0 - p_calc->prob[ORDER_IS][0] - p_calc->prob[ORDER_IS][1];
}
else{
p_calc->prob[ORDER_IS][0] = 0.0;
p_calc->prob[ORDER_IS][1] = 0.0;
p_calc->prob[ORDER_IS][2] = 1.0;
}

sum = _r[9]*dt+_r[4]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_SI][0] = one_minus_exp_sum*((_r[9]*dt)/sum);
p_calc->prob[ORDER_SI][1] = one_minus_exp_sum*((_r[4]*dt)/sum);
p_calc->prob[ORDER_SI][2] = 1.0 - p_calc->prob[ORDER_SI][0] - p_calc->prob[ORDER_SI][1];
}
else{
p_calc->prob[ORDER_SI][0] = 0.0;
p_calc->prob[ORDER_SI][1] = 0.0;
p_calc->prob[ORDER_SI][2] = 1.0;
}

sum = _r[9]*dt+_r[1]*dt+_r[7]*dt+_r[10]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_RS][0] = one_minus_exp_sum*((_r[9]*dt)/sum);
p_calc->prob[ORDER_RS][1] = one_minus_exp_sum*((_r[1]*dt)/sum);
p_calc->prob[ORDER_RS][2] = one_minus_exp_sum*((_r[7]*dt)/sum);
p_calc->prob[ORDER_RS][3] = one_minus_exp_sum*((_r[10]*dt)/sum);
p_calc->prob[ORDER_RS][4] = 1.0 - p_calc->prob[ORDER_RS][0] - p_calc->prob[ORDER_RS][1] - p_calc->prob[ORDER_RS][2] - p_calc->prob[ORDER_RS][3];
}
else{
p_calc->prob[ORDER_RS][0] = 0.0;
p_calc->prob[ORDER_RS][1] = 0.0;
p_calc->prob[ORDER_RS][2] = 0.0;
p_calc->prob[ORDER_RS][3] = 0.0;
p_calc->prob[ORDER_RS][4] = 1.0;
}

sum = _r[9]*dt+_r[3]*dt+_r[5]*dt+_r[8]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_SR][0] = one_minus_exp_sum*((_r[9]*dt)/sum);
p_calc->prob[ORDER_SR][1] = one_minus_exp_sum*((_r[3]*dt)/sum);
p_calc->prob[ORDER_SR][2] = one_minus_exp_sum*((_r[5]*dt)/sum);
p_calc->prob[ORDER_SR][3] = one_minus_exp_sum*((_r[8]*dt)/sum);
p_calc->prob[ORDER_SR][4] = 1.0 - p_calc->prob[ORDER_SR][0] - p_calc->prob[ORDER_SR][1] - p_calc->prob[ORDER_SR][2] - p_calc->prob[ORDER_SR][3];
}
else{
p_calc->prob[ORDER_SR][0] = 0.0;
p_calc->prob[ORDER_SR][1] = 0.0;
p_calc->prob[ORDER_SR][2] = 0.0;
p_calc->prob[ORDER_SR][3] = 0.0;
p_calc->prob[ORDER_SR][4] = 1.0;
}

sum = _r[9]*dt+_r[4]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_IR][0] = one_minus_exp_sum*((_r[9]*dt)/sum);
p_calc->prob[ORDER_IR][1] = one_minus_exp_sum*((_r[4]*dt)/sum);
p_calc->prob[ORDER_IR][2] = 1.0 - p_calc->prob[ORDER_IR][0] - p_calc->prob[ORDER_IR][1];
}
else{
p_calc->prob[ORDER_IR][0] = 0.0;
p_calc->prob[ORDER_IR][1] = 0.0;
p_calc->prob[ORDER_IR][2] = 1.0;
}

sum = _r[9]*dt+_r[4]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_RI][0] = one_minus_exp_sum*((_r[9]*dt)/sum);
p_calc->prob[ORDER_RI][1] = one_minus_exp_sum*((_r[4]*dt)/sum);
p_calc->prob[ORDER_RI][2] = 1.0 - p_calc->prob[ORDER_RI][0] - p_calc->prob[ORDER_RI][1];
}
else{
p_calc->prob[ORDER_RI][0] = 0.0;
p_calc->prob[ORDER_RI][1] = 0.0;
p_calc->prob[ORDER_RI][2] = 1.0;
}

sum = _r[9]*dt+_r[0]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_SQ][0] = one_minus_exp_sum*((_r[9]*dt)/sum);
p_calc->prob[ORDER_SQ][1] = one_minus_exp_sum*((_r[0]*dt)/sum);
p_calc->prob[ORDER_SQ][2] = 1.0 - p_calc->prob[ORDER_SQ][0] - p_calc->prob[ORDER_SQ][1];
}
else{
p_calc->prob[ORDER_SQ][0] = 0.0;
p_calc->prob[ORDER_SQ][1] = 0.0;
p_calc->prob[ORDER_SQ][2] = 1.0;
}

sum = _r[9]*dt+_r[0]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_QS][0] = one_minus_exp_sum*((_r[9]*dt)/sum);
p_calc->prob[ORDER_QS][1] = one_minus_exp_sum*((_r[0]*dt)/sum);
p_calc->prob[ORDER_QS][2] = 1.0 - p_calc->prob[ORDER_QS][0] - p_calc->prob[ORDER_QS][1];
}
else{
p_calc->prob[ORDER_QS][0] = 0.0;
p_calc->prob[ORDER_QS][1] = 0.0;
p_calc->prob[ORDER_QS][2] = 1.0;
}

sum = _r[9]*dt+_r[0]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_QR][0] = one_minus_exp_sum*((_r[9]*dt)/sum);
p_calc->prob[ORDER_QR][1] = one_minus_exp_sum*((_r[0]*dt)/sum);
p_calc->prob[ORDER_QR][2] = 1.0 - p_calc->prob[ORDER_QR][0] - p_calc->prob[ORDER_QR][1];
}
else{
p_calc->prob[ORDER_QR][0] = 0.0;
p_calc->prob[ORDER_QR][1] = 0.0;
p_calc->prob[ORDER_QR][2] = 1.0;
}

sum = _r[9]*dt+_r[0]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_RQ][0] = one_minus_exp_sum*((_r[9]*dt)/sum);
p_calc->prob[ORDER_RQ][1] = one_minus_exp_sum*((_r[0]*dt)/sum);
p_calc->prob[ORDER_RQ][2] = 1.0 - p_calc->prob[ORDER_RQ][0] - p_calc->prob[ORDER_RQ][1];
}
else{
p_calc->prob[ORDER_RQ][0] = 0.0;
p_calc->prob[ORDER_RQ][1] = 0.0;
p_calc->prob[ORDER_RQ][2] = 1.0;
}



            /*3-multinomial drawn (automaticaly generated code)*/
            
            plom_ran_multinomial(p_calc->randgsl, 4, (unsigned int) X[ORDER_SS*N_CAC+cac], p_calc->prob[ORDER_SS], p_calc->inc[ORDER_SS][cac]);
            plom_ran_multinomial(p_calc->randgsl, 3, (unsigned int) X[ORDER_IS*N_CAC+cac], p_calc->prob[ORDER_IS], p_calc->inc[ORDER_IS][cac]);
            plom_ran_multinomial(p_calc->randgsl, 3, (unsigned int) X[ORDER_SI*N_CAC+cac], p_calc->prob[ORDER_SI], p_calc->inc[ORDER_SI][cac]);
            plom_ran_multinomial(p_calc->randgsl, 5, (unsigned int) X[ORDER_RS*N_CAC+cac], p_calc->prob[ORDER_RS], p_calc->inc[ORDER_RS][cac]);
            plom_ran_multinomial(p_calc->randgsl, 5, (unsigned int) X[ORDER_SR*N_CAC+cac], p_calc->prob[ORDER_SR], p_calc->inc[ORDER_SR][cac]);
            plom_ran_multinomial(p_calc->randgsl, 3, (unsigned int) X[ORDER_IR*N_CAC+cac], p_calc->prob[ORDER_IR], p_calc->inc[ORDER_IR][cac]);
            plom_ran_multinomial(p_calc->randgsl, 3, (unsigned int) X[ORDER_RI*N_CAC+cac], p_calc->prob[ORDER_RI], p_calc->inc[ORDER_RI][cac]);
            plom_ran_multinomial(p_calc->randgsl, 3, (unsigned int) X[ORDER_SQ*N_CAC+cac], p_calc->prob[ORDER_SQ], p_calc->inc[ORDER_SQ][cac]);
            plom_ran_multinomial(p_calc->randgsl, 3, (unsigned int) X[ORDER_QS*N_CAC+cac], p_calc->prob[ORDER_QS], p_calc->inc[ORDER_QS][cac]);
            plom_ran_multinomial(p_calc->randgsl, 3, (unsigned int) X[ORDER_QR*N_CAC+cac], p_calc->prob[ORDER_QR], p_calc->inc[ORDER_QR][cac]);
            plom_ran_multinomial(p_calc->randgsl, 3, (unsigned int) X[ORDER_RQ*N_CAC+cac], p_calc->prob[ORDER_RQ], p_calc->inc[ORDER_RQ][cac]);

            /*4-update state variables (automaticaly generated code)*/
            //use inc to cache the Poisson draw as thew might be re-used for the incidence computation
            
            p_calc->inc[ORDER_U][cac][0] = gsl_ran_poisson(p_calc->randgsl, (gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]))*dt);
            p_calc->inc[ORDER_RR][cac][0] = gsl_ran_poisson(p_calc->randgsl, ((par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR))+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]))*white_noise__0)*dt);
            p_calc->inc[ORDER_RR][cac][1] = gsl_ran_poisson(p_calc->randgsl, ((par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR))+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]))*white_noise__1)*dt);
            p_calc->inc[ORDER_RR][cac][2] = gsl_ran_poisson(p_calc->randgsl, ((par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR))+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]))*white_noise__0)*dt);
            p_calc->inc[ORDER_RR][cac][3] = gsl_ran_poisson(p_calc->randgsl, ((par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR))+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]))*white_noise__1)*dt);

            X[ORDER_SS*N_CAC+cac] = p_calc->inc[ORDER_SS][cac][3] + p_calc->inc[ORDER_U][cac][0];
X[ORDER_IS*N_CAC+cac] = p_calc->inc[ORDER_IS][cac][2] + p_calc->inc[ORDER_SS][cac][1] + p_calc->inc[ORDER_RS][cac][2];
X[ORDER_SI*N_CAC+cac] = p_calc->inc[ORDER_SI][cac][2] + p_calc->inc[ORDER_SS][cac][2] + p_calc->inc[ORDER_SR][cac][2];
X[ORDER_RS*N_CAC+cac] = p_calc->inc[ORDER_RS][cac][4] + p_calc->inc[ORDER_QS][cac][1];
X[ORDER_SR*N_CAC+cac] = p_calc->inc[ORDER_SR][cac][4] + p_calc->inc[ORDER_SQ][cac][1];
X[ORDER_IR*N_CAC+cac] = p_calc->inc[ORDER_IR][cac][2] + p_calc->inc[ORDER_SR][cac][1] + p_calc->inc[ORDER_RR][cac][0];
X[ORDER_RI*N_CAC+cac] = p_calc->inc[ORDER_RI][cac][2] + p_calc->inc[ORDER_RS][cac][1] + p_calc->inc[ORDER_RR][cac][1];
X[ORDER_SQ*N_CAC+cac] = p_calc->inc[ORDER_SQ][cac][2] + p_calc->inc[ORDER_SI][cac][1] + p_calc->inc[ORDER_SR][cac][3];
X[ORDER_QS*N_CAC+cac] = p_calc->inc[ORDER_QS][cac][2] + p_calc->inc[ORDER_IS][cac][1] + p_calc->inc[ORDER_RS][cac][3];
X[ORDER_QR*N_CAC+cac] = p_calc->inc[ORDER_QR][cac][2] + p_calc->inc[ORDER_IR][cac][1] + p_calc->inc[ORDER_RR][cac][2];
X[ORDER_RQ*N_CAC+cac] = p_calc->inc[ORDER_RQ][cac][2] + p_calc->inc[ORDER_RI][cac][1] + p_calc->inc[ORDER_RR][cac][3];


        }/*end for on ac*/
    } /*end for on c*/

    /*compute incidence:integral between t and t+1 (automaticaly generated code)*/

    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;
    
    o = 0;

    for(ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for(n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc +=  + p_calc->inc[ORDER_IS][cac][1] + p_calc->inc[ORDER_IR][cac][1];
        }
        X[offset] += sum_inc;
        offset++;
    }
    
    o = 1;

    for(ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for(n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc +=  + p_calc->inc[ORDER_SI][cac][1] + p_calc->inc[ORDER_RI][cac][1];
        }
        X[offset] += sum_inc;
        offset++;
    }
    
}



//stepping functions for ODE and SDEs



int step_ode(double t, const double X[], double f[], void *params)

{

    
    struct s_calc *p_calc = (struct s_calc *) params;
    struct s_data *p_data = p_calc->p_data;
    struct s_par *p_par = p_calc->p_par;
    

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;

    int i, c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;
    
    double **par = p_par->natural;

    double _r[N_CAC][16];

    
    double _sf[N_CAC][1];

    

    


    for(cac=0;cac<N_CAC;cac++){
        


        
        _sf[cac][0] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR));

        
        _r[cac][0] = (par[ORDER_q][routers[ORDER_q]->map[cac]]);
        _r[cac][1] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][3] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][4] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][5] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][6] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][7] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);
        _r[cac][8] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][9] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][10] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][11] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][12] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][13] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][14] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][15] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        
    }

    for(c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            /*automaticaly generated code:*/
            /*ODE system*/
            
            f[0*N_CAC+cac] =  - (_r[cac][13]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][4]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][2]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][8]);
            f[1*N_CAC+cac] =  - (_r[cac][13]*X[ORDER_IS*N_CAC+cac]) - (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][4]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][14]*X[ORDER_RS*N_CAC+cac]);
            f[2*N_CAC+cac] =  - (_r[cac][13]*X[ORDER_SI*N_CAC+cac]) - (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][2]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][15]*X[ORDER_SR*N_CAC+cac]);
            f[3*N_CAC+cac] =  - (_r[cac][13]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][3]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][14]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][12]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][0]*X[ORDER_QS*N_CAC+cac]);
            f[4*N_CAC+cac] =  - (_r[cac][13]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][6]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][15]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][1]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]);
            f[5*N_CAC+cac] =  - (_r[cac][13]*X[ORDER_IR*N_CAC+cac]) - (_r[cac][7]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][6]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][9]);
            f[6*N_CAC+cac] =  - (_r[cac][13]*X[ORDER_RI*N_CAC+cac]) - (_r[cac][7]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][3]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][11]);
            f[7*N_CAC+cac] =  - (_r[cac][13]*X[ORDER_SQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]) + (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][1]*X[ORDER_SR*N_CAC+cac]);
            f[8*N_CAC+cac] =  - (_r[cac][13]*X[ORDER_QS*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QS*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][12]*X[ORDER_RS*N_CAC+cac]);
            f[9*N_CAC+cac] =  - (_r[cac][13]*X[ORDER_QR*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QR*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][5]);
            f[10*N_CAC+cac] =  - (_r[cac][13]*X[ORDER_RQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_RQ*N_CAC+cac]) + (_r[cac][7]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][10]);
        }
    }


    //TODO: drift of the diffusion
    //for(i=N_PAR_SV*N_CAC; i<(N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot); i++){
    //    f[i] = 0.0;
    //}

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/
    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    
    o = 0;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc +=  + (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IR*N_CAC+cac]);
        }

        f[offset] = sum_inc;
        offset++;
    }
    
    o = 1;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc +=  + (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][7]*X[ORDER_RI*N_CAC+cac]);
        }

        f[offset] = sum_inc;
        offset++;
    }
    

    
    return GSL_SUCCESS;
    

}


void step_sde_no_env_sto(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)

{

    
    double *X = p_X->proj;
    double dt = p_X->dt;
    double *f = p_calc->y_pred;
    

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;

    int i, c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;
    
    double **par = p_par->natural;

    double _r[N_CAC][16];

    
    double _sf[N_CAC][1];

    
    double dem_sto__0[N_CAC];
    double dem_sto__1[N_CAC];
    double dem_sto__2[N_CAC];
    double dem_sto__3[N_CAC];
    double dem_sto__4[N_CAC];
    double dem_sto__5[N_CAC];
    double dem_sto__6[N_CAC];
    double dem_sto__7[N_CAC];
    double dem_sto__8[N_CAC];
    double dem_sto__9[N_CAC];
    double dem_sto__10[N_CAC];
    double dem_sto__11[N_CAC];
    double dem_sto__12[N_CAC];
    double dem_sto__13[N_CAC];
    double dem_sto__14[N_CAC];
    double dem_sto__15[N_CAC];
    double dem_sto__16[N_CAC];
    double dem_sto__17[N_CAC];
    double dem_sto__18[N_CAC];
    double dem_sto__19[N_CAC];
    double dem_sto__20[N_CAC];
    double dem_sto__21[N_CAC];
    double dem_sto__22[N_CAC];
    double dem_sto__23[N_CAC];
    double dem_sto__24[N_CAC];
    double dem_sto__25[N_CAC];
    double dem_sto__26[N_CAC];
    double dem_sto__27[N_CAC];
    double dem_sto__28[N_CAC];
    double dem_sto__29[N_CAC];
    double dem_sto__30[N_CAC];
    double dem_sto__31[N_CAC];

    


    for(cac=0;cac<N_CAC;cac++){
        


        
        _sf[cac][0] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR));

        
        _r[cac][0] = (par[ORDER_q][routers[ORDER_q]->map[cac]]);
        _r[cac][1] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][3] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][4] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][5] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][6] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][7] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);
        _r[cac][8] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][9] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][10] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][11] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][12] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][13] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][14] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][15] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        
        dem_sto__0[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__1[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__2[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__3[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__4[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__5[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__6[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__7[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__8[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__9[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__10[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__11[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__12[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__13[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__14[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__15[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__16[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__17[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__18[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__19[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__20[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__21[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__22[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__23[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__24[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__25[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__26[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__27[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__28[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__29[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__30[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__31[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
    }

    for(c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            /*automaticaly generated code:*/
            /*ODE system*/
            
            f[0*N_CAC+cac] = X[0*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][4]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][2]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][8]))*dt + - sqrt((_r[cac][13]*X[ORDER_SS*N_CAC+cac]))*dem_sto__1[cac]- sqrt((_r[cac][4]*X[ORDER_SS*N_CAC+cac]))*dem_sto__12[cac]- sqrt((_r[cac][2]*X[ORDER_SS*N_CAC+cac]))*dem_sto__13[cac]+ sqrt((_r[cac][8]))*dem_sto__0[cac];
            f[1*N_CAC+cac] = X[1*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_IS*N_CAC+cac]) - (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][4]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][14]*X[ORDER_RS*N_CAC+cac]))*dt + - sqrt((_r[cac][13]*X[ORDER_IS*N_CAC+cac]))*dem_sto__2[cac]- sqrt((_r[cac][7]*X[ORDER_IS*N_CAC+cac]))*dem_sto__16[cac]+ sqrt((_r[cac][4]*X[ORDER_SS*N_CAC+cac]))*dem_sto__12[cac]+ sqrt((_r[cac][14]*X[ORDER_RS*N_CAC+cac]))*dem_sto__24[cac];
            f[2*N_CAC+cac] = X[2*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SI*N_CAC+cac]) - (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][2]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][15]*X[ORDER_SR*N_CAC+cac]))*dt + - sqrt((_r[cac][13]*X[ORDER_SI*N_CAC+cac]))*dem_sto__3[cac]- sqrt((_r[cac][7]*X[ORDER_SI*N_CAC+cac]))*dem_sto__17[cac]+ sqrt((_r[cac][2]*X[ORDER_SS*N_CAC+cac]))*dem_sto__13[cac]+ sqrt((_r[cac][15]*X[ORDER_SR*N_CAC+cac]))*dem_sto__25[cac];
            f[3*N_CAC+cac] = X[3*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][3]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][14]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][12]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][0]*X[ORDER_QS*N_CAC+cac]))*dt + - sqrt((_r[cac][13]*X[ORDER_RS*N_CAC+cac]))*dem_sto__4[cac]- sqrt((_r[cac][3]*X[ORDER_RS*N_CAC+cac]))*dem_sto__15[cac]- sqrt((_r[cac][14]*X[ORDER_RS*N_CAC+cac]))*dem_sto__24[cac]- sqrt((_r[cac][12]*X[ORDER_RS*N_CAC+cac]))*dem_sto__28[cac]+ sqrt((_r[cac][0]*X[ORDER_QS*N_CAC+cac]))*dem_sto__20[cac];
            f[4*N_CAC+cac] = X[4*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][6]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][15]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][1]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]))*dt + - sqrt((_r[cac][13]*X[ORDER_SR*N_CAC+cac]))*dem_sto__5[cac]- sqrt((_r[cac][6]*X[ORDER_SR*N_CAC+cac]))*dem_sto__14[cac]- sqrt((_r[cac][15]*X[ORDER_SR*N_CAC+cac]))*dem_sto__25[cac]- sqrt((_r[cac][1]*X[ORDER_SR*N_CAC+cac]))*dem_sto__29[cac]+ sqrt((_r[cac][0]*X[ORDER_SQ*N_CAC+cac]))*dem_sto__21[cac];
            f[5*N_CAC+cac] = X[5*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_IR*N_CAC+cac]) - (_r[cac][7]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][6]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][9]))*dt + - sqrt((_r[cac][13]*X[ORDER_IR*N_CAC+cac]))*dem_sto__6[cac]- sqrt((_r[cac][7]*X[ORDER_IR*N_CAC+cac]))*dem_sto__18[cac]+ sqrt((_r[cac][6]*X[ORDER_SR*N_CAC+cac]))*dem_sto__14[cac]+ sqrt((_r[cac][9]))*dem_sto__26[cac];
            f[6*N_CAC+cac] = X[6*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_RI*N_CAC+cac]) - (_r[cac][7]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][3]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][11]))*dt + - sqrt((_r[cac][13]*X[ORDER_RI*N_CAC+cac]))*dem_sto__7[cac]- sqrt((_r[cac][7]*X[ORDER_RI*N_CAC+cac]))*dem_sto__19[cac]+ sqrt((_r[cac][3]*X[ORDER_RS*N_CAC+cac]))*dem_sto__15[cac]+ sqrt((_r[cac][11]))*dem_sto__27[cac];
            f[7*N_CAC+cac] = X[7*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]) + (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][1]*X[ORDER_SR*N_CAC+cac]))*dt + - sqrt((_r[cac][13]*X[ORDER_SQ*N_CAC+cac]))*dem_sto__8[cac]- sqrt((_r[cac][0]*X[ORDER_SQ*N_CAC+cac]))*dem_sto__21[cac]+ sqrt((_r[cac][7]*X[ORDER_SI*N_CAC+cac]))*dem_sto__17[cac]+ sqrt((_r[cac][1]*X[ORDER_SR*N_CAC+cac]))*dem_sto__29[cac];
            f[8*N_CAC+cac] = X[8*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_QS*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QS*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][12]*X[ORDER_RS*N_CAC+cac]))*dt + - sqrt((_r[cac][13]*X[ORDER_QS*N_CAC+cac]))*dem_sto__9[cac]- sqrt((_r[cac][0]*X[ORDER_QS*N_CAC+cac]))*dem_sto__20[cac]+ sqrt((_r[cac][7]*X[ORDER_IS*N_CAC+cac]))*dem_sto__16[cac]+ sqrt((_r[cac][12]*X[ORDER_RS*N_CAC+cac]))*dem_sto__28[cac];
            f[9*N_CAC+cac] = X[9*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_QR*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QR*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][5]))*dt + - sqrt((_r[cac][13]*X[ORDER_QR*N_CAC+cac]))*dem_sto__10[cac]- sqrt((_r[cac][0]*X[ORDER_QR*N_CAC+cac]))*dem_sto__22[cac]+ sqrt((_r[cac][7]*X[ORDER_IR*N_CAC+cac]))*dem_sto__18[cac]+ sqrt((_r[cac][5]))*dem_sto__30[cac];
            f[10*N_CAC+cac] = X[10*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_RQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_RQ*N_CAC+cac]) + (_r[cac][7]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][10]))*dt + - sqrt((_r[cac][13]*X[ORDER_RQ*N_CAC+cac]))*dem_sto__11[cac]- sqrt((_r[cac][0]*X[ORDER_RQ*N_CAC+cac]))*dem_sto__23[cac]+ sqrt((_r[cac][7]*X[ORDER_RI*N_CAC+cac]))*dem_sto__19[cac]+ sqrt((_r[cac][10]))*dem_sto__31[cac];
        }
    }


    //TODO: drift of the diffusion
    //for(i=N_PAR_SV*N_CAC; i<(N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot); i++){
    //    f[i] = 0.0;
    //}

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/
    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    
    o = 0;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += ( + (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IR*N_CAC+cac]))*dt + sqrt((_r[cac][7]*X[ORDER_IS*N_CAC+cac]))*dem_sto__16[cac]+ sqrt((_r[cac][7]*X[ORDER_IR*N_CAC+cac]))*dem_sto__18[cac];
        }

        f[offset] = X[offset] +  sum_inc;
        offset++;
    }
    
    o = 1;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += ( + (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][7]*X[ORDER_RI*N_CAC+cac]))*dt + sqrt((_r[cac][7]*X[ORDER_SI*N_CAC+cac]))*dem_sto__17[cac]+ sqrt((_r[cac][7]*X[ORDER_RI*N_CAC+cac]))*dem_sto__19[cac];
        }

        f[offset] = X[offset] +  sum_inc;
        offset++;
    }
    

    
    //y_pred (f) -> X (and we ensure that X is > 0.0)
    for(i=0; i<N_PAR_SV*N_CAC; i++){
        X[i] =  (f[i] < 0.0) ? 0.0 : f[i];
    }

    for(i=N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot; i<N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot +N_TS_INC_UNIQUE; i++){
        X[i] = (f[i] < 0.0) ? 0.0 : f[i];
    }
    

}


void step_sde_full(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)

{

    
    double *X = p_X->proj;
    double dt = p_X->dt;
    double *f = p_calc->y_pred;
    

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;

    int i, c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;
    
    double **par = p_par->natural;

    double _r[N_CAC][16];

    
    double _sf[N_CAC][1];

    
    double dem_sto__0[N_CAC];
    double dem_sto__1[N_CAC];
    double dem_sto__2[N_CAC];
    double dem_sto__3[N_CAC];
    double dem_sto__4[N_CAC];
    double dem_sto__5[N_CAC];
    double dem_sto__6[N_CAC];
    double dem_sto__7[N_CAC];
    double dem_sto__8[N_CAC];
    double dem_sto__9[N_CAC];
    double dem_sto__10[N_CAC];
    double dem_sto__11[N_CAC];
    double dem_sto__12[N_CAC];
    double dem_sto__13[N_CAC];
    double dem_sto__14[N_CAC];
    double dem_sto__15[N_CAC];
    double dem_sto__16[N_CAC];
    double dem_sto__17[N_CAC];
    double dem_sto__18[N_CAC];
    double dem_sto__19[N_CAC];
    double dem_sto__20[N_CAC];
    double dem_sto__21[N_CAC];
    double dem_sto__22[N_CAC];
    double dem_sto__23[N_CAC];
    double dem_sto__24[N_CAC];
    double dem_sto__25[N_CAC];
    double dem_sto__26[N_CAC];
    double dem_sto__27[N_CAC];
    double dem_sto__28[N_CAC];
    double dem_sto__29[N_CAC];
    double dem_sto__30[N_CAC];
    double dem_sto__31[N_CAC];
    double white_noise__0[N_CAC];
    double white_noise__1[N_CAC];

    


    for(cac=0;cac<N_CAC;cac++){
        


        
        _sf[cac][0] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR));

        
        _r[cac][0] = (par[ORDER_q][routers[ORDER_q]->map[cac]]);
        _r[cac][1] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][3] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][4] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][5] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][6] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][7] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);
        _r[cac][8] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][9] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][10] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][11] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][12] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][13] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][14] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][15] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        
        dem_sto__0[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__1[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__2[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__3[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__4[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__5[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__6[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__7[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__8[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__9[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__10[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__11[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__12[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__13[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__14[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__15[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__16[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__17[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__18[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__19[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__20[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__21[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__22[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__23[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__24[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__25[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__26[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__27[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__28[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__29[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__30[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__31[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        white_noise__0[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        white_noise__1[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
    }

    for(c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            /*automaticaly generated code:*/
            /*ODE system*/
            
            f[0*N_CAC+cac] = X[0*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][4]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][2]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][8]))*dt + - sqrt((_r[cac][13]*X[ORDER_SS*N_CAC+cac]))*dem_sto__1[cac]- sqrt((_r[cac][4]*X[ORDER_SS*N_CAC+cac]))*dem_sto__12[cac]- sqrt((_r[cac][2]*X[ORDER_SS*N_CAC+cac]))*dem_sto__13[cac]+ sqrt((_r[cac][8]))*dem_sto__0[cac] + - (_r[cac][4]*X[ORDER_SS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac]- (_r[cac][2]*X[ORDER_SS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac];
            f[1*N_CAC+cac] = X[1*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_IS*N_CAC+cac]) - (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][4]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][14]*X[ORDER_RS*N_CAC+cac]))*dt + - sqrt((_r[cac][13]*X[ORDER_IS*N_CAC+cac]))*dem_sto__2[cac]- sqrt((_r[cac][7]*X[ORDER_IS*N_CAC+cac]))*dem_sto__16[cac]+ sqrt((_r[cac][4]*X[ORDER_SS*N_CAC+cac]))*dem_sto__12[cac]+ sqrt((_r[cac][14]*X[ORDER_RS*N_CAC+cac]))*dem_sto__24[cac] + + (_r[cac][4]*X[ORDER_SS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac]+ (_r[cac][14]*X[ORDER_RS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
            f[2*N_CAC+cac] = X[2*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SI*N_CAC+cac]) - (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][2]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][15]*X[ORDER_SR*N_CAC+cac]))*dt + - sqrt((_r[cac][13]*X[ORDER_SI*N_CAC+cac]))*dem_sto__3[cac]- sqrt((_r[cac][7]*X[ORDER_SI*N_CAC+cac]))*dem_sto__17[cac]+ sqrt((_r[cac][2]*X[ORDER_SS*N_CAC+cac]))*dem_sto__13[cac]+ sqrt((_r[cac][15]*X[ORDER_SR*N_CAC+cac]))*dem_sto__25[cac] + + (_r[cac][2]*X[ORDER_SS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac]+ (_r[cac][15]*X[ORDER_SR*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac];
            f[3*N_CAC+cac] = X[3*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][3]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][14]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][12]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][0]*X[ORDER_QS*N_CAC+cac]))*dt + - sqrt((_r[cac][13]*X[ORDER_RS*N_CAC+cac]))*dem_sto__4[cac]- sqrt((_r[cac][3]*X[ORDER_RS*N_CAC+cac]))*dem_sto__15[cac]- sqrt((_r[cac][14]*X[ORDER_RS*N_CAC+cac]))*dem_sto__24[cac]- sqrt((_r[cac][12]*X[ORDER_RS*N_CAC+cac]))*dem_sto__28[cac]+ sqrt((_r[cac][0]*X[ORDER_QS*N_CAC+cac]))*dem_sto__20[cac] + - (_r[cac][3]*X[ORDER_RS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac]- (_r[cac][14]*X[ORDER_RS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac]- (_r[cac][12]*X[ORDER_RS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
            f[4*N_CAC+cac] = X[4*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][6]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][15]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][1]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]))*dt + - sqrt((_r[cac][13]*X[ORDER_SR*N_CAC+cac]))*dem_sto__5[cac]- sqrt((_r[cac][6]*X[ORDER_SR*N_CAC+cac]))*dem_sto__14[cac]- sqrt((_r[cac][15]*X[ORDER_SR*N_CAC+cac]))*dem_sto__25[cac]- sqrt((_r[cac][1]*X[ORDER_SR*N_CAC+cac]))*dem_sto__29[cac]+ sqrt((_r[cac][0]*X[ORDER_SQ*N_CAC+cac]))*dem_sto__21[cac] + - (_r[cac][6]*X[ORDER_SR*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac]- (_r[cac][15]*X[ORDER_SR*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac]- (_r[cac][1]*X[ORDER_SR*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac];
            f[5*N_CAC+cac] = X[5*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_IR*N_CAC+cac]) - (_r[cac][7]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][6]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][9]))*dt + - sqrt((_r[cac][13]*X[ORDER_IR*N_CAC+cac]))*dem_sto__6[cac]- sqrt((_r[cac][7]*X[ORDER_IR*N_CAC+cac]))*dem_sto__18[cac]+ sqrt((_r[cac][6]*X[ORDER_SR*N_CAC+cac]))*dem_sto__14[cac]+ sqrt((_r[cac][9]))*dem_sto__26[cac] + + (_r[cac][6]*X[ORDER_SR*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac]+ (_r[cac][9])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
            f[6*N_CAC+cac] = X[6*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_RI*N_CAC+cac]) - (_r[cac][7]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][3]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][11]))*dt + - sqrt((_r[cac][13]*X[ORDER_RI*N_CAC+cac]))*dem_sto__7[cac]- sqrt((_r[cac][7]*X[ORDER_RI*N_CAC+cac]))*dem_sto__19[cac]+ sqrt((_r[cac][3]*X[ORDER_RS*N_CAC+cac]))*dem_sto__15[cac]+ sqrt((_r[cac][11]))*dem_sto__27[cac] + + (_r[cac][3]*X[ORDER_RS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac]+ (_r[cac][11])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac];
            f[7*N_CAC+cac] = X[7*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]) + (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][1]*X[ORDER_SR*N_CAC+cac]))*dt + - sqrt((_r[cac][13]*X[ORDER_SQ*N_CAC+cac]))*dem_sto__8[cac]- sqrt((_r[cac][0]*X[ORDER_SQ*N_CAC+cac]))*dem_sto__21[cac]+ sqrt((_r[cac][7]*X[ORDER_SI*N_CAC+cac]))*dem_sto__17[cac]+ sqrt((_r[cac][1]*X[ORDER_SR*N_CAC+cac]))*dem_sto__29[cac] + + (_r[cac][1]*X[ORDER_SR*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac];
            f[8*N_CAC+cac] = X[8*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_QS*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QS*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][12]*X[ORDER_RS*N_CAC+cac]))*dt + - sqrt((_r[cac][13]*X[ORDER_QS*N_CAC+cac]))*dem_sto__9[cac]- sqrt((_r[cac][0]*X[ORDER_QS*N_CAC+cac]))*dem_sto__20[cac]+ sqrt((_r[cac][7]*X[ORDER_IS*N_CAC+cac]))*dem_sto__16[cac]+ sqrt((_r[cac][12]*X[ORDER_RS*N_CAC+cac]))*dem_sto__28[cac] + + (_r[cac][12]*X[ORDER_RS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
            f[9*N_CAC+cac] = X[9*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_QR*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QR*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][5]))*dt + - sqrt((_r[cac][13]*X[ORDER_QR*N_CAC+cac]))*dem_sto__10[cac]- sqrt((_r[cac][0]*X[ORDER_QR*N_CAC+cac]))*dem_sto__22[cac]+ sqrt((_r[cac][7]*X[ORDER_IR*N_CAC+cac]))*dem_sto__18[cac]+ sqrt((_r[cac][5]))*dem_sto__30[cac] + + (_r[cac][5])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
            f[10*N_CAC+cac] = X[10*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_RQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_RQ*N_CAC+cac]) + (_r[cac][7]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][10]))*dt + - sqrt((_r[cac][13]*X[ORDER_RQ*N_CAC+cac]))*dem_sto__11[cac]- sqrt((_r[cac][0]*X[ORDER_RQ*N_CAC+cac]))*dem_sto__23[cac]+ sqrt((_r[cac][7]*X[ORDER_RI*N_CAC+cac]))*dem_sto__19[cac]+ sqrt((_r[cac][10]))*dem_sto__31[cac] + + (_r[cac][10])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac];
        }
    }


    //TODO: drift of the diffusion
    //for(i=N_PAR_SV*N_CAC; i<(N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot); i++){
    //    f[i] = 0.0;
    //}

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/
    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    
    o = 0;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += ( + (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IR*N_CAC+cac]))*dt + + sqrt((_r[cac][7]*X[ORDER_IS*N_CAC+cac]))*dem_sto__16[cac]+ sqrt((_r[cac][7]*X[ORDER_IR*N_CAC+cac]))*dem_sto__18[cac] ;
        }

        f[offset] = X[offset] +  sum_inc;
        offset++;
    }
    
    o = 1;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += ( + (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][7]*X[ORDER_RI*N_CAC+cac]))*dt + + sqrt((_r[cac][7]*X[ORDER_SI*N_CAC+cac]))*dem_sto__17[cac]+ sqrt((_r[cac][7]*X[ORDER_RI*N_CAC+cac]))*dem_sto__19[cac] ;
        }

        f[offset] = X[offset] +  sum_inc;
        offset++;
    }
    

    
    //y_pred (f) -> X (and we ensure that X is > 0.0)
    for(i=0; i<N_PAR_SV*N_CAC; i++){
        X[i] =  (f[i] < 0.0) ? 0.0 : f[i];
    }

    for(i=N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot; i<N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot +N_TS_INC_UNIQUE; i++){
        X[i] = (f[i] < 0.0) ? 0.0 : f[i];
    }
    

}


void step_sde_no_dem_sto(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)

{

    
    double *X = p_X->proj;
    double dt = p_X->dt;
    double *f = p_calc->y_pred;
    

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;

    int i, c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;
    
    double **par = p_par->natural;

    double _r[N_CAC][16];

    
    double _sf[N_CAC][1];

    
    double white_noise__0[N_CAC];
    double white_noise__1[N_CAC];

    


    for(cac=0;cac<N_CAC;cac++){
        


        
        _sf[cac][0] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR));

        
        _r[cac][0] = (par[ORDER_q][routers[ORDER_q]->map[cac]]);
        _r[cac][1] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][3] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][4] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][5] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][6] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][7] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);
        _r[cac][8] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][9] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][10] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][11] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][12] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][13] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][14] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][15] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        
        white_noise__0[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        white_noise__1[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
    }

    for(c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            /*automaticaly generated code:*/
            /*ODE system*/
            
            f[0*N_CAC+cac] = X[0*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][4]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][2]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][8]))*dt + - (_r[cac][4]*X[ORDER_SS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac]- (_r[cac][2]*X[ORDER_SS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac];
            f[1*N_CAC+cac] = X[1*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_IS*N_CAC+cac]) - (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][4]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][14]*X[ORDER_RS*N_CAC+cac]))*dt + + (_r[cac][4]*X[ORDER_SS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac]+ (_r[cac][14]*X[ORDER_RS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
            f[2*N_CAC+cac] = X[2*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SI*N_CAC+cac]) - (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][2]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][15]*X[ORDER_SR*N_CAC+cac]))*dt + + (_r[cac][2]*X[ORDER_SS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac]+ (_r[cac][15]*X[ORDER_SR*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac];
            f[3*N_CAC+cac] = X[3*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][3]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][14]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][12]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][0]*X[ORDER_QS*N_CAC+cac]))*dt + - (_r[cac][3]*X[ORDER_RS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac]- (_r[cac][14]*X[ORDER_RS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac]- (_r[cac][12]*X[ORDER_RS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
            f[4*N_CAC+cac] = X[4*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][6]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][15]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][1]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]))*dt + - (_r[cac][6]*X[ORDER_SR*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac]- (_r[cac][15]*X[ORDER_SR*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac]- (_r[cac][1]*X[ORDER_SR*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac];
            f[5*N_CAC+cac] = X[5*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_IR*N_CAC+cac]) - (_r[cac][7]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][6]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][9]))*dt + + (_r[cac][6]*X[ORDER_SR*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac]+ (_r[cac][9])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
            f[6*N_CAC+cac] = X[6*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_RI*N_CAC+cac]) - (_r[cac][7]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][3]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][11]))*dt + + (_r[cac][3]*X[ORDER_RS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac]+ (_r[cac][11])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac];
            f[7*N_CAC+cac] = X[7*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]) + (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][1]*X[ORDER_SR*N_CAC+cac]))*dt + + (_r[cac][1]*X[ORDER_SR*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac];
            f[8*N_CAC+cac] = X[8*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_QS*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QS*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][12]*X[ORDER_RS*N_CAC+cac]))*dt + + (_r[cac][12]*X[ORDER_RS*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
            f[9*N_CAC+cac] = X[9*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_QR*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QR*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][5]))*dt + + (_r[cac][5])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
            f[10*N_CAC+cac] = X[10*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_RQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_RQ*N_CAC+cac]) + (_r[cac][7]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][10]))*dt + + (_r[cac][10])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__1[cac];
        }
    }


    //TODO: drift of the diffusion
    //for(i=N_PAR_SV*N_CAC; i<(N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot); i++){
    //    f[i] = 0.0;
    //}

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/
    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    
    o = 0;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += ( + (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IR*N_CAC+cac]))*dt ;
        }

        f[offset] = X[offset] +  sum_inc;
        offset++;
    }
    
    o = 1;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += ( + (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][7]*X[ORDER_RI*N_CAC+cac]))*dt ;
        }

        f[offset] = X[offset] +  sum_inc;
        offset++;
    }
    

    
    //y_pred (f) -> X (and we ensure that X is > 0.0)
    for(i=0; i<N_PAR_SV*N_CAC; i++){
        X[i] =  (f[i] < 0.0) ? 0.0 : f[i];
    }

    for(i=N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot; i<N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot +N_TS_INC_UNIQUE; i++){
        X[i] = (f[i] < 0.0) ? 0.0 : f[i];
    }
    

}


void step_sde_no_dem_sto_no_env_sto(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)

{

    
    double *X = p_X->proj;
    double dt = p_X->dt;
    double *f = p_calc->y_pred;
    

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;

    int i, c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;
    
    double **par = p_par->natural;

    double _r[N_CAC][16];

    
    double _sf[N_CAC][1];

    

    


    for(cac=0;cac<N_CAC;cac++){
        


        
        _sf[cac][0] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR));

        
        _r[cac][0] = (par[ORDER_q][routers[ORDER_q]->map[cac]]);
        _r[cac][1] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][3] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][4] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][5] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][6] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][7] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);
        _r[cac][8] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][9] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][10] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][11] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][12] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][13] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][14] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][15] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        
    }

    for(c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            /*automaticaly generated code:*/
            /*ODE system*/
            
            f[0*N_CAC+cac] = X[0*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][4]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][2]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][8]))*dt;
            f[1*N_CAC+cac] = X[1*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_IS*N_CAC+cac]) - (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][4]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][14]*X[ORDER_RS*N_CAC+cac]))*dt;
            f[2*N_CAC+cac] = X[2*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SI*N_CAC+cac]) - (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][2]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][15]*X[ORDER_SR*N_CAC+cac]))*dt;
            f[3*N_CAC+cac] = X[3*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][3]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][14]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][12]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][0]*X[ORDER_QS*N_CAC+cac]))*dt;
            f[4*N_CAC+cac] = X[4*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][6]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][15]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][1]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]))*dt;
            f[5*N_CAC+cac] = X[5*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_IR*N_CAC+cac]) - (_r[cac][7]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][6]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][9]))*dt;
            f[6*N_CAC+cac] = X[6*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_RI*N_CAC+cac]) - (_r[cac][7]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][3]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][11]))*dt;
            f[7*N_CAC+cac] = X[7*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_SQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]) + (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][1]*X[ORDER_SR*N_CAC+cac]))*dt;
            f[8*N_CAC+cac] = X[8*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_QS*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QS*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][12]*X[ORDER_RS*N_CAC+cac]))*dt;
            f[9*N_CAC+cac] = X[9*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_QR*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QR*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][5]))*dt;
            f[10*N_CAC+cac] = X[10*N_CAC+cac] +  ( - (_r[cac][13]*X[ORDER_RQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_RQ*N_CAC+cac]) + (_r[cac][7]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][10]))*dt;
        }
    }


    //TODO: drift of the diffusion
    //for(i=N_PAR_SV*N_CAC; i<(N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot); i++){
    //    f[i] = 0.0;
    //}

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/
    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    
    o = 0;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += ( + (_r[cac][7]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][7]*X[ORDER_IR*N_CAC+cac]))*dt;
        }

        f[offset] = X[offset] +  sum_inc;
        offset++;
    }
    
    o = 1;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += ( + (_r[cac][7]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][7]*X[ORDER_RI*N_CAC+cac]))*dt;
        }

        f[offset] = X[offset] +  sum_inc;
        offset++;
    }
    

    
    //y_pred (f) -> X (and we ensure that X is > 0.0)
    for(i=0; i<N_PAR_SV*N_CAC; i++){
        X[i] =  (f[i] < 0.0) ? 0.0 : f[i];
    }

    for(i=N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot; i<N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot +N_TS_INC_UNIQUE; i++){
        X[i] = (f[i] < 0.0) ? 0.0 : f[i];
    }
    

}



double likelihood(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
    /*x is the predicted value from the model that we contrast with a time serie ts.
      Note: user should not use this function but get_log_likelihood
    */

    struct s_router **routers = p_data->routers;

    double like; /* likelihood value */

    double y = p_data->data[n][ts];

    double **par = p_par->natural;

    /*automaticaly generated code*/
    double gsl_mu = par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])*gsl_spline_eval(p_calc->spline[ORDER_ptest][ts],t,p_calc->acc[ORDER_ptest][ts]);
    double gsl_sd = sqrt( par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])*gsl_spline_eval(p_calc->spline[ORDER_ptest][ts],t,p_calc->acc[ORDER_ptest][ts])*(-GSL_MIN(par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])/gsl_spline_eval(p_calc->spline[ORDER_hfmd][ts],t,p_calc->acc[ORDER_hfmd][ts]),1.0)+1.0) );

    if (y > 0.0) {
        like=gsl_cdf_gaussian_P(y+0.5-gsl_mu, gsl_sd)-gsl_cdf_gaussian_P(y-0.5-gsl_mu, gsl_sd);
    } else {
        like=gsl_cdf_gaussian_P(y+0.5-gsl_mu, gsl_sd);
    }

    return sanitize_likelihood(like);
}

double obs_mean(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
  /*x is the predicted value from the model that we contrast with a time serie ts*/
  struct s_router **routers = p_data->routers;

  double **par = p_par->natural;
  
  /*automaticaly generated code*/
  double mu = par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])*gsl_spline_eval(p_calc->spline[ORDER_ptest][ts],t,p_calc->acc[ORDER_ptest][ts]);

  return mu;
}

double obs_var(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
  /*x is the predicted value from the model that we contrast with a time serie ts*/
  struct s_router **routers = p_data->routers;

  double **par = p_par->natural;

  /*automaticaly generated code*/
  double var = par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])*gsl_spline_eval(p_calc->spline[ORDER_ptest][ts],t,p_calc->acc[ORDER_ptest][ts])*(-GSL_MIN(par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])/gsl_spline_eval(p_calc->spline[ORDER_hfmd][ts],t,p_calc->acc[ORDER_hfmd][ts]),1.0)+1.0);

  return var;
}


double observation(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
  /*x is the predicted value from the model that we contrast with a time serie ts*/
  struct s_router **routers = p_data->routers;

  double **par = p_par->natural;  

  /*return an observation of the process model*/

  /*automaticaly generated code*/
  double gsl_mu= par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])*gsl_spline_eval(p_calc->spline[ORDER_ptest][ts],t,p_calc->acc[ORDER_ptest][ts]);
  double gsl_sd=sqrt(par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])*gsl_spline_eval(p_calc->spline[ORDER_ptest][ts],t,p_calc->acc[ORDER_ptest][ts])*(-GSL_MIN(par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])/gsl_spline_eval(p_calc->spline[ORDER_hfmd][ts],t,p_calc->acc[ORDER_hfmd][ts]),1.0)+1.0));

  double yobs= gsl_mu+gsl_ran_gaussian(p_calc->randgsl, gsl_sd);

  return (yobs >0) ? yobs : 0.0;
}
