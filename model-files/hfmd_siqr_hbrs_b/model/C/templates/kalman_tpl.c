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

#include "kalman.h"

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
#define ORDER_rep 21


#define ORDER_U 11
#define ORDER_RR 12




#define ORDER_N 0
#define ORDER_prop 1
#define ORDER_ptest 2
#define ORDER_hfmd 3
#define ORDER_mu_b 4
#define ORDER_mu_d 5


/**
 * the function used by f_prediction_ode_rk:
 * dX/dt = f(t, X, params)
 *
 * it is an expanded copy of core/prediction.c/func
 */
int step_ode_ekf(double t, const double X[], double f[], void *params)
{
    struct s_calc *p_calc = (struct s_calc *) params;
    struct s_par *p_par = p_calc->p_par;  /* syntaxic shortcut */
    struct s_data *p_data = p_calc->p_data;
    struct s_obs2ts **obs2ts = p_data->obs2ts;		/* syntaxic shortcut */
    struct s_router **routers = p_data->routers;	/* syntaxic shortcut */

    struct s_kalman_specific_data *p_kalman_specific_data = (struct s_kalman_specific_data *) p_calc->method_specific_thread_safe_data;

    gsl_matrix *Ft = p_kalman_specific_data->Ft;
    gsl_matrix *Q = p_kalman_specific_data->Q;
    struct s_group ***compo_groups_drift_par_proc = p_kalman_specific_data->compo_groups_drift_par_proc;
    gsl_matrix *FtCt = p_kalman_specific_data->FtCt;

    gsl_matrix_const_view Ct   = gsl_matrix_const_view_array(&X[N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE], N_KAL, N_KAL);
    gsl_matrix_view ff = gsl_matrix_view_array(&f[N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE], N_KAL, N_KAL);

    int i, c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;

    double **par = p_par->natural;

    

    double _r[N_CAC][10];

    
    double _sf[N_CAC][1];

    

    for(cac=0;cac<N_CAC;cac++) {
	

        
        _sf[cac][0] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR));

        
        _r[cac][0] = (par[ORDER_q][routers[ORDER_q]->map[cac]]);
        _r[cac][1] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][3] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][4] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][5] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);
        _r[cac][6] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][7] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][8] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][9] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
    }

    for (c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

	    
	    f[0*N_CAC+cac] =  - (_r[cac][8]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][3]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][1]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][6]);
	    f[1*N_CAC+cac] =  - (_r[cac][8]*X[ORDER_IS*N_CAC+cac]) - (_r[cac][5]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][3]*X[ORDER_SS*N_CAC+cac]);
	    f[2*N_CAC+cac] =  - (_r[cac][8]*X[ORDER_SI*N_CAC+cac]) - (_r[cac][5]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][1]*X[ORDER_SS*N_CAC+cac]);
	    f[3*N_CAC+cac] =  - (_r[cac][8]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][2]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][3]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][0]*X[ORDER_QS*N_CAC+cac]);
	    f[4*N_CAC+cac] =  - (_r[cac][8]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][4]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][1]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]);
	    f[5*N_CAC+cac] =  - (_r[cac][8]*X[ORDER_IR*N_CAC+cac]) - (_r[cac][5]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][4]*X[ORDER_SR*N_CAC+cac]);
	    f[6*N_CAC+cac] =  - (_r[cac][8]*X[ORDER_RI*N_CAC+cac]) - (_r[cac][5]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][2]*X[ORDER_RS*N_CAC+cac]);
	    f[7*N_CAC+cac] =  - (_r[cac][8]*X[ORDER_SQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]) + (_r[cac][5]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][1]*X[ORDER_SR*N_CAC+cac]);
	    f[8*N_CAC+cac] =  - (_r[cac][8]*X[ORDER_QS*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QS*N_CAC+cac]) + (_r[cac][5]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][3]*X[ORDER_RS*N_CAC+cac]);
	    f[9*N_CAC+cac] =  - (_r[cac][8]*X[ORDER_QR*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QR*N_CAC+cac]) + (_r[cac][5]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][7]);
	    f[10*N_CAC+cac] =  - (_r[cac][8]*X[ORDER_RQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_RQ*N_CAC+cac]) + (_r[cac][5]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][9]);
        }
    }

    

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/

    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;
    
    o = 0;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for(n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc +=  + (_r[cac][5]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][5]*X[ORDER_IR*N_CAC+cac]);
        }

        f[offset] = sum_inc;
        offset++;
    }
    
    o = 1;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for(n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc +=  + (_r[cac][5]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][5]*X[ORDER_RI*N_CAC+cac]);
        }

        f[offset] = sum_inc;
        offset++;
    }
    

    ////////////////
    // covariance //
    ////////////////

    // evaluate Q and jacobian
    p_kalman_specific_data->eval_Q(Q, X, p_par, p_data, p_calc, p_kalman_specific_data, t);
    eval_jac(Ft, X, p_par, p_data, p_calc, compo_groups_drift_par_proc, t);

    // compute Ft*Ct+Ct*Ft'+Q 
    //here Ct is symmetrical and transpose(FtCt) == transpose(Ct)transpose(Ft) == Ct transpose(Ft)
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ft, &Ct.matrix, 0.0, FtCt);
    for(i=0; i< ff.matrix.size1; i++){
	for(c=0; c< ff.matrix.size2; c++){
	    gsl_matrix_set(&ff.matrix, 
			   i,
			   c, 
			   gsl_matrix_get(FtCt, i, c) + gsl_matrix_get(FtCt, c, i) + gsl_matrix_get(Q, i, c));

	}
    }

    return GSL_SUCCESS;
}



int cac_drift_in_cac_ts(int cac_drift, int o, int ts_unique, struct s_obs2ts **obs2ts)
{
    /* helper function for fourth part of eval_jac computation:
       return 1 if cac_drift is in cac of the considered time serie (o, ts_unique) */

    int n_cac_ts, c_ts, ac_ts, cac_ts;

    for(n_cac_ts=0; n_cac_ts< obs2ts[o]->n_cac[ts_unique]; n_cac_ts++) {
        c_ts = obs2ts[o]->cac[ts_unique][n_cac_ts][0];
        ac_ts = obs2ts[o]->cac[ts_unique][n_cac_ts][1];
        cac_ts = c_ts*N_AC+ac_ts;

        if(cac_ts == cac_drift) {
            return 1;
        }
    }

    return 0;
}

void eval_jac(gsl_matrix *Ft, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_group ***compo_groups_drift_par_proc, double t)
{
    //X is p_X->proj

    int c, ac, cac;
    int ts, ts_unique, stream, n_cac;

    /*t is the current time in the unit of the data */

    //syntaxic shortcut
    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;  /* syntaxic shortcut */

    

    //the automaticaly generated code may need these variables   
    double **par = p_par->natural;

    //some terms are always 0: derivative of the ODE (excluding the observed variable) against the observed variable, derivative of the dynamic of the observed variable against the observed variables, derivative of the drift eq.
    gsl_matrix_set_zero(Ft);


    double _rj[N_CAC][37];

    
    double _sf[N_CAC][1];


    for(cac=0; cac<N_CAC; cac++){
	


        
        _sf[cac][0] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR));

        
        _rj[cac][0] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);
        _rj[cac][1] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][2] = -par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][3] = -X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][4] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][5] = X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][6] = (par[ORDER_q][routers[ORDER_q]->map[cac]]);
        _rj[cac][7] = X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][8] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][9] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][10] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][11] = 0;
        _rj[cac][12] = (par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][13] = (par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][14] = X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][15] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][16] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][17] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][18] = -par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][19] = -X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][20] = (par[ORDER_v][routers[ORDER_v]->map[cac]])-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][21] = -par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][22] = -X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][23] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][24] = X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][25] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][26] = -par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][27] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-(par[ORDER_q][routers[ORDER_q]->map[cac]]);
        _rj[cac][28] = X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][29] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][30] = -X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][31] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-(par[ORDER_q][routers[ORDER_q]->map[cac]])-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][32] = -X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][33] = X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][34] = -X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][35] = (par[ORDER_v][routers[ORDER_v]->map[cac]])-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][36] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-(par[ORDER_q][routers[ORDER_q]->map[cac]])-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
    }


    //first non null part of the jacobian matrix: derivative of the ODE (excluding the observed variable) against the state variable only ( automaticaly generated code )
    for(c=0; c<N_C; c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;
            
            
            gsl_matrix_set(Ft, 0*N_CAC+cac, 0*N_CAC+cac, _rj[cac][8]);
            
            gsl_matrix_set(Ft, 0*N_CAC+cac, 1*N_CAC+cac, _rj[cac][3]);
            
            gsl_matrix_set(Ft, 0*N_CAC+cac, 2*N_CAC+cac, _rj[cac][19]);
            
            gsl_matrix_set(Ft, 0*N_CAC+cac, 3*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 0*N_CAC+cac, 4*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 0*N_CAC+cac, 5*N_CAC+cac, _rj[cac][3]);
            
            gsl_matrix_set(Ft, 0*N_CAC+cac, 6*N_CAC+cac, _rj[cac][19]);
            
            gsl_matrix_set(Ft, 0*N_CAC+cac, 7*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 0*N_CAC+cac, 8*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 0*N_CAC+cac, 9*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 0*N_CAC+cac, 10*N_CAC+cac, _rj[cac][11]);
            
            
            
            gsl_matrix_set(Ft, 1*N_CAC+cac, 0*N_CAC+cac, _rj[cac][4]);
            
            gsl_matrix_set(Ft, 1*N_CAC+cac, 1*N_CAC+cac, _rj[cac][9]);
            
            gsl_matrix_set(Ft, 1*N_CAC+cac, 2*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 1*N_CAC+cac, 3*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 1*N_CAC+cac, 4*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 1*N_CAC+cac, 5*N_CAC+cac, _rj[cac][33]);
            
            gsl_matrix_set(Ft, 1*N_CAC+cac, 6*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 1*N_CAC+cac, 7*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 1*N_CAC+cac, 8*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 1*N_CAC+cac, 9*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 1*N_CAC+cac, 10*N_CAC+cac, _rj[cac][11]);
            
            
            
            gsl_matrix_set(Ft, 2*N_CAC+cac, 0*N_CAC+cac, _rj[cac][29]);
            
            gsl_matrix_set(Ft, 2*N_CAC+cac, 1*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 2*N_CAC+cac, 2*N_CAC+cac, _rj[cac][1]);
            
            gsl_matrix_set(Ft, 2*N_CAC+cac, 3*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 2*N_CAC+cac, 4*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 2*N_CAC+cac, 5*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 2*N_CAC+cac, 6*N_CAC+cac, _rj[cac][5]);
            
            gsl_matrix_set(Ft, 2*N_CAC+cac, 7*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 2*N_CAC+cac, 8*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 2*N_CAC+cac, 9*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 2*N_CAC+cac, 10*N_CAC+cac, _rj[cac][11]);
            
            
            
            gsl_matrix_set(Ft, 3*N_CAC+cac, 0*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 3*N_CAC+cac, 1*N_CAC+cac, _rj[cac][30]);
            
            gsl_matrix_set(Ft, 3*N_CAC+cac, 2*N_CAC+cac, _rj[cac][22]);
            
            gsl_matrix_set(Ft, 3*N_CAC+cac, 3*N_CAC+cac, _rj[cac][16]);
            
            gsl_matrix_set(Ft, 3*N_CAC+cac, 4*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 3*N_CAC+cac, 5*N_CAC+cac, _rj[cac][30]);
            
            gsl_matrix_set(Ft, 3*N_CAC+cac, 6*N_CAC+cac, _rj[cac][22]);
            
            gsl_matrix_set(Ft, 3*N_CAC+cac, 7*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 3*N_CAC+cac, 8*N_CAC+cac, _rj[cac][6]);
            
            gsl_matrix_set(Ft, 3*N_CAC+cac, 9*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 3*N_CAC+cac, 10*N_CAC+cac, _rj[cac][11]);
            
            
            
            gsl_matrix_set(Ft, 4*N_CAC+cac, 0*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 4*N_CAC+cac, 1*N_CAC+cac, _rj[cac][32]);
            
            gsl_matrix_set(Ft, 4*N_CAC+cac, 2*N_CAC+cac, _rj[cac][34]);
            
            gsl_matrix_set(Ft, 4*N_CAC+cac, 3*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 4*N_CAC+cac, 4*N_CAC+cac, _rj[cac][15]);
            
            gsl_matrix_set(Ft, 4*N_CAC+cac, 5*N_CAC+cac, _rj[cac][32]);
            
            gsl_matrix_set(Ft, 4*N_CAC+cac, 6*N_CAC+cac, _rj[cac][34]);
            
            gsl_matrix_set(Ft, 4*N_CAC+cac, 7*N_CAC+cac, _rj[cac][6]);
            
            gsl_matrix_set(Ft, 4*N_CAC+cac, 8*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 4*N_CAC+cac, 9*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 4*N_CAC+cac, 10*N_CAC+cac, _rj[cac][11]);
            
            
            
            gsl_matrix_set(Ft, 5*N_CAC+cac, 0*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 5*N_CAC+cac, 1*N_CAC+cac, _rj[cac][14]);
            
            gsl_matrix_set(Ft, 5*N_CAC+cac, 2*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 5*N_CAC+cac, 3*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 5*N_CAC+cac, 4*N_CAC+cac, _rj[cac][23]);
            
            gsl_matrix_set(Ft, 5*N_CAC+cac, 5*N_CAC+cac, _rj[cac][17]);
            
            gsl_matrix_set(Ft, 5*N_CAC+cac, 6*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 5*N_CAC+cac, 7*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 5*N_CAC+cac, 8*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 5*N_CAC+cac, 9*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 5*N_CAC+cac, 10*N_CAC+cac, _rj[cac][11]);
            
            
            
            gsl_matrix_set(Ft, 6*N_CAC+cac, 0*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 6*N_CAC+cac, 1*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 6*N_CAC+cac, 2*N_CAC+cac, _rj[cac][7]);
            
            gsl_matrix_set(Ft, 6*N_CAC+cac, 3*N_CAC+cac, _rj[cac][10]);
            
            gsl_matrix_set(Ft, 6*N_CAC+cac, 4*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 6*N_CAC+cac, 5*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 6*N_CAC+cac, 6*N_CAC+cac, _rj[cac][25]);
            
            gsl_matrix_set(Ft, 6*N_CAC+cac, 7*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 6*N_CAC+cac, 8*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 6*N_CAC+cac, 9*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 6*N_CAC+cac, 10*N_CAC+cac, _rj[cac][11]);
            
            
            
            gsl_matrix_set(Ft, 7*N_CAC+cac, 0*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 7*N_CAC+cac, 1*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 7*N_CAC+cac, 2*N_CAC+cac, _rj[cac][12]);
            
            gsl_matrix_set(Ft, 7*N_CAC+cac, 3*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 7*N_CAC+cac, 4*N_CAC+cac, _rj[cac][29]);
            
            gsl_matrix_set(Ft, 7*N_CAC+cac, 5*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 7*N_CAC+cac, 6*N_CAC+cac, _rj[cac][28]);
            
            gsl_matrix_set(Ft, 7*N_CAC+cac, 7*N_CAC+cac, _rj[cac][27]);
            
            gsl_matrix_set(Ft, 7*N_CAC+cac, 8*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 7*N_CAC+cac, 9*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 7*N_CAC+cac, 10*N_CAC+cac, _rj[cac][11]);
            
            
            
            gsl_matrix_set(Ft, 8*N_CAC+cac, 0*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 8*N_CAC+cac, 1*N_CAC+cac, _rj[cac][13]);
            
            gsl_matrix_set(Ft, 8*N_CAC+cac, 2*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 8*N_CAC+cac, 3*N_CAC+cac, _rj[cac][4]);
            
            gsl_matrix_set(Ft, 8*N_CAC+cac, 4*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 8*N_CAC+cac, 5*N_CAC+cac, _rj[cac][24]);
            
            gsl_matrix_set(Ft, 8*N_CAC+cac, 6*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 8*N_CAC+cac, 7*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 8*N_CAC+cac, 8*N_CAC+cac, _rj[cac][27]);
            
            gsl_matrix_set(Ft, 8*N_CAC+cac, 9*N_CAC+cac, _rj[cac][11]);
            
            gsl_matrix_set(Ft, 8*N_CAC+cac, 10*N_CAC+cac, _rj[cac][11]);
            
            
            
            gsl_matrix_set(Ft, 9*N_CAC+cac, 0*N_CAC+cac, _rj[cac][21]);
            
            gsl_matrix_set(Ft, 9*N_CAC+cac, 1*N_CAC+cac, _rj[cac][2]);
            
            gsl_matrix_set(Ft, 9*N_CAC+cac, 2*N_CAC+cac, _rj[cac][21]);
            
            gsl_matrix_set(Ft, 9*N_CAC+cac, 3*N_CAC+cac, _rj[cac][21]);
            
            gsl_matrix_set(Ft, 9*N_CAC+cac, 4*N_CAC+cac, _rj[cac][21]);
            
            gsl_matrix_set(Ft, 9*N_CAC+cac, 5*N_CAC+cac, _rj[cac][20]);
            
            gsl_matrix_set(Ft, 9*N_CAC+cac, 6*N_CAC+cac, _rj[cac][21]);
            
            gsl_matrix_set(Ft, 9*N_CAC+cac, 7*N_CAC+cac, _rj[cac][21]);
            
            gsl_matrix_set(Ft, 9*N_CAC+cac, 8*N_CAC+cac, _rj[cac][21]);
            
            gsl_matrix_set(Ft, 9*N_CAC+cac, 9*N_CAC+cac, _rj[cac][36]);
            
            gsl_matrix_set(Ft, 9*N_CAC+cac, 10*N_CAC+cac, _rj[cac][21]);
            
            
            
            gsl_matrix_set(Ft, 10*N_CAC+cac, 0*N_CAC+cac, _rj[cac][18]);
            
            gsl_matrix_set(Ft, 10*N_CAC+cac, 1*N_CAC+cac, _rj[cac][18]);
            
            gsl_matrix_set(Ft, 10*N_CAC+cac, 2*N_CAC+cac, _rj[cac][26]);
            
            gsl_matrix_set(Ft, 10*N_CAC+cac, 3*N_CAC+cac, _rj[cac][18]);
            
            gsl_matrix_set(Ft, 10*N_CAC+cac, 4*N_CAC+cac, _rj[cac][18]);
            
            gsl_matrix_set(Ft, 10*N_CAC+cac, 5*N_CAC+cac, _rj[cac][18]);
            
            gsl_matrix_set(Ft, 10*N_CAC+cac, 6*N_CAC+cac, _rj[cac][35]);
            
            gsl_matrix_set(Ft, 10*N_CAC+cac, 7*N_CAC+cac, _rj[cac][18]);
            
            gsl_matrix_set(Ft, 10*N_CAC+cac, 8*N_CAC+cac, _rj[cac][18]);
            
            gsl_matrix_set(Ft, 10*N_CAC+cac, 9*N_CAC+cac, _rj[cac][18]);
            
            gsl_matrix_set(Ft, 10*N_CAC+cac, 10*N_CAC+cac, _rj[cac][31]);
            
            
        }
    }

    //second non null part of the jacobian matrix: derivative of the dynamic of the observed variable against the state variable only ( automaticaly generated code )
    ts = 0;
    
    for(ts_unique=0; ts_unique < obs2ts[0]->n_ts_unique; ts_unique++) {
        for(stream=0; stream < obs2ts[0]->n_stream[ts_unique]; stream++) {

            for(n_cac=0; n_cac< obs2ts[0]->n_cac[ts_unique]; n_cac++) {
                c = obs2ts[0]->cac[ts_unique][n_cac][0];
                ac = obs2ts[0]->cac[ts_unique][n_cac][1];
                cac = c*N_AC+ac;

                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 0*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 1*N_CAC+cac, _rj[cac][0]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 2*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 3*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 4*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 5*N_CAC+cac, _rj[cac][0]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 6*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 7*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 8*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 9*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 10*N_CAC+cac, _rj[cac][11]);
                

            }

            ts++;
        }
    }
    
    for(ts_unique=0; ts_unique < obs2ts[1]->n_ts_unique; ts_unique++) {
        for(stream=0; stream < obs2ts[1]->n_stream[ts_unique]; stream++) {

            for(n_cac=0; n_cac< obs2ts[1]->n_cac[ts_unique]; n_cac++) {
                c = obs2ts[1]->cac[ts_unique][n_cac][0];
                ac = obs2ts[1]->cac[ts_unique][n_cac][1];
                cac = c*N_AC+ac;

                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 0*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 1*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 2*N_CAC+cac, _rj[cac][0]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 3*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 4*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 5*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 6*N_CAC+cac, _rj[cac][0]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 7*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 8*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 9*N_CAC+cac, _rj[cac][11]);
                
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, 10*N_CAC+cac, _rj[cac][11]);
                

            }

            ts++;
        }
    }
    


    

}

/**
 * derivative of the mean of the observation process against state
 * variables and observed variables the derivative are templated
 * (automaticaly generated code) and are a function of "x", the
 * observed variable (xk_t_ts).
 */
void eval_ht(struct s_kalman_update * p, double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
    struct s_router **routers = p_data->routers;

    double **par = p_par->natural;    

    //derivative against state variable are always nul so we focus on the derivative against the observed variable
    gsl_vector_set(p->ht, N_PAR_SV*N_CAC +ts, par[ORDER_rep][routers[ORDER_rep]->map[ts]]*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])*gsl_spline_eval(p_calc->spline[ORDER_ptest][ts],t,p_calc->acc[ORDER_ptest][ts]));
}


/**
 * Approximation of the variance of a function of one random variable
 * Second order Taylor expansion
 * Var(f(X))=[f'(EX)]^2Var(X)+\frac{[f''(EX)]^2}{4}Var^2(X)
 */
double var_f_x(double varx, double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
    struct s_router **routers = p_data->routers;

    double **par = p_par->natural;   

    double derf = par[ORDER_rep][routers[ORDER_rep]->map[ts]]*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])*gsl_spline_eval(p_calc->spline[ORDER_ptest][ts],t,p_calc->acc[ORDER_ptest][ts]);
    double der2f = 0;

    return pow(derf, 2)*varx + (pow(der2f, 2)/4.0) * pow(varx, 2);
}





/**
 * evaluate Q (no_env_sto)
 */
void eval_Q_no_env_sto(gsl_matrix *Q, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t)
{

    struct s_router **routers = p_data->routers;
    //the automaticaly generated code may need these variables
    double **par = p_par->natural;    
    int i, k, cac, offset;

    

    //////////////////////////////////////////////////////////////
    // demographic stochasticity and white noise terms (if any) //
    //////////////////////////////////////////////////////////////

    
    int j;
    double term;

    

    
    double _sf[N_CAC][1];
    

    /*
      Q_proc contains only term involving state variables. We just
      replicate those term for every cac
     */

    for (cac=0; cac<N_CAC; cac++) {

	

        
        _sf[cac][0] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR));

        
        i = 0 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_SS*N_CAC+cac]+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 1 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 1 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_IS*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 2 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 2 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_SI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 3 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_RS*N_CAC+cac]+X[ORDER_QS*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 4 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_SR*N_CAC+cac]+X[ORDER_SQ*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 5 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 5 * N_CAC + cac;
        j = 5 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_IR*N_CAC+cac]+X[ORDER_IR*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 6 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 6 * N_CAC + cac;
        j = 6 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_RI*N_CAC+cac]+X[ORDER_RI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 7 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = -X[ORDER_SI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 7 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -X[ORDER_SQ*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])-X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 7 * N_CAC + cac;
        j = 7 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_SQ*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SQ*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 8 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = -X[ORDER_IS*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 8 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -X[ORDER_QS*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])-X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 8 * N_CAC + cac;
        j = 8 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_QS*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_QS*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 9 * N_CAC + cac;
        j = 5 * N_CAC + cac;
        term = -X[ORDER_IR*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 9 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_QR*N_CAC+cac]+X[ORDER_IR*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_QR*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])+par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 10 * N_CAC + cac;
        j = 6 * N_CAC + cac;
        term = -X[ORDER_RI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 10 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_RQ*N_CAC+cac]+X[ORDER_RI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_RQ*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])+par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);

        gsl_matrix_set(Q, i, j, term);

	

        
    }

    



    
    /*
      Q_obs contains only term involving at least one observed
      variable. The expansion is difficult: Q is expressed in terms of
      time series and not observed variable we only know the
      relationship between state variable and observed variable
      (that's what Q_obs gives us). We need to recreate time series
      from index of observed variable (x.to.ind and x.from.ind).

      for one given observed variable, s_obs2ts gives us the time
      series involved (and the cac involved)
    */

    int ts_unique, ts_unique_j;
    int stream, stream_j;
    int ts, ts_j;
    int cac_j;
    int n_cac, n_cac_j;

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_obs2ts *p_obs2ts, *p_obs2ts_j;

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[0];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 1 * N_CAC + cac;

		term = -X[ORDER_IS*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[0];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 5 * N_CAC + cac;

		term = -X[ORDER_IR*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[0];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 8 * N_CAC + cac;

		term = X[ORDER_IS*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[0];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 9 * N_CAC + cac;

		term = X[ORDER_IR*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains 2 observed variables: complicated case, we have to
      find the common cac in between x.i.ind and x.j.ind

      Note that we need to sum the terms is these cases...
    */


    ts=0;
    p_obs2ts = obs2ts[0];
    p_obs2ts_j = obs2ts[0];

    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {

            ts_j=0;
            for(ts_unique_j=0; ts_unique_j < p_obs2ts_j->n_ts_unique; ts_unique_j++) {
                for(stream_j=0; stream_j < p_obs2ts_j->n_stream[ts_unique_j]; stream_j++) {

                    i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                    j = N_PAR_SV*N_CAC + p_obs2ts_j->offset + ts_j;

		    term = 0.0;

                    for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                        cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

			
			term += X[ORDER_IR*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_IS*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);
			

		    }

		    gsl_matrix_set(Q, i, j, term);
		    

                    ts_j++;
                }
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[1];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 2 * N_CAC + cac;

		term = -X[ORDER_SI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[1];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 6 * N_CAC + cac;

		term = -X[ORDER_RI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[1];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 7 * N_CAC + cac;

		term = X[ORDER_SI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[1];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 10 * N_CAC + cac;

		term = X[ORDER_RI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains 2 observed variables: complicated case, we have to
      find the common cac in between x.i.ind and x.j.ind

      Note that we need to sum the terms is these cases...
    */


    ts=0;
    p_obs2ts = obs2ts[1];
    p_obs2ts_j = obs2ts[1];

    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {

            ts_j=0;
            for(ts_unique_j=0; ts_unique_j < p_obs2ts_j->n_ts_unique; ts_unique_j++) {
                for(stream_j=0; stream_j < p_obs2ts_j->n_stream[ts_unique_j]; stream_j++) {

                    i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                    j = N_PAR_SV*N_CAC + p_obs2ts_j->offset + ts_j;

		    term = 0.0;

                    for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                        cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

			
			term += X[ORDER_RI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);
			

		    }

		    gsl_matrix_set(Q, i, j, term);
		    

                    ts_j++;
                }
            }
            ts++;
        }
    }

    

    
    

    

}


/**
 * evaluate Q (full)
 */
void eval_Q_full(gsl_matrix *Q, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t)
{

    struct s_router **routers = p_data->routers;
    //the automaticaly generated code may need these variables
    double **par = p_par->natural;    
    int i, k, cac, offset;

    

    //////////////////////////////////////////////////////////////
    // demographic stochasticity and white noise terms (if any) //
    //////////////////////////////////////////////////////////////

    
    int j;
    double term;

    

    
    double _sf[N_CAC][1];
    

    /*
      Q_proc contains only term involving state variables. We just
      replicate those term for every cac
     */

    for (cac=0; cac<N_CAC; cac++) {

	

        
        _sf[cac][0] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR));

        
        i = 0 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_SS*N_CAC+cac]+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+pow(X[ORDER_SS*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2)+pow(X[ORDER_SS*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 1 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-pow(X[ORDER_SS*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 1 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_IS*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+pow(X[ORDER_SS*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 2 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-pow(X[ORDER_SS*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 2 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_SI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+pow(X[ORDER_SS*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 3 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2)+X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 3 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 3 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 3 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_RS*N_CAC+cac]+X[ORDER_QS*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+pow(X[ORDER_RS*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2)+pow(X[ORDER_RS*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 4 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2)+X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 4 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 4 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 4 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2)+X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 4 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_SR*N_CAC+cac]+X[ORDER_SQ*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+pow(X[ORDER_SR*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2)+pow(X[ORDER_SR*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 5 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 5 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 5 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 5 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-pow(X[ORDER_SR*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 5 * N_CAC + cac;
        j = 5 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_IR*N_CAC+cac]+X[ORDER_IR*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+pow(X[ORDER_SR*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 6 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 6 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 6 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-pow(X[ORDER_RS*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 6 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 6 * N_CAC + cac;
        j = 6 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_RI*N_CAC+cac]+X[ORDER_RI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+pow(X[ORDER_RS*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 7 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 7 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = -X[ORDER_SI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 7 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 7 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -X[ORDER_SQ*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])-X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-pow(X[ORDER_SR*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 7 * N_CAC + cac;
        j = 6 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 7 * N_CAC + cac;
        j = 7 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_SQ*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SQ*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+pow(X[ORDER_SR*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 8 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 8 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = -X[ORDER_IS*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 8 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -X[ORDER_QS*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])-X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-pow(X[ORDER_RS*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 8 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 8 * N_CAC + cac;
        j = 5 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 8 * N_CAC + cac;
        j = 8 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_QS*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_QS*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+pow(X[ORDER_RS*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 9 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 5 * N_CAC + cac;
        term = -X[ORDER_IR*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 8 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 9 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_QR*N_CAC+cac]+X[ORDER_IR*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_QR*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])+par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*pow(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 10 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 6 * N_CAC + cac;
        term = -X[ORDER_RI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_RS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 7 * N_CAC + cac;
        term = X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 10 * N_CAC + cac;
        term = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])*X[ORDER_RQ*N_CAC+cac]+X[ORDER_RI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_RQ*N_CAC+cac]*(par[ORDER_q][routers[ORDER_q]->map[cac]])+par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*pow(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
    }

    



    
    /*
      Q_obs contains only term involving at least one observed
      variable. The expansion is difficult: Q is expressed in terms of
      time series and not observed variable we only know the
      relationship between state variable and observed variable
      (that's what Q_obs gives us). We need to recreate time series
      from index of observed variable (x.to.ind and x.from.ind).

      for one given observed variable, s_obs2ts gives us the time
      series involved (and the cac involved)
    */

    int ts_unique, ts_unique_j;
    int stream, stream_j;
    int ts, ts_j;
    int cac_j;
    int n_cac, n_cac_j;

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_obs2ts *p_obs2ts, *p_obs2ts_j;

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[0];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 1 * N_CAC + cac;

		term = -X[ORDER_IS*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[0];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 5 * N_CAC + cac;

		term = -X[ORDER_IR*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[0];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 8 * N_CAC + cac;

		term = X[ORDER_IS*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[0];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 9 * N_CAC + cac;

		term = X[ORDER_IR*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains 2 observed variables: complicated case, we have to
      find the common cac in between x.i.ind and x.j.ind

      Note that we need to sum the terms is these cases...
    */


    ts=0;
    p_obs2ts = obs2ts[0];
    p_obs2ts_j = obs2ts[0];

    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {

            ts_j=0;
            for(ts_unique_j=0; ts_unique_j < p_obs2ts_j->n_ts_unique; ts_unique_j++) {
                for(stream_j=0; stream_j < p_obs2ts_j->n_stream[ts_unique_j]; stream_j++) {

                    i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                    j = N_PAR_SV*N_CAC + p_obs2ts_j->offset + ts_j;

		    term = 0.0;

                    for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                        cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

			
			term += X[ORDER_IR*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_IS*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);
			

		    }

		    gsl_matrix_set(Q, i, j, term);
		    

                    ts_j++;
                }
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[1];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 2 * N_CAC + cac;

		term = -X[ORDER_SI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[1];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 6 * N_CAC + cac;

		term = -X[ORDER_RI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[1];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 7 * N_CAC + cac;

		term = X[ORDER_SI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[1];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                j = 10 * N_CAC + cac;

		term = X[ORDER_RI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);

                gsl_matrix_set(Q, i, j, term);
		gsl_matrix_set(Q, j, i, term);
            }
            ts++;
        }
    }

    

    

    

    /*
      x contains 2 observed variables: complicated case, we have to
      find the common cac in between x.i.ind and x.j.ind

      Note that we need to sum the terms is these cases...
    */


    ts=0;
    p_obs2ts = obs2ts[1];
    p_obs2ts_j = obs2ts[1];

    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {

            ts_j=0;
            for(ts_unique_j=0; ts_unique_j < p_obs2ts_j->n_ts_unique; ts_unique_j++) {
                for(stream_j=0; stream_j < p_obs2ts_j->n_stream[ts_unique_j]; stream_j++) {

                    i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                    j = N_PAR_SV*N_CAC + p_obs2ts_j->offset + ts_j;

		    term = 0.0;

                    for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                        cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

			
			term += X[ORDER_RI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SI*N_CAC+cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]);
			

		    }

		    gsl_matrix_set(Q, i, j, term);
		    

                    ts_j++;
                }
            }
            ts++;
        }
    }

    

    
    

    

}


/**
 * evaluate Q (no_dem_sto)
 */
void eval_Q_no_dem_sto(gsl_matrix *Q, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t)
{

    struct s_router **routers = p_data->routers;
    //the automaticaly generated code may need these variables
    double **par = p_par->natural;    
    int i, k, cac, offset;

    

    //////////////////////////////////////////////////////////////
    // demographic stochasticity and white noise terms (if any) //
    //////////////////////////////////////////////////////////////

    
    int j;
    double term;

    

    
    double _sf[N_CAC][1];
    

    /*
      Q_proc contains only term involving state variables. We just
      replicate those term for every cac
     */

    for (cac=0; cac<N_CAC; cac++) {

	

        
        _sf[cac][0] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR));

        
        i = 0 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = pow(X[ORDER_SS*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2)+pow(X[ORDER_SS*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 1 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -pow(X[ORDER_SS*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 1 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = pow(X[ORDER_SS*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 2 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -pow(X[ORDER_SS*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 2 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = pow(X[ORDER_SS*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 3 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2)+X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 3 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 3 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 3 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = pow(X[ORDER_RS*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2)+pow(X[ORDER_RS*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 4 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2)+X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 4 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 4 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 4 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2)+X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 4 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = pow(X[ORDER_SR*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2)+pow(X[ORDER_SR*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 5 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 5 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 5 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 5 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -pow(X[ORDER_SR*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 5 * N_CAC + cac;
        j = 5 * N_CAC + cac;
        term = pow(X[ORDER_SR*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 6 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 6 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 6 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -pow(X[ORDER_RS*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 6 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 6 * N_CAC + cac;
        j = 6 * N_CAC + cac;
        term = pow(X[ORDER_RS*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 7 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 7 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = X[ORDER_SR*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 7 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 7 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -pow(X[ORDER_SR*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 7 * N_CAC + cac;
        j = 6 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 7 * N_CAC + cac;
        j = 7 * N_CAC + cac;
        term = pow(X[ORDER_SR*N_CAC+cac],2)*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 8 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 8 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 8 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -pow(X[ORDER_RS*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 8 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 8 * N_CAC + cac;
        j = 5 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 8 * N_CAC + cac;
        j = 8 * N_CAC + cac;
        term = pow(X[ORDER_RS*N_CAC+cac],2)*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 9 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 1 * N_CAC + cac;
        term = X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 5 * N_CAC + cac;
        term = X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 8 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 9 * N_CAC + cac;
        j = 9 * N_CAC + cac;
        term = pow(par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]],2)*pow(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
        i = 10 * N_CAC + cac;
        j = 0 * N_CAC + cac;
        term = -X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 2 * N_CAC + cac;
        term = X[ORDER_SS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 3 * N_CAC + cac;
        term = -X[ORDER_RS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 4 * N_CAC + cac;
        term = -X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 6 * N_CAC + cac;
        term = X[ORDER_RS*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 7 * N_CAC + cac;
        term = X[ORDER_SR*N_CAC+cac]*pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	
        gsl_matrix_set(Q, j, i, term);
	

        
        i = 10 * N_CAC + cac;
        j = 10 * N_CAC + cac;
        term = pow(par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]],2)*pow(par[ORDER_sigma][routers[ORDER_sigma]->map[cac]],2)*pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]],2)*pow(par[ORDER_v][routers[ORDER_v]->map[cac]],2)*pow(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0,2)*pow(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]],2)*pow(-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]+gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac],2)/pow(gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]),2);

        gsl_matrix_set(Q, i, j, term);

	

        
    }

    



    

    

}


/**
 * evaluate Q (no_dem_sto_no_env_sto)
 */
void eval_Q_no_dem_sto_no_env_sto(gsl_matrix *Q, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t)
{

    struct s_router **routers = p_data->routers;
    //the automaticaly generated code may need these variables
    double **par = p_par->natural;    
    int i, k, cac, offset;

    

    //////////////////////////////////////////////////////////////
    // demographic stochasticity and white noise terms (if any) //
    //////////////////////////////////////////////////////////////

    



    

    

}


