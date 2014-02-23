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

#include "simulation.h"

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
#define ORDER_p 22
#define ORDER_rep 23


#define ORDER_U 11
#define ORDER_DU 12




#define ORDER_ptest 0
#define ORDER_mu_b 1
#define ORDER_mu_d 2

void ensure_cst_pop_size(struct s_data *p_data)
{
    
    int nn, cac; 
    print_warning("variable birth and death rate (mu_b and mu_d) detected in covariates. mu_d have been set to mu_b to ensure a constant population size to analyze the attractor");

    for (nn=0; nn < N_DATA_PAR_FIXED; nn++) {
        for (cac=0; cac < N_CAC; cac++) {
            p_data->par_fixed[ORDER_mu_d][nn][cac] = p_data->par_fixed[ORDER_mu_b][nn][cac];
        }
    }
    
}

int step_lyap (double t, const double X[], double f[], void *params)
{
    struct s_calc *p_calc = (struct s_calc *) params;
    struct s_par *p_par = p_calc->p_par;  /* syntaxic shortcut */
    struct s_data *p_data = p_calc->p_data;

    int i;
    int c, ac, cac;
    const int nn = p_calc->current_nn;
    double **par = p_par->natural;
    double ***covar = p_data->par_fixed;

    struct s_router **routers = p_data->routers;  /* syntaxic shortcut */

    /* non linear system (automatically generated code)*/
    double _r[N_CAC][17];

    
    double _sf[N_CAC][1];

    

    for(cac=0;cac<N_CAC;cac++) {
	

        
        _sf[cac][0] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR_IN_DATA_UNIT));

        
        _r[cac][0] = (par[ORDER_q][routers[ORDER_q]->map[cac]]);
        _r[cac][1] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _r[cac][2] = covar[ORDER_mu_b][nn][cac]*p_data->pop_size_t0[cac]*par[ORDER_p][routers[ORDER_p]->map[cac]];
        _r[cac][3] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _r[cac][4] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);
        _r[cac][5] = covar[ORDER_mu_b][nn][cac]*p_data->pop_size_t0[cac]*(-par[ORDER_p][routers[ORDER_p]->map[cac]]+1.0);
        _r[cac][6] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac];
        _r[cac][7] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(p_data->pop_size_t0[cac]-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/p_data->pop_size_t0[cac];
        _r[cac][8] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac];
        _r[cac][9] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _r[cac][10] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])*(p_data->pop_size_t0[cac]-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/p_data->pop_size_t0[cac];
        _r[cac][11] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(p_data->pop_size_t0[cac]-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/p_data->pop_size_t0[cac];
        _r[cac][12] = covar[ORDER_mu_d][nn][cac];
        _r[cac][13] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _r[cac][14] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac];
        _r[cac][15] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac];
        _r[cac][16] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])*(p_data->pop_size_t0[cac]-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/p_data->pop_size_t0[cac];
    }

    for (c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

	    
	    f[0*N_CAC+cac] =  - (_r[cac][12]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][6]*X[ORDER_SS*N_CAC+cac]) - (_r[cac][1]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][5]);
	    f[1*N_CAC+cac] =  - (_r[cac][12]*X[ORDER_IS*N_CAC+cac]) - (_r[cac][4]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][6]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][14]*X[ORDER_RS*N_CAC+cac]);
	    f[2*N_CAC+cac] =  - (_r[cac][12]*X[ORDER_SI*N_CAC+cac]) - (_r[cac][4]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][1]*X[ORDER_SS*N_CAC+cac]) + (_r[cac][9]*X[ORDER_SR*N_CAC+cac]);
	    f[3*N_CAC+cac] =  - (_r[cac][12]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][13]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][14]*X[ORDER_RS*N_CAC+cac]) - (_r[cac][8]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][2]) + (_r[cac][0]*X[ORDER_QS*N_CAC+cac]);
	    f[4*N_CAC+cac] =  - (_r[cac][12]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][15]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][9]*X[ORDER_SR*N_CAC+cac]) - (_r[cac][3]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]);
	    f[5*N_CAC+cac] =  - (_r[cac][12]*X[ORDER_IR*N_CAC+cac]) - (_r[cac][4]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][15]*X[ORDER_SR*N_CAC+cac]) + (_r[cac][16]);
	    f[6*N_CAC+cac] =  - (_r[cac][12]*X[ORDER_RI*N_CAC+cac]) - (_r[cac][4]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][13]*X[ORDER_RS*N_CAC+cac]) + (_r[cac][10]);
	    f[7*N_CAC+cac] =  - (_r[cac][12]*X[ORDER_SQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_SQ*N_CAC+cac]) + (_r[cac][4]*X[ORDER_SI*N_CAC+cac]) + (_r[cac][3]*X[ORDER_SR*N_CAC+cac]);
	    f[8*N_CAC+cac] =  - (_r[cac][12]*X[ORDER_QS*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QS*N_CAC+cac]) + (_r[cac][4]*X[ORDER_IS*N_CAC+cac]) + (_r[cac][8]*X[ORDER_RS*N_CAC+cac]);
	    f[9*N_CAC+cac] =  - (_r[cac][12]*X[ORDER_QR*N_CAC+cac]) - (_r[cac][0]*X[ORDER_QR*N_CAC+cac]) + (_r[cac][4]*X[ORDER_IR*N_CAC+cac]) + (_r[cac][11]);
	    f[10*N_CAC+cac] =  - (_r[cac][12]*X[ORDER_RQ*N_CAC+cac]) - (_r[cac][0]*X[ORDER_RQ*N_CAC+cac]) + (_r[cac][4]*X[ORDER_RI*N_CAC+cac]) + (_r[cac][7]);
        }
    }

    /* linear system: product of jacobian matrix (DIM*DIM) per

       | y[1*DIM+0]       y[1*DIM+1] ...     y[1*DIM+(DIM-1)]   |
       | y[2*DIM+0]       y[2*DIM+1] ...     y[2*DIM+(DIM-1)]   |
       | ...                                                    |
       | y[DIM*DIM+0]     y[DIM*DIM+1] ...  y[DIM*DIM+(DIM-1)]  |

       (automaticaly generated code)
    */

    double _rj[N_CAC][42];
    for(cac=0; cac<N_CAC; cac++){
        
        _rj[cac][0] = X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac]+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][1] = -covar[ORDER_mu_d][nn][cac]-(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac]+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][2] = -par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac]+par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(p_data->pop_size_t0[cac]-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/p_data->pop_size_t0[cac];
        _rj[cac][3] = -covar[ORDER_mu_d][nn][cac]-(par[ORDER_q][routers[ORDER_q]->map[cac]])-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][4] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][5] = (par[ORDER_q][routers[ORDER_q]->map[cac]]);
        _rj[cac][6] = X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac]+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][7] = -par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][8] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][9] = -par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][10] = -X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac]-X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][11] = -par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][12] = -par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][13] = -covar[ORDER_mu_d][nn][cac]-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac]-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac]-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][14] = -covar[ORDER_mu_d][nn][cac]-(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac]-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac]+par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(p_data->pop_size_t0[cac]-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/p_data->pop_size_t0[cac];
        _rj[cac][15] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][16] = (par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][17] = -X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac]-X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][18] = -X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][19] = 0;
        _rj[cac][20] = (par[ORDER_v][routers[ORDER_v]->map[cac]])-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac]+par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(p_data->pop_size_t0[cac]-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/p_data->pop_size_t0[cac];
        _rj[cac][21] = -par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac]+par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][22] = -par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac]+par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(p_data->pop_size_t0[cac]-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/p_data->pop_size_t0[cac];
        _rj[cac][23] = -covar[ORDER_mu_d][nn][cac]-(par[ORDER_q][routers[ORDER_q]->map[cac]]);
        _rj[cac][24] = X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac]-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac]+par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(p_data->pop_size_t0[cac]-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/p_data->pop_size_t0[cac];
        _rj[cac][25] = X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][26] = -X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][27] = -covar[ORDER_mu_d][nn][cac]-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac]-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][28] = -covar[ORDER_mu_d][nn][cac]-(par[ORDER_q][routers[ORDER_q]->map[cac]])-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][29] = X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][30] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][31] = par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][32] = (par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][33] = -covar[ORDER_mu_d][nn][cac]-(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac]+X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][34] = -covar[ORDER_mu_d][nn][cac]-(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_RS*N_CAC+cac]*par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac]-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac]+par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(p_data->pop_size_t0[cac]-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/p_data->pop_size_t0[cac];
        _rj[cac][35] = par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][36] = X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac]-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac]+par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(p_data->pop_size_t0[cac]-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/p_data->pop_size_t0[cac];
        _rj[cac][37] = -par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac]+par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
        _rj[cac][38] = (par[ORDER_v][routers[ORDER_v]->map[cac]])-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac]+par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(p_data->pop_size_t0[cac]-X[ORDER_IR*N_CAC+cac]-X[ORDER_IS*N_CAC+cac]-X[ORDER_QR*N_CAC+cac]-X[ORDER_QS*N_CAC+cac]-X[ORDER_RI*N_CAC+cac]-X[ORDER_RQ*N_CAC+cac]-X[ORDER_RS*N_CAC+cac]-X[ORDER_SI*N_CAC+cac]-X[ORDER_SQ*N_CAC+cac]-X[ORDER_SR*N_CAC+cac]-X[ORDER_SS*N_CAC+cac])/p_data->pop_size_t0[cac];
        _rj[cac][39] = -X[ORDER_SR*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][40] = -X[ORDER_SS*N_CAC+cac]*par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)/p_data->pop_size_t0[cac];
        _rj[cac][41] = -covar[ORDER_mu_d][nn][cac]-par[ORDER_r0_1][routers[ORDER_r0_1]->map[cac]]*par[ORDER_sigma][routers[ORDER_sigma]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_IR*N_CAC+cac]+X[ORDER_IS*N_CAC+cac]+par[ORDER_iota_1][routers[ORDER_iota_1]->map[cac]])/p_data->pop_size_t0[cac]-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*par[ORDER_z][routers[ORDER_z]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac]-par[ORDER_r0_2][routers[ORDER_r0_2]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(-par[ORDER_z][routers[ORDER_z]->map[cac]]+1.0)*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1.0)*(X[ORDER_RI*N_CAC+cac]+X[ORDER_SI*N_CAC+cac]+par[ORDER_iota_2][routers[ORDER_iota_2]->map[cac]])/p_data->pop_size_t0[cac];
    }


    for(c=0; c<N_C; c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;
            for(i=0; i<(N_PAR_SV*N_CAC); i++) {
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (0*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][27])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][40])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][26])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][40])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][26])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (1*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][8])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][33])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][35])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][6])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (2*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][15])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][1])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][30])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][0])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (3*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][17])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][18])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][17])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][18])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][5])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (4*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][39])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][10])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][41])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][39])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][10])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][5])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (5*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][12])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][36])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][12])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][12])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][21])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][14])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][12])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][12])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][12])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][12])*X[N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][12])*X[N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (6*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][11])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][11])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][24])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][37])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][11])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][11])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][34])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][11])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][11])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][11])*X[N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][11])*X[N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (7*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][16])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][31])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][29])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][23])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (8*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][32])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][4])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][25])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][23])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (9*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][9])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][22])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][9])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][9])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][9])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][38])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][9])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][9])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][9])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][28])*X[N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][9])*X[N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (10*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][7])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][7])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][2])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][7])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][7])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][7])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][20])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][7])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][7])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][7])*X[N_PAR_SV*N_CAC+ (9*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][3])*X[N_PAR_SV*N_CAC+ (10*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
            }
        }
    }

    return GSL_SUCCESS;
}
