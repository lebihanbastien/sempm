#include "poincare.h"
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>

#define MAX_REF_ITER 500	//constant for refinement iterations in Poincare maps

/**
 * \file poincare.cpp
 * \brief Computation of Poincaré maps.
 */

//----------------------------------------------------------------------------------------
//
//  Poincaré map
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Computes a Poincaré map
 *   \param pmap a reference to the Poincaré map parameters
 *
 *    Requires init_inv_man and init_coc
 **/
void pmap_build_random(Pmap& pmap)
{
	cout << "---------------------------------------------------" << endl;
	cout << "                                                   " << endl;
	cout << "               Poincare map computation            " << endl;
	cout << "                                                   " << endl;
	cout << "---------------------------------------------------" << endl;
	cout << std::showpos << setiosflags(ios::scientific) << setprecision(15);

    //------------------------------------------------------------------------------------
    //Filename to save outputs
    //------------------------------------------------------------------------------------
	string filename = pmap_filename(pmap);


	//------------------------------------------------------------------------------------
	//Plot variables
	//------------------------------------------------------------------------------------
	int isPlot = 1;
	gnuplot_ctrl* h1;
	h1 = gnuplot_init();
	char ch;
	int color;
	if (isPlot)
		color = 1;

    cout << "pmap_build_random. plot saved in " << filename << ".eps" << endl;


	//------------------------------------------------------------------------------------
	//Energy of the libration point at t = t0
	//------------------------------------------------------------------------------------
	double st0[4];
	for (int i = 0; i < 4; i++)
		st0[i] = 0.0;
	pmap.h_nc_li = h_nc_pmap(pmap, st0);

	//------------------------------------------------------------------------------------
	//Set the Hamiltonian h_nc to h_nc_li + dh_nc_t0
	//------------------------------------------------------------------------------------
	pmap.h_nc = pmap.h_nc_li + pmap.dh_nc_t0;

	//------------------------------------------------------------------------------------
	//Print
	//------------------------------------------------------------------------------------
	//The right header is written in the txt file
	cout << "pmap_build_random. data saved in " << filename << ".txt" << endl;
	header_fprint_bin(filename);

	//------------------------------------------------------------------------------------
	//Generating the quasi-random sequence of initial condition (s1, s3)
	//Inside a grid of the form  [pmap.si_min pmap.si_max] x [pmap.si_min pmap.si_max]
	//------------------------------------------------------------------------------------
	//Building necessary objects
	gsl_qrng* q = gsl_qrng_alloc(gsl_qrng_sobol, 4);
	double** randseq = dmatrix(0, pmap.n_sol - 1, 0, 3);

	//Storing sequence, all variables in [0, 1]
	for (int i = 0; i < pmap.n_sol; i++)
		gsl_qrng_get(q, randseq[i]);

	//Normalisation, all variables in [pmap.si_min pmap.si_max]
	for (int i = 0; i < pmap.n_sol; i++)
	{
		//s1 & s3
		randseq[i][0] = pmap.si_min * (1 - randseq[i][0])
				+ pmap.si_max * randseq[i][0];
		randseq[i][2] = pmap.si_min * (1 - randseq[i][2])
				+ pmap.si_max * randseq[i][2];
	}

	//------------------------------------------------------------------------------------
	//Loop
	//------------------------------------------------------------------------------------
	int label = 1;
	cout << std::noshowpos << resetiosflags(ios::scientific) << setprecision(4);
    #pragma omp parallel for if(pmap.is_par) shared(label)
	for (int ii = 0; ii < pmap.n_sol; ii++)
	{
		//--------------------------------------------------------------------------------
		//Integration tools
		//--------------------------------------------------------------------------------
		OdeStruct ode_s_6;
		pmap_init_ode(ode_s_6, pmap);
		//For root finding
		OdeStruct ode_s_6_root;
		pmap_init_ode(ode_s_6_root, pmap);

		//--------------------------------------------------------------------------------
		//Event structures
		//--------------------------------------------------------------------------------
		value_params val_par;
		val_par.dim = 2;                   //event on z component
		val_par.direction = 0;                   //all zeros are detected
		val_par.max_events = pmap.max_events;     //maximum of events
		val_par.value = 0.0;                 //event when z = 0
		value_function fvalue;
		fvalue.val_par = &val_par;
		fvalue.value = &linear_intersection;

		//--------------------------------------------------------------------------------
		//Orbit structure
		//--------------------------------------------------------------------------------
		double time;
		int status;
		Ofsc orbit_ofs(OFS_ORDER);
		Orbit orbit;
		init_orbit(&orbit, &CM, &CMh, &JCM, &Mcoc, &Vcoc, &val_par, &fvalue,
				&ode_s_6, &ode_s_6_root, &pmap, &SEML, &orbit_ofs,
				pmap.var_dim_h, label);

		//--------------------------------------------------------------------------------
		// Move on the grid of initial conditions
		//--------------------------------------------------------------------------------
		double* sti = dvector(0, 3);
		sti[0] = randseq[ii][0];
		sti[1] = 0.0;
		sti[2] = randseq[ii][2];
		sti[3] = 0.0;

		//--------------------------------------------------------------------------------
		//Update the IC so that the initial energy is equal to pmap.dh_nc_t0
		//--------------------------------------------------------------------------------
		status = orbit_init_pmap(orbit, sti);

		//Between the lines: if we do not want to correct the energy
		//------------------------------------------------------------------------------//
		//status = GSL_SUCCESS;
		//orbit_update_ic(orbit, sti);
		//cout << "Energy = " << h_nc_pmap(pmap, sti) << endl;
		//------------------------------------------------------------------------------//

		//--------------------------------------------------------------------------------
		//If the correction is succesful, we can go on with the computation of the return map
		//--------------------------------------------------------------------------------
		if (status == GSL_SUCCESS)
		{
			//----------------------------------------------------------------------------
			// Computation
			//----------------------------------------------------------------------------
			tic();
			status = orbit_compute_pmap(&orbit, pmap.int_method);
			time = toc();

			//----------------------------------------------------------------------------
			// Postprocess
			//----------------------------------------------------------------------------
			if (status == GSL_SUCCESS)
			{
				cout << "Return map " << ii+1 << "/" << pmap.n_sol
						<< " computed in " << time << "s. " << endl;

				//------------------------------------------------------------------------
				// Print in file
				//------------------------------------------------------------------------
                #pragma omp critical
				{
					orbit.label = label++;
					orbit.tf = pmap.tmax_nc;
					pmap_orbit_fprint_bin(&orbit, filename, 1);
				}

				//------------------------------------------------------------------------
				// Plot
				//------------------------------------------------------------------------
                #pragma omp critical
				{
					orbit_poincare_plot(&orbit, h1, color++);
				}

			}
		}else
		{
            cout << "Return map " << ii+1 << "/" << pmap.n_sol
						<< " is not computed since no IC with the desired energy have been found" << endl;


		}

		//Memory release
		free_orbit(&orbit);
		free_dvector(sti, 0, 3);
	}

	//Memory release
	gsl_qrng_free(q);
	free_dmatrix(randseq, 0, pmap.n_sol - 1, 0, 3);

	//Plot handling at the end
	if (isPlot)
	{
		printf("Press ENTER to close the gnuplot window(s)\n");
		scanf("%c", &ch);

        string set_output = "set output \""+ filename+".eps\"";
		//Save in EPS format
		gnuplot_cmd(h1, "set terminal postscript eps solid color enhanced");
		gnuplot_cmd(h1, (char*) set_output.c_str());
		gnuplot_cmd(h1, "replot");
	}
	gnuplot_close(h1);

}

//----------------------------------------------------------------------------------------
//
//  Init ode routines
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Initialize one ode structure ode_s_6 with parameters of the Poincaré map pmap.
 **/
void pmap_init_ode(OdeStruct& ode_s_6, Pmap& pmap)
{
	//Root-finding
	const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
	//Stepper
	const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
	//Init ode structure
	init_ode_structure(&ode_s_6, T,            //stepper: int
			T_root,       //stepper: root
			pmap.pabs,    //precision: int_abs
			pmap.prel,    //precision: int_rel
			pmap.proot,   //precision: root
			6,            //dimension
			Config::configManager().G_PREC_HSTART(), //initial int step
			qbcp_vfn_novar, //vector field
			NULL,    //jacobian
			&SEML);  //parameters
}

//----------------------------------------------------------------------------------------
//
//  Print & Read
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Returns the output filename associated with the current pmap
 **/
string pmap_filename(Pmap& pmap)
{
	string ss_h_nc = num_to_string(pmap.dh_nc_t0);
	string ss_order = num_to_string(pmap.ofts_order);
	string ss_ofs_order = num_to_string(pmap.ofs_order);
	string ss_max_event = num_to_string(pmap.max_events);
	string ss_type, ss_method;

	//Get the type
	ss_type = pmap.t0_nc == 0 ? "Serv_pm_t0_" : "Serv_pm_";

	//Get the method
	switch (pmap.int_method)
	{
	case DUAL_INT:
		ss_method = "_DUAL_INT";
		break;

	case DUAL_INT_NO_RESET:
		ss_method = "_DUAL_INT_NOT_RESET";
		break;

	case DUAL_INT_STEPPED:
		ss_method = "_DUAL_INT_STEPPED";
		break;

	case SINGLE_INT:
		ss_method = "_SINGLE_INT";
		break;
	}

	//Frequency of projection per period
	string ssProj = num_to_string(pmap.f_proj_T);

	//Get the PM style
	string filename;
	switch (SEML.param_style)
	{
	case Csts::GRAPH:
		filename = SEML.cs.F_PRINT + ss_type; //default case, so no additionnal notations
		break;
	case Csts::NORMFORM:
		filename = SEML.cs.F_PRINT + ss_type + "NF_"; //Normal form style
		break;
	case Csts::MIXED:
		filename = SEML.cs.F_PRINT + ss_type + "MX_";  //Normal form style
		break;
	}

	//Final name
	filename = filename + "Energy_" + ss_h_nc + "_order_" + ss_order + "_ofs_"
			+ ss_ofs_order + "_proj_" + ssProj + "_max_events_" + ss_max_event
			+ ss_method;

	//Return
	return filename;
}

/**
 *   \brief Reset a given binary file
 **/
void header_fprint_bin(string filename)
{
	fstream myfile;
	myfile.open((filename + ".bin").c_str(), std::ios::binary | ios::out);
	myfile.close();
}

/**
 *   \brief Print the poincare map of an orbit in a bin file.
 *          Only the points for which pz > 0 are kept.
 **/
void pmap_orbit_fprint_bin(Orbit* orbit, string filename, int append)
{
	if (orbit->int_method == -1)
	{
		cout << "orbit_print. The orbit was not previously computed. No print."
				<< endl;
	}
	else
	{
		//ifstream readStream;
		fstream myfile;
		string ss1;
		double zEM[6], hz;

		//Open stream
		if (append)
			myfile.open((filename + ".bin").c_str(),
					std::ios::app | std::ios::binary | ios::out); //for appending at the end of the file
		else
			myfile.open((filename + ".bin").c_str(),
					std::ios::binary | ios::out);

		//First line
		ncsys_m_to_sys_m(orbit->pmap->t0_nc, orbit->z0, zEM,
				(FBPL*) orbit->ode_s_6->d->sys->params);
		hz = qbcp_H(orbit->pmap->t0_nc, zEM, orbit->ode_s_6->d->sys->params);

		double res;

		//1. label
		res = orbit->label;
		myfile.write((char*) &res, sizeof(double));

		//2. NC
		for (int k = 0; k < 6; k++)
		{
			res = creal(orbit->z0[k]);
			myfile.write((char*) &res, sizeof(double));
		}
		//3. RCM
		for (int k = 0; k < 4; k++)
		{
			res = creal(orbit->si[k]);
			myfile.write((char*) &res, sizeof(double));
		}

		//4. time
		res = creal(orbit->pmap->t0_nc);
		myfile.write((char*) &res, sizeof(double));

		//5. dHz
		res = creal(hz - orbit->pmap->h_nc_li);
		myfile.write((char*) &res, sizeof(double));

		//6. dHw
		res = creal(hz - orbit->pmap->h_nc_li);
		myfile.write((char*) &res, sizeof(double));

		//7. dHz - dHw
		res = creal(0.0);
		myfile.write((char*) &res, sizeof(double));

		//8. pmap.dh_nc_t0
		res = creal(orbit->pmap->dh_nc_t0);
		myfile.write((char*) &res, sizeof(double));

		//9.ePm
		res = creal(0.0);
		myfile.write((char*) &res, sizeof(double));

		//10. number
		res = creal(0.0);
		myfile.write((char*) &res, sizeof(double));

		//Loop on all events
		for (int i = 0; i <= orbit->last_indix; i++)
		{
			//if pz > 0, we save the point
			if (orbit->z0_mat[5][i] > 0)
			{
				//1. label
				res = orbit->label;
				myfile.write((char*) &res, sizeof(double));
				//2. NC
				for (int k = 0; k < 6; k++)
				{
					res = orbit->z0_mat[k][i];
					myfile.write((char*) &res, sizeof(double));
				}
				//3. RCM
				for (int k = 0; k < 4; k++)
				{
					res = orbit->s0_mat[k][i];
					myfile.write((char*) &res, sizeof(double));
				}

				//4. time
				res = orbit->te_mat[i];
				myfile.write((char*) &res, sizeof(double));

				//5. dHz
				res = orbit->hz[i] - orbit->pmap->h_nc_li;
				myfile.write((char*) &res, sizeof(double));

				//6. dHw
				res = orbit->hw[i] - orbit->pmap->h_nc_li;
				myfile.write((char*) &res, sizeof(double));

				//7. dHz - dHw
				res = orbit->hz[i] - orbit->hw[i];
				myfile.write((char*) &res, sizeof(double));

				//8. pmap.dh_nc_t0
				res = orbit->pmap->dh_nc_t0;
				myfile.write((char*) &res, sizeof(double));

				//9.ePm
				res = orbit->ePm[i];
				myfile.write((char*) &res, sizeof(double));

				//10. number
				res = orbit->nevent[i];
				myfile.write((char*) &res, sizeof(double));
			}
		}

		myfile.close();
	}
}

//----------------------------------------------------------------------------------------
//
//  Init of one orbit
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Update the reduced coordinates st0 with values from orbit.s0 and the value sr.
 *          Used in orbit_init_pmap and delta_h_nc_orbit.
 *
 *          The following patterns are followed (examples):
 *              - if SEML.model = CRTBP and orbit.dim = 2: s2 = sr, s4 = 0
 *                   st0[0] = orbit->si[0];;
 *                   st0[1] = sr;
 *                   st0[2] = orbit->si[2];
 *                   st0[3] = 0.0;
 *
 *             - if SEML.model = QBCP and orbit.dim = 2: s2 = s4 = sr
 *                   st0[0] = sr;
 *                   st0[1] = orbit->si[1];
 *                   st0[2] = sr;
 *                   st0[3] = orbit->si[3];
 **/
void orbit_update_s0(Orbit* orbit, double st0[], double sr)
{
	//Set the current state
	if (SEML.model == Csts::CRTBP) //if model = CRTBP
	{
		switch (orbit->var_dim_h)
		{
		case 1:
		{
			st0[0] = sr;
			st0[1] = orbit->si[1];
			st0[2] = orbit->si[2];
			st0[3] = 0.0;

			break;
		}

		case 2:
		{
			st0[0] = orbit->si[0];
			st0[1] = sr;
			st0[2] = orbit->si[2];
			st0[3] = 0;

			break;
		}

		case 3:
		{
			st0[0] = orbit->si[0];
			st0[1] = orbit->si[1];
			st0[2] = sr;
			st0[3] = 0.0;

			break;
		}

		default:  //as 3
		{
			st0[0] = orbit->si[0];
			st0[1] = orbit->si[1];
			st0[2] = sr;
			st0[3] = 0.0;
			break;
		}
		}
	}
	else //if model = QBCP
	{
		//The starting condition z(t0) = 0 is s2 == s4
		switch (orbit->var_dim_h)
		{
		case 1:
		{
			st0[0] = sr;
			st0[1] = orbit->si[1];
			st0[2] = sr;
			st0[3] = orbit->si[3];
			break;
		}
		case 2:
		{
			st0[0] = orbit->si[0];
			st0[1] = sr;
			st0[2] = orbit->si[2];
			st0[3] = sr;
			break;
		}
		case 3:
		{
			st0[0] = sr;
			st0[1] = orbit->si[1];
			st0[2] = sr;
			st0[3] = orbit->si[3];
			break;
		}
		case 4:
		{
			st0[0] = orbit->si[0];
			st0[1] = sr;
			st0[2] = orbit->si[2];
			st0[3] = sr;
			break;
		}
		default: //as case 3
		{
			st0[0] = sr;
			st0[1] = orbit->si[1];
			st0[2] = sr;
			st0[3] = orbit->si[3];
			break;
		}
		}
	}
}

/**
 *   \brief Initialize an orbit with respect to a Poincaré map so that:
 *          H(orbit.s0, t0) - H(0, t0) = pmap.dh_nc_t0
 **/
int orbit_init_pmap(Orbit& orbit, double st0[])
{
	//------------------------------------------------------------------------------------
	//Starting by updating the IC
	//------------------------------------------------------------------------------------
	orbit_update_ic(orbit, st0, orbit.pmap->t0_nc);

	//------------------------------------------------------------------------------------
	//Root finding
	//------------------------------------------------------------------------------------
	gsl_function F;                //called function in the root finding routine
	double s_low, s_high;            //variables for root bracketing
	double Hup, Hdown, r, fy;        //Energy values and root
	int status;
	int itermax = 200;

	//Initialization of F
	F.function = &delta_h_nc_orbit; //Energy delta with respect to target
	F.params = &orbit;   //orbit structure contains the parameters

	//Bracketing the root
	s_low = 0.0;
	s_high = 1e-6 * orbit.pmap->si_max;  //arbitrary small value

	//Increasing s_high until a root is found (Hup*Hdown < 0)
	int iter = 0;
	do
	{
		//Evaluating the bracket
		Hup = delta_h_nc_orbit(s_high, &orbit);
		Hdown = delta_h_nc_orbit(s_low, &orbit);
		s_high *= 1.1;
		//iter++; --> commented: quit a bit dangerous but works in practice
	} while (Hup * Hdown > 0 && iter < itermax);

	//Swap the two variables if needed
	if (s_low > s_high)
	{
		s_low = s_low + s_high;
		s_high = s_low - s_high;
		s_low = s_low - s_high;
	}

	// TURN OFF GSL ERROR ERROR HANDLER
	// (GSL JUST PRINT ERROR MESSAGE AND KILL THE PROGRAM IF FLAG IS ON)
	gsl_set_error_handler_off();

	//If a root is found, refine root
	if (Hup * Hdown < 0)
	{
		//Setting the solver
		status = gsl_root_fsolver_set(orbit.ode_s_6->s_root, &F, s_low, s_high);
		//Loop
		iter = 0;
		do
		{
			status = gsl_root_fsolver_iterate(orbit.ode_s_6->s_root); //updating the solver
			r = gsl_root_fsolver_root(orbit.ode_s_6->s_root); //updating the root
			fy = delta_h_nc_orbit(r, &orbit);             //Checking convergence
			status = gsl_root_test_residual(fy, orbit.ode_s_6->eps_root); //Checking convergence
		} while (status == GSL_CONTINUE && (++iter) < 50);

		if (status == GSL_SUCCESS)
		{
			//Update st0
			orbit_update_s0(&orbit, st0, r);
			//Update the orbit once st0 is good
			orbit_update_ic(orbit, st0, orbit.pmap->t0_nc);
			return GSL_SUCCESS;
		}
		else
		{
			//cout <<  "orbit_init_pmap: No refined root was found (1). The orbit is unchanged" << endl;
			return GSL_FAILURE;
		}
	}
	else
	{
		//cout <<  "orbit_init_pmap: No root was found (2). The orbit is unchanged" << endl;
		return GSL_FAILURE;
	}
}

/**
 \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit
 *      given an array of initial TFC conditions si
 **/
void orbit_update_ic(Orbit& orbit, const double si[], double t0)
{
	//------------------------------------------------------------------------------------
	// 1. Update si
	//------------------------------------------------------------------------------------
	for (int p = 0; p < REDUCED_NV; p++)
		orbit.si[p] = si[p];

	//------------------------------------------------------------------------------------
	// 2. Update s0
	//------------------------------------------------------------------------------------
	rcm_to_ccm(si, orbit.s0, REDUCED_NV);

	//------------------------------------------------------------------------------------
	// 2. Update s0d
	//------------------------------------------------------------------------------------
	rcm_to_cmm8(si, orbit.s0d);

	//------------------------------------------------------------------------------------
	// 4. Update z0
	//------------------------------------------------------------------------------------
	//z0 = W(si, 0.0)
	rcm_to_nc_by_tfc(si, t0, orbit.fbpl->us.n, orbit.pmap->ofts_order,
			orbit.pmap->ofs_order, *orbit.Wh, *orbit.PC, *orbit.V, orbit.z0,
			orbit.pmap->graph_style);

	//------------------------------------------------------------------------------------
	// 5. Update zh0
	//------------------------------------------------------------------------------------
	//zh0 = wh(si, 0.0)
	rcm_to_tf(si, t0, orbit.fbpl->us.n, orbit.pmap->ofts_order,
			orbit.pmap->ofs_order, *orbit.Wh, *orbit.ofs, orbit.zh0,
			orbit.pmap->graph_style);
}

/**
 *   \brief Computes the hamiltonian at the position st0, in system coordinates and units.
 *   \param pmap a reference to the pmap that carries a set of useful parameters
 *   \param st0 the input state
 *
 *   WARNING: Direct use of the higher-order parameterization via rcm_to_nc_by_tfc
 **/
double h_nc_pmap(Pmap& pmap, double st0[])
{
	//------------------------------------------------------------------------------------
	// Inner variables (NC, TFC)
	//------------------------------------------------------------------------------------
	double zEM[6];
	double z0d[6];
	Ofsc AUX(OFS_ORDER);

	//------------------------------------------------------------------------------------
	// RCM to NC
	//------------------------------------------------------------------------------------
	rcm_to_nc_by_tfc(st0, pmap.t0_nc, SEML.us.n, pmap.ofts_order, pmap.ofs_order,
			CMh, Mcoc, Vcoc, z0d, pmap.graph_style);

	//------------------------------------------------------------------------------------
	// Computation
	//------------------------------------------------------------------------------------
	ncsys_m_to_sys_m(pmap.t0_nc, z0d, zEM, &SEML);
	return qbcp_H(pmap.t0_nc, zEM, &SEML);
}

/**
 \brief Computes the difference between an given Ham value and the state configuration defined in the routine orbit_update_s0 (see comments therein).
 \param  sr a double to complete the current tested configuration
 \param  params a pointer to the orbit with a given Ham value (why void? the idea is to a have a generic function, but might be useless at this point)
 \return the difference between the two hamiltonians
 **/
double delta_h_nc_orbit(double sr, void* params)
{
	Orbit* orbit = (Orbit*) params;
	double st0[4];

	//Update st0
	orbit_update_s0(orbit, st0, sr);

	//Return the hamiltonian value
	return (h_nc_pmap(*orbit->pmap, st0) - orbit->pmap->h_nc);
}

//----------------------------------------------------------------------------------------
//
// Computation of one orbit
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Computes the poincare map of one given orbit, with the following method:
 *          - Integration of the NC vector field and reset of the state with the
 *            semi-analytical approximation of the center manifold each time the
 *            xy-plane is crossed.
 *          - If orbit->pmap->f_proj_T > -1, there is no additional projection on the
 *            center manifold, other than when the xy-plane is crossed. Else, additional
 *            projections are performed between the crossings to help stabilize the
 *            computation, with frequency orbit->pmap->f_proj_T.
 *          - Note that all directions of crossings are considered.
 **/
int orbit_compute_pmap(Orbit* orbit, int method)
{
	//------------------------------------------------------------------------------------
	//Init
	//------------------------------------------------------------------------------------
	int status = GSL_FAILURE;

	//Starting in the right direction
	set_dir_ode_structure(orbit->ode_s_6, orbit->pmap->t0_nc,
			orbit->pmap->tmax_nc);

	//------------------------------------------------------------------------------------
	//Integration with NC vector field and reset of the state @ crossing of the xy-plane
	//------------------------------------------------------------------------------------
	status = single_pmap_proj(orbit);

	//------------------------------------------------------------------------------------
	//Reset
	//------------------------------------------------------------------------------------
	reset_ode_structure(orbit->ode_s_6);

	return status;
}

/**
 *   \brief Computes the Poincare map of one given orbit, with a single method:
 *          only the NC vector field is computed. The state is reset with the
 *          semi-analytical approximation of the center manifold each time the xy-plane
 *          is crossed.
 *
 *          Moreover, the state is also reset between the crossings, with a frequency
 *          pmap.f_proj_T given as a fraction of the SEM period T. If pmap.f_proj_T == -1,
 *          no additional projection is performed.
 **/
int single_pmap_proj(Orbit* orbit)
{
	//------------------------------------------------------------------------------------
	//Initialization of inner variables
	//------------------------------------------------------------------------------------
	int n_events = 0; 			        //number of events during integration
	double previous_zv_nc[6], zv_nc[6]; //for event detection purposes: state at step n and n-1
	double previous_t_nc, t_nc, tf_nc;  //times (event and integration)
	double sv_rcm[4];                   //for integration purposes, real TFC current configuration
	double zv_nc_mat[6];                //buffer to store the events (configuration right before the event)
	double t_nc_mat[2];                 //buffer to store the time of the events by pair of 2 (time before and after the event)
	value_output previous_sv_rcm;       //value output structure before the event
	value_output new_sv_rcm;            //value output structure after  the event
	double radius; 			            //current distance to origin

	//------------------------------------------------------------------------------------
	//Initialization of the state & time
	//------------------------------------------------------------------------------------
	for (int i = 0; i < 6; i++)
		zv_nc[i] = orbit->z0[i];
	for (int i = 0; i < 4; i++)
		sv_rcm[i] = orbit->si[i];
	orbit->reset_number = 0;
	tf_nc = orbit->pmap->tmax_nc;
	t_nc = orbit->pmap->t0_nc;

	//------------------------------------------------------------------------------------
	//Update of the method
	//------------------------------------------------------------------------------------
	orbit->int_method = SINGLE_INT;

	//------------------------------------------------------------------------------------
	//Find the projection parameters (in the long run: make it more general)
	//------------------------------------------------------------------------------------
	double omega1 = cimag(Fh[0].get_coef(1, 0)->ofs_get_coef(0));
	double omega3 = cimag(Fh[1].get_coef(1, 1)->ofs_get_coef(0));

	//------------------------------------------------------------------------------------
	//Integration and event detection up to t = tf or the number of events is maximal
	//------------------------------------------------------------------------------------
	//Projection tools
	double ePm;
	int nreset = 1;
	double status;
	do
	{
		//Keep track of step n-1
		previous_t_nc = t_nc;
		for (int i = 0; i < 6; i++)
			previous_zv_nc[i] = zv_nc[i];

		//Evolve one step of z(t)
		status = gslc_proj_step(orbit, zv_nc, &t_nc, orbit->pmap->t0_nc, tf_nc,
				&ePm, &nreset, omega1, omega3, (orbit->pmap->f_proj_T > 0));

		//Return if necessary
		if (status != GSL_SUCCESS)
			return GSL_FAILURE;

		//Event (crossing) detection
		previous_sv_rcm = orbit->fvalue->value(previous_t_nc, previous_zv_nc,
				orbit->fvalue->val_par);
		new_sv_rcm = orbit->fvalue->value(t_nc, zv_nc, orbit->fvalue->val_par);

		//A zero of the value function has been crossed
		if (previous_sv_rcm.val * new_sv_rcm.val < 0)
		{
			//Storage of values
			for (int i = 0; i < 6; i++)
				zv_nc_mat[i] = previous_zv_nc[i];
			t_nc_mat[0] = previous_t_nc;
			t_nc_mat[1] = t_nc;

			//Refine root so that yv is really on the z=0 plane
			status = refine_root(orbit, zv_nc, &t_nc, sv_rcm, zv_nc_mat,
					t_nc_mat, omega1, omega3, n_events);
			if (status != GSL_SUCCESS)
				return GSL_FAILURE;

			//Storage in orbit
			for (int i = 0; i < 6; i++)
				orbit->z0_mat[i][n_events] = zv_nc[i]; //Store position (NC)
			for (int i = 0; i < 4; i++)
				orbit->s0_mat[i][n_events] = sv_rcm[i]; //Store position (TFC)
			orbit->te_mat[n_events] = t_nc; //Store time
			orbit->nevent[n_events] = n_events + 1;  //Store the number of the event

			//Reset ode structure for next step
			reset_ode_structure(orbit->ode_s_6);

			//Add one event
			n_events++;
		}

		//Check if the system is diverging
		radius = sqrt(zv_nc[0] * zv_nc[0] + zv_nc[1] * zv_nc[1]+ zv_nc[2] * zv_nc[2]);
		if (radius > orbit->pmap->max_rad_div)
		{
			cout << "single_pmap_proj. the trajectory is diverging: radius = "
					<< radius << ". the solution is discarded." << endl;
			return GSL_FAILURE;
		}
	} while (status == GSL_SUCCESS && fabs(t_nc) < fabs(tf_nc)
			&& n_events < orbit->pmap->max_events);

	if (fabs(t_nc) >= fabs(tf_nc))
	{
		printf(
				"single_pmap_proj. Final time was reached, last state is returned at the end of ye.\n");
		for (int i = 0; i < 6; i++)
		{
			orbit->z0_mat[i][n_events] = zv_nc[i];
			orbit->te_mat[n_events] = t_nc;
			//WARNING: no TFC storage at this step!
		}
	}

	//------------------------------------------------------------------------------------
	//Last update of the orbit
	//------------------------------------------------------------------------------------
	//Last indix in the event storage
	if (n_events > 0)
		orbit->last_indix = n_events - 1;
	else
		orbit->last_indix = 0;

	return GSL_SUCCESS;
}

/**
 *   \brief Refine the root z = 0 for a single integration in NC coordinates. Used in single_pmap_proj.
 **/
int refine_root(Orbit* orbit, double* zv_nc, double* t_nc, double* sv_rcm, double* zv_nc_mat,
		double* t_nc_mat, double omega1, double omega3, int events)
{
	//------------------------------------------------------------------------------------
	//Initialization of inner variables
	//------------------------------------------------------------------------------------
	gsl_function F;             //called function in the root finding routine
	struct OdeParams params;    //parameters for F
	double zvs_nc[6];           //new start for root finding
	double t_nc_low, t_nc_high; //bracketing the root
	int status;                 //status of the search
	double root_gsl;            //root for GSL routine
	double res_gsl;             //error on the zero for GSL routine
	int iter;                   //accumulated iterations during root finding
	double ti;                  //inner routine time

	//------------------------------------------------------------------------------------
	// Init inner vectors
	//------------------------------------------------------------------------------------
	//Copy of initial values for storage in ye after the root finding
	for (int i = 0; i < 6; i++)
		zv_nc[i] = zv_nc_mat[i];

	//New start for root finding
	for (int i = 0; i < 6; i++)
		zvs_nc[i] = zv_nc_mat[i];

	//------------------------------------------------------------------------------------
	//Initialization of objects for the root finding
	//------------------------------------------------------------------------------------
	params.d = orbit->ode_s_6->d;
	params.fvalue = orbit->fvalue;
	F.function = &odezero_event;
	F.params = &params;
	params.t0 = t_nc_mat[0] - 1e-15;     //new initial time is previous_t - epsilon
	params.y0 = zvs_nc;

	//Bracketing the root
	t_nc_low = (t_nc_mat[0] < t_nc_mat[1]) ? t_nc_mat[0] : t_nc_mat[1];
	t_nc_high = (t_nc_mat[0] < t_nc_mat[1]) ? t_nc_mat[1] : t_nc_mat[0];

	//Setting the solver
	gsl_root_fsolver_set(orbit->ode_s_6->s_root, &F, t_nc_low, t_nc_high);

	//------------------------------------------------------------------------------------
	//Loop for root refinement
	//------------------------------------------------------------------------------------
	res_gsl = zvs_nc[1];
	iter = 0;
	do
	{
		status = gsl_root_fsolver_iterate(orbit->ode_s_6->s_root); //updating the solver
		root_gsl = gsl_root_fsolver_root(orbit->ode_s_6->s_root);   //updating the root
		//Checking convergence
		res_gsl = odezero_event(root_gsl, &params);
		status = gsl_root_test_residual(res_gsl, orbit->ode_s_6->eps_root);
	} while (status == GSL_CONTINUE && (++iter) < MAX_REF_ITER);

	if (iter >= MAX_REF_ITER)
	{
		cout
				<< "refine_root. WARNING: number of iter max exceeded in custom_odezero, with precision: "
				<< res_gsl << "Premature ending." << endl;
		return GSL_FAILURE;
	}

	//------------------------------------------------------------------------------------
	//Updating the outputs @time t = r
	//------------------------------------------------------------------------------------
	ti = t_nc_mat[0]; //set the time to previous step in order to integrate from ti to r
	status = gsl_odeiv2_driver_apply(orbit->ode_s_6->d, &ti, root_gsl, zv_nc); //updating zv_nc
	*t_nc = root_gsl;                                                          //updating t_nc

	//------------------------------------------------------------------------------------
	//Energy @ zv_nc(t) before projection
	//------------------------------------------------------------------------------------
	ncsys_m_to_sys_m(*t_nc, zv_nc, zvs_nc, (FBPL*) orbit->ode_s_6->d->sys->params);
	orbit->hz[events] = qbcp_H(*t_nc, zvs_nc, orbit->ode_s_6->d->sys->params);

	//------------------------------------------------------------------------------------
	// Project the state on the center manifold
	//------------------------------------------------------------------------------------
	//For comparison, the state before projection is stored in yvi
	double yvi[6];
	for (int i = 0; i < 6; i++)
		yvi[i] = zv_nc[i];

	//Project on the center manifold
	cdouble scp[4];
	nc_proj_cmm(zv_nc, root_gsl, SEML.us_em.n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3, scp,
			4);
	ccm_to_nc_by_tfc(scp, root_gsl, orbit->fbpl->us.n, orbit->pmap->ofts_order,
			orbit->pmap->ofs_order, *orbit->Wh, Mcoc, Vcoc, zv_nc,
			orbit->pmap->graph_style);

	//Compute projected state in RCM coordinates
	ccm_to_rcm(scp, sv_rcm, REDUCED_NV);

	//Get the current error
	double ePm = fabs(yvi[0] - zv_nc[0]);
	for (int i = 1; i < 6; i++)
	{
		if (fabs(yvi[i] - zv_nc[i]) > ePm)
			ePm = fabs(yvi[i] - zv_nc[i]);
	}

	//------------------------------------------------------------------------------------
	//Update the error in ePm
	//------------------------------------------------------------------------------------
	orbit->ePm[events] = ePm;

	//------------------------------------------------------------------------------------
	//Important: forcing exactly z = 0, otherwise a false event may be triggered @next step
	//------------------------------------------------------------------------------------
	zv_nc[2] = 0.0;

	//------------------------------------------------------------------------------------
	//Energy @ zv_nc(t) after projection
	//------------------------------------------------------------------------------------
	ncsys_m_to_sys_m(*t_nc, zv_nc, zvs_nc, (FBPL*) orbit->ode_s_6->d->sys->params);
	orbit->hw[events] = qbcp_H(*t_nc, zvs_nc, orbit->ode_s_6->d->sys->params);

	return GSL_SUCCESS;
}

//----------------------------------------------------------------------------------------
//
// Steppers with projection
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Integrates one step the current state yv using projection on the CM if necessary
 **/
int gslc_proj_step(Orbit* orbit, double zv_nc[], double* t_nc, double t0_nc, double t1_nc,
		double* ePm, int* n_reset, double omega1, double omega3, int is_reset_on)
{
	int status;
	double zv_nc_p[6], zv_nc_i[6];
	cdouble scp[REDUCED_NV];

	//------------------------------------------------------------------------------------
	//Evolve one step of z(t)
	//------------------------------------------------------------------------------------
	status = gsl_odeiv2_evolve_apply(orbit->ode_s_6->e, orbit->ode_s_6->c,
			orbit->ode_s_6->s, &orbit->ode_s_6->sys, t_nc, t1_nc, &orbit->ode_s_6->h,
			zv_nc);
	if (status != GSL_SUCCESS)
	{
		cout
				<< "gslc_proj_step. error: integration of z(t) has gone wrong. break."
				<< endl;
		return GSL_FAILURE;
	}

	//------------------------------------------------------------------------------------
	//Projection if necessary
	//------------------------------------------------------------------------------------
	if (is_reset_on && fabs(*t_nc - t0_nc) > fabs(*n_reset * orbit->pmap->t_proj_nc))
	{
		//------------------------------------------------------------------------------------
		// Projection on the center manifold
		//------------------------------------------------------------------------------------
		//Get the closest point on the center manifold
		nc_proj_cmm(zv_nc, *t_nc, SEML.us_em.n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3,
				scp, REDUCED_NV);
		//Update the state
		ccm_to_nc_by_tfc(scp, *t_nc, orbit->fbpl->us.n, orbit->pmap->ofts_order,
				orbit->pmap->ofs_order, *orbit->Wh, Mcoc, Vcoc, zv_nc_p,
				orbit->pmap->graph_style);
		//For comparison
		for (int i = 0; i < 6; i++)
			zv_nc_i[i] = zv_nc[i];
		// Copy of zv_nc_p in current state
		for (int i = 0; i < 6; i++)
			zv_nc[i] = zv_nc_p[i];

		//------------------------------------------------------------------------------------
		// Get the current projection error
		//------------------------------------------------------------------------------------
		//Get the current error (infinity norm)
		*ePm = fabs(zv_nc_i[0] - zv_nc[0]);
		for (int i = 1; i < 6; i++)
		{
			if (fabs(zv_nc_i[i] - zv_nc[i]) > *ePm)
				*ePm = fabs(zv_nc_i[i] - zv_nc[i]);
		}

		//------------------------------------------------------------------------------------
		//Reset ode structure for next step
		//------------------------------------------------------------------------------------
		reset_ode_structure(orbit->ode_s_6);

		//------------------------------------------------------------------------------------
		//One additional reset
		//------------------------------------------------------------------------------------
		*n_reset = *n_reset + 1;
	}

	return GSL_SUCCESS;
}

/**
 *   \brief Integrates the current state yv up to t = t1, using projection on the CM if necessary
 **/
int gslc_proj_evolve(Orbit* orbit, double zv_nc[], double* t_nc, double t0_nc, double t1_nc,
		double* ePm, int* n_reset, double omega1, double omega3, int is_reset_on)
{
	//Reset ode structure
	reset_ode_structure(orbit->ode_s_6);

	//Use gslc_proj_step to integrate the state from t0 to t1
	int status;
	do
	{
		status = gslc_proj_step(orbit, zv_nc, t_nc, t0_nc, t1_nc, ePm, n_reset, omega1,
				omega3, is_reset_on);

	} while (fabs(*t_nc) < fabs(t1_nc));

	return status;
}

//----------------------------------------------------------------------------------------
//
// Plotting
//
//----------------------------------------------------------------------------------------
/**
 *  \brief Plot the Poincaré map results stored in orbit on a 2-dim xy plot.
 *         Note that, since the orbit contains all points on the xy-plane (and not just the
 *         points for which pz > 0), only half of the results are plotted.
 **/
void orbit_poincare_plot(Orbit* orbit, gnuplot_ctrl* h1, int color)
{
	//Checking initialization of h1 & h2
	if (h1 == NULL)
		h1 = gnuplot_init();

	//Init & store
	double xe1[orbit->last_indix + 1];
	double xe2[orbit->last_indix + 1];

	int n_real_sol = 0;
	for (int i = 0; i <= orbit->last_indix; i++)
	{
		if (orbit->z0_mat[5][i] > 0)
		{
			xe1[n_real_sol] = orbit->z0_mat[0][i];
			xe2[n_real_sol] = orbit->z0_mat[1][i];
			n_real_sol++;
		}
	}

	//Plot
	gnuplot_set_xlabel(h1, (char*) "X");
	gnuplot_set_ylabel(h1, (char*) "Y");
	gnuplot_cmd(h1, "set title \"Poincare map (NC coordinates) \"");
	gnuplot_plot_xy(h1, xe1, xe2, n_real_sol, (char*) "", "points", "7", "2",
			color);
}

//----------------------------------------------------------------------------------------
//
//  Orbit C structure handling
//
//----------------------------------------------------------------------------------------
/**
 * \brief Initialize one orbit structure
 **/
void init_orbit(Orbit* orbit, vector<Oftsc>* W, vector<Oftsc>* Wh,
		matrix<Oftsc>* DW, matrix<Ofsc>* PC, vector<Ofsc>* V,
		value_params* val_par, value_function* fvalue, OdeStruct* ode_s_6,
		OdeStruct* ode_s_6_root, Pmap* pmap, FBPL* fbpl, Ofsc* orbit_ofs,
		int var_dim_h, int label)
{
	//------------------------------------------------------------------------------------
	//Parent
	//------------------------------------------------------------------------------------
	orbit->pmap = pmap;    //Poincare map (parent)
	orbit->fbpl = fbpl;    //QBCP around a given libration point (parent)
	orbit->n = fbpl->us.n; //Pulsation of the SEM system
	//------------------------------------------------------------------------------------
	//Parameterization (common to all orbits)
	//------------------------------------------------------------------------------------
	orbit->W = W;            //z(t) = W(s(t), t)
	orbit->Wh = Wh;          //zh(t) = Wh(s(t), t)
	orbit->DW = DW;          //Jacobian of W
	orbit->ofs = orbit_ofs;  //Auxiliary Ofs object

	//------------------------------------------------------------------------------------
	//COC (common to all orbits)
	//------------------------------------------------------------------------------------
	orbit->PC = PC; //COC matrix
	orbit->V = V;   //COC vector

	//------------------------------------------------------------------------------------
	//For event detection
	//------------------------------------------------------------------------------------
	orbit->val_par = val_par;  //Event parameters
	orbit->fvalue = fvalue;    //Value function for event detection

	orbit->z0_mat = dmatrix(0, 5, 0, val_par->max_events); //Matrix of states on the Poincaré map in NC coordinates
	orbit->zh0_mat = dmatrix(0, 5, 0, val_par->max_events); //Matrix of states on the Poincaré map in TFC coordinates
	orbit->s0_mat = dmatrix(0, REDUCED_NV - 1, 0, val_par->max_events); //Matrix of states on the Poincaré map in RCM coordinates
	orbit->te_mat = dvector(0, val_par->max_events); //Matrix of times on the Poincaré map
	orbit->hz = dvector(0, val_par->max_events); //Energy gap between z(t) and W(s(t),t) at each point on the map
	orbit->hw = dvector(0, val_par->max_events); //Energy gap between W(s(t),t) and the initial energy at each point on the map
	orbit->ePm = dvector(0, val_par->max_events); //Projection error at each event
	orbit->nevent = ivector(0, val_par->max_events); //Number of the event

	//------------------------------------------------------------------------------------
	//Characteristics
	//------------------------------------------------------------------------------------
	orbit->z0 = dvector(0, 5); //Initial state in NC coordinates dim = 6
	orbit->zh0 = dvector(0, 5);  //Initial state in TFC coordinates dim = 6
	orbit->si = dvector(0, REDUCED_NV - 1);  //Initial state in RCM coordinates dim = 4
	orbit->s0 = dcvector(0, REDUCED_NV - 1); //Initial state in CCM4 coordinates dim = 4
	orbit->s0d = dvector(0, 2 * REDUCED_NV - 1); //Initial state in CCM8 coordinates (real+imag part) dim = 8
	orbit->xf = dvector(0, 5); //Final state NC dim = 6
	orbit->tf = pmap->tmax_nc; //Final time after computation
	orbit->int_method = -1; //Integration method. -1 if not computed
	orbit->label = label; //Label of the orbit in the map
	orbit->var_dim_h = var_dim_h; //Modified dimension during energy adjustment

	//------------------------------------------------------------------------------------
	//ODE integration
	//------------------------------------------------------------------------------------
	orbit->ode_s_6 = ode_s_6;   //NC ode structure
	orbit->ode_s_6_root = ode_s_6_root; //NC ode structure (root finding)
}

/**
 \brief Free one orbit
 **/
void free_orbit(Orbit* orbit)
{
	//------------------------------------------------------------------------------------
	//Characteristics
	//------------------------------------------------------------------------------------
	free_dvector(orbit->z0, 0, 5);
	free_dvector(orbit->zh0, 0, 5);
	free_dvector(orbit->si, 0, REDUCED_NV - 1);
	free_dvector(orbit->s0d, 0, 2 * REDUCED_NV - 1);
	free_dvector(orbit->xf, 0, 5);
	free_dcvector(orbit->s0, 0, REDUCED_NV - 1);

	//------------------------------------------------------------------------------------
	//For event detection
	//------------------------------------------------------------------------------------
	free_dmatrix(orbit->z0_mat, 0, 5, 0, orbit->val_par->max_events);
	free_dmatrix(orbit->zh0_mat, 0, 5, 0, orbit->val_par->max_events);
	free_dmatrix(orbit->s0_mat, 0, REDUCED_NV - 1, 0,
			orbit->val_par->max_events);

	free_dvector(orbit->te_mat, 0, orbit->val_par->max_events);
	free_dvector(orbit->ePm, 0, orbit->val_par->max_events);
	free_dvector(orbit->hz, 0, orbit->val_par->max_events);
	free_dvector(orbit->hw, 0, orbit->val_par->max_events);
	free_ivector(orbit->nevent, 0, orbit->val_par->max_events);

	free_ode_structure(orbit->ode_s_6);
	free_ode_structure(orbit->ode_s_6_root);
}

