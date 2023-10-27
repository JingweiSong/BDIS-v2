

#include <iostream>
#include <string>
#include <vector>

#include <thread>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>


#include <sys/time.h>
#include <stdio.h>

#include "oflow.h"
#include "patchgrid.h"

using std::cout;
using std::endl;
using std::vector;

namespace OFC
{

	OFClass::OFClass(const float** im_ao_in, const float** im_ao_dx_in, const float** im_ao_dy_in,
			const float** im_bo_in, const float** im_bo_dx_in, const float** im_bo_dy_in,
			const int imgpadding_in,
			float* outflow,
			float* outprobability,
			const float* initflow,
			const int width_in, const int height_in,
			const int sc_f_in, const int sc_l_in,
			const int max_iter_in, const int min_iter_in,
			const float dp_thresh_in,
			const float dr_thresh_in,
			const float res_thresh_in,
			const int p_samp_s_in,
			const float patove_in,
			const bool usefbcon_in,
			const int costfct_in,
			const int noc_in,
			const int patnorm_in,
			const int verbosity_in,
			const float ratio_patch_valid,
			const int num_window,
			const float unit_disturb,
			const bool bool_var,
			const float min_prob)
			: im_ao(im_ao_in), im_ao_dx(im_ao_dx_in), im_ao_dy(im_ao_dy_in),
			  im_bo(im_bo_in), im_bo_dx(im_bo_dx_in), im_bo_dy(im_bo_dy_in)
	{


#ifdef WITH_OPENMP
		if (verbosity_in>1)
		  cout <<  "OPENMP is ON - used in pconst, pinit, potim";
#ifdef USE_PARALLEL_ON_FLOWAGGR
		if (verbosity_in>1)
		  cout << ", cflow ";
#endif
		if (verbosity_in>1) cout << endl;
#endif

		op.nop = 1;
		op.p_samp_s = p_samp_s_in;
		op.outlierthresh = (float)op.p_samp_s / 2;
		op.patove = patove_in;
		op.sc_f = sc_f_in;
		op.sc_l = sc_l_in;
		op.max_iter = max_iter_in;
		op.min_iter = min_iter_in;
		op.dp_thresh = dp_thresh_in * dp_thresh_in;
		op.dr_thresh = dr_thresh_in;
		op.res_thresh = res_thresh_in;
		op.steps = std::max(1, (int)floor(op.p_samp_s * (1 - op.patove)));
		op.novals = noc_in * (p_samp_s_in) * (p_samp_s_in);
		op.usefbcon = usefbcon_in;
		op.costfct = costfct_in;
		op.noc = noc_in;
		op.patnorm = patnorm_in;
		op.verbosity = verbosity_in;
		op.noscales = op.sc_f - op.sc_l + 1;


		op.normoutlier_tmpbsq = (v4sf){ op.normoutlier * op.normoutlier, op.normoutlier * op.normoutlier,
										op.normoutlier * op.normoutlier, op.normoutlier * op.normoutlier };
		op.normoutlier_tmp2bsq = __builtin_ia32_mulps(op.normoutlier_tmpbsq, op.twos);
		op.normoutlier_tmp4bsq = __builtin_ia32_mulps(op.normoutlier_tmpbsq, op.fours);


		op.ratio_patch_valid = ratio_patch_valid;
		op.num_window = num_window;
		op.unit_disturb = unit_disturb;
		op.bool_var = bool_var;
		op.min_prob = min_prob;


		double sigma_0_squared_2_spatial = 2 * 4 * 4;
		spatial_prob = new float[op.p_samp_s * op.p_samp_s];
		for (int i = 0; i < op.p_samp_s; i++)
		{
			for (int j = 0; j < op.p_samp_s; j++)
			{
				float length_square =
						std::pow(op.p_samp_s / 2 - 0.5 - i, 2) + std::pow(op.p_samp_s / 2 - 0.5 - j, 2);
				*(spatial_prob + i * op.p_samp_s + j) = fast_exp(-length_square / sigma_0_squared_2_spatial);
			}
		}


		struct timeval tv_start_all, tv_end_all, tv_start_all_global, tv_end_all_global;
		if (op.verbosity > 0)
			gettimeofday(&tv_start_all_global, nullptr);


		double tt_patconstr[op.noscales], tt_patinit[op.noscales], tt_patoptim[op.noscales], tt_compflow[op.noscales], tt_tvopt[op.noscales], tt_all[op.noscales];
		for (int sl = op.sc_f; sl >= op.sc_l; --sl)
		{
			tt_patconstr[sl - op.sc_l] = 0;
			tt_patinit[sl - op.sc_l] = 0;
			tt_patoptim[sl - op.sc_l] = 0;
			tt_compflow[sl - op.sc_l] = 0;
			tt_tvopt[sl - op.sc_l] = 0;
			tt_all[sl - op.sc_l] = 0;
		}

		if (op.verbosity > 1) gettimeofday(&tv_start_all, nullptr);


		vector<OFC::PatGridClass*> grid_fw(op.noscales);
		vector<OFC::PatGridClass*> grid_bw(op.noscales);
		vector<float*> flow_fw(op.noscales);
		vector<float*> prob_fw(op.noscales);
		vector<float*> flow_bw(op.noscales);
		cpl.resize(op.noscales);
		cpr.resize(op.noscales);
		for (int sl = op.sc_f; sl >= op.sc_l; --sl)
		{
			int i = sl - op.sc_l;

			float sc_fct = pow(2, -sl);
			cpl[i].sc_fct = sc_fct;
			cpl[i].height = height_in * sc_fct;
			cpl[i].width = width_in * sc_fct;
			cpl[i].imgpadding = imgpadding_in;
			cpl[i].tmp_lb = -(float)op.p_samp_s / 2;
			cpl[i].tmp_ubw = (float)(cpl[i].width + op.p_samp_s / 2 - 2);
			cpl[i].tmp_ubh = (float)(cpl[i].height + op.p_samp_s / 2 - 2);
			cpl[i].tmp_w = cpl[i].width + 2 * imgpadding_in;
			cpl[i].tmp_h = cpl[i].height + 2 * imgpadding_in;
			cpl[i].curr_lv = sl;
			cpl[i].camlr = 0;


			cpr[i] = cpl[i];
			cpr[i].camlr = 1;

			flow_fw[i] = new float[op.nop * cpl[i].width * cpl[i].height];
			prob_fw[i] = new float[op.nop * cpl[i].width * cpl[i].height];


			grid_fw[i] = new OFC::PatGridClass(&(cpl[i]), &(cpr[i]), &op);

			if (op.usefbcon)
			{
				flow_bw[i] = new float[op.nop * cpr[i].width * cpr[i].height];
				grid_bw[i] = new OFC::PatGridClass(&(cpr[i]), &(cpl[i]), &op);


				grid_fw[i]->SetComplGrid(grid_bw[i]);
				grid_bw[i]->SetComplGrid(grid_fw[i]);
			}
		}


		if (op.verbosity > 1)
		{
			gettimeofday(&tv_end_all, nullptr);
			double tt_gridconst = (tv_end_all.tv_sec - tv_start_all.tv_sec) * 1000.0f +
								  (tv_end_all.tv_usec - tv_start_all.tv_usec) / 1000.0f;
			printf("TIME (Grid Memo. Alloc. ) (ms): %3g\n", tt_gridconst);
		}


		float* likelihood = new float[op.sc_f - op.sc_l + 1];
		float likelihood_sum = 0;
		for (int sl = op.sc_f; sl >= op.sc_l; --sl)
		{
			*(likelihood + op.sc_f - sl) = 1 / float(sl);
			likelihood_sum = likelihood_sum + *(likelihood + op.sc_f - sl);
		}
		for (int sl = op.sc_f; sl >= op.sc_l; --sl)
		{
			*(likelihood + op.sc_f - sl) = *(likelihood + op.sc_f - sl) / likelihood_sum;
		}


		for (int sl = op.sc_f; sl >= op.sc_l; --sl)
		{
			int ii = sl - op.sc_l;

			if (op.verbosity > 1) gettimeofday(&tv_start_all, nullptr);


			grid_fw[ii]->InitializeGrid(im_ao[sl], im_ao_dx[sl], im_ao_dy[sl]);
			grid_fw[ii]->SetTargetImage(im_bo[sl], im_bo_dx[sl], im_bo_dy[sl]);
			if (op.usefbcon)
			{
				grid_bw[ii]->InitializeGrid(im_bo[sl], im_bo_dx[sl], im_bo_dy[sl]);
				grid_bw[ii]->SetTargetImage(im_ao[sl], im_ao_dx[sl], im_ao_dy[sl]);
			}


			if (op.verbosity > 1)
			{
				gettimeofday(&tv_end_all, nullptr);
				tt_patconstr[ii] = (tv_end_all.tv_sec - tv_start_all.tv_sec) * 1000.0f +
								   (tv_end_all.tv_usec - tv_start_all.tv_usec) / 1000.0f;
				tt_all[ii] += tt_patconstr[ii];
				gettimeofday(&tv_start_all, nullptr);
			}


			if (sl < op.sc_f)
			{
				grid_fw[ii]->InitializeFromCoarserOF(flow_fw[ii + 1], prob_fw[ii + 1]);


				if (op.usefbcon)
					grid_bw[ii]->InitializeFromCoarserOF(flow_bw[ii + 1], prob_fw[ii + 1]);
			}
			else if (sl == op.sc_f && initflow != nullptr)
			{
				grid_fw[ii]->InitializeFromCoarserOF(initflow, prob_fw[ii + 1]);
			}


			if (op.verbosity > 1)
			{
				gettimeofday(&tv_end_all, nullptr);
				tt_patinit[ii] = (tv_end_all.tv_sec - tv_start_all.tv_sec) * 1000.0f +
								 (tv_end_all.tv_usec - tv_start_all.tv_usec) / 1000.0f;
				tt_all[ii] += tt_patinit[ii];
				gettimeofday(&tv_start_all, nullptr);
			}


			if (sl == op.sc_f)
				grid_fw[ii]->Optimize(true);
			else
				grid_fw[ii]->Optimize(false);

			if (op.verbosity > 1)
			{
				gettimeofday(&tv_end_all, nullptr);
				tt_patoptim[ii] = (tv_end_all.tv_sec - tv_start_all.tv_sec) * 1000.0f +
								  (tv_end_all.tv_usec - tv_start_all.tv_usec) / 1000.0f;
				tt_all[ii] += tt_patoptim[ii];

				gettimeofday(&tv_start_all, nullptr);
			}

			float* tmp_ptr = flow_fw[ii];
			float* tmp_ptr_prob = prob_fw[ii];
			if (sl == op.sc_l)
			{
				tmp_ptr = outflow;
				tmp_ptr_prob = outprobability;
			}

			if (sl == op.sc_l)
				grid_fw[ii]->AggregateFlowDense(tmp_ptr, tmp_ptr_prob, false, true, spatial_prob,
						*(likelihood + op.sc_f - sl));
			else if (sl == op.sc_f)
				grid_fw[ii]->AggregateFlowDense(tmp_ptr, tmp_ptr_prob, true, false, spatial_prob,
						*(likelihood + op.sc_f - sl));
			else
				grid_fw[ii]->AggregateFlowDense(tmp_ptr, tmp_ptr_prob, false, false, spatial_prob,
						*(likelihood + op.sc_f - sl));

			grid_fw[ii]->removeSmallSegments(pow(4, op.sc_f - sl), 1, tmp_ptr, tmp_ptr_prob);


			if (op.usefbcon && sl > op.sc_l)
				grid_bw[ii]->AggregateFlowDense(flow_bw[ii]);


			if (op.verbosity > 1)
			{
				gettimeofday(&tv_end_all, nullptr);
				tt_compflow[ii] = (tv_end_all.tv_sec - tv_start_all.tv_sec) * 1000.0f +
								  (tv_end_all.tv_usec - tv_start_all.tv_usec) / 1000.0f;
				tt_all[ii] += tt_compflow[ii];

				gettimeofday(&tv_start_all, nullptr);
			}


			if (op.verbosity > 1)
			{
				gettimeofday(&tv_end_all, nullptr);
				tt_tvopt[ii] = (tv_end_all.tv_sec - tv_start_all.tv_sec) * 1000.0f +
							   (tv_end_all.tv_usec - tv_start_all.tv_usec) / 1000.0f;
				tt_all[ii] += tt_tvopt[ii];
				printf("TIME (Sc: %i, #p:%6i, pconst, pinit, poptim, cflow, tvopt, total): %8.2f %8.2f %8.2f %8.2f %8.2f -> %8.2f ms.\n",
						sl, grid_fw[ii]->GetNoPatches(), tt_patconstr[ii], tt_patinit[ii], tt_patoptim[ii],
						tt_compflow[ii], tt_tvopt[ii], tt_all[ii]);
			}


		}


		delete[] likelihood;
		for (int sl = op.sc_f; sl >= op.sc_l; --sl)
		{

			delete[] flow_fw[sl - op.sc_l];
			delete grid_fw[sl - op.sc_l];

			delete[] prob_fw[sl - op.sc_l];


			if (op.usefbcon)
			{
				delete[] flow_bw[sl - op.sc_l];
				delete grid_bw[sl - op.sc_l];
			}
		}


		if (op.verbosity > 0)
		{
			gettimeofday(&tv_end_all_global, nullptr);
			double tt = (tv_end_all_global.tv_sec - tv_start_all_global.tv_sec) * 1000.0f +
						(tv_end_all_global.tv_usec - tv_start_all_global.tv_usec) / 1000.0f;
			printf("TIME (O.Flow Run-Time   ) (ms): %3g\n", tt);
		}


	}

	OFClass::~OFClass()
	{
		delete[] spatial_prob;
	}


}














