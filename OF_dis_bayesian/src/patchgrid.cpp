


#include <iostream>
#include <string>
#include <vector>
#include <valarray>

#include <thread>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <stdio.h>

#include "patch.h"
#include "patchgrid.h"


using std::cout;
using std::endl;
using std::vector;


namespace OFC
{

	PatGridClass::PatGridClass(
			const camparam* cpt_in,
			const camparam* cpo_in,
			const optparam* op_in)
			:
			cpt(cpt_in),
			cpo(cpo_in),
			op(op_in)
	{


		steps = op->steps;
		nopw = ceil((float)cpt->width / (float)steps);
		noph = ceil((float)cpt->height / (float)steps);
		const int offsetw = floor((cpt->width - (nopw - 1) * steps) / 2);
		const int offseth = floor((cpt->height - (noph - 1) * steps) / 2);

		nopatches = nopw * noph;
		pt_ref.resize(nopatches);
		p_init.resize(nopatches);
		prob_init.resize(nopatches);
		pat.reserve(nopatches);

		im_ao_eg = new Eigen::Map<const Eigen::MatrixXf>(nullptr, cpt->height, cpt->width);
		im_ao_dx_eg = new Eigen::Map<const Eigen::MatrixXf>(nullptr, cpt->height, cpt->width);
		im_ao_dy_eg = new Eigen::Map<const Eigen::MatrixXf>(nullptr, cpt->height, cpt->width);

		im_bo_eg = new Eigen::Map<const Eigen::MatrixXf>(nullptr, cpt->height, cpt->width);
		im_bo_dx_eg = new Eigen::Map<const Eigen::MatrixXf>(nullptr, cpt->height, cpt->width);
		im_bo_dy_eg = new Eigen::Map<const Eigen::MatrixXf>(nullptr, cpt->height, cpt->width);

		int patchid = 0;
		for (int x = 0; x < nopw; ++x)
		{
			for (int y = 0; y < noph; ++y)
			{
				int i = x * noph + y;

				pt_ref[i][0] = x * steps + offsetw;
				pt_ref[i][1] = y * steps + offseth;
				p_init[i].setZero();
				prob_init[i].setZero();

				pat.push_back(new OFC::PatClass(cpt, cpo, op, patchid));
				patchid++;
			}
		}


	}

	PatGridClass::~PatGridClass()
	{
		delete im_ao_eg;
		delete im_ao_dx_eg;
		delete im_ao_dy_eg;

		delete im_bo_eg;
		delete im_bo_dx_eg;
		delete im_bo_dy_eg;

		for (int i = 0; i < nopatches; ++i)
			delete pat[i];

	}

	void PatGridClass::SetComplGrid(PatGridClass* cg_in)
	{
		cg = cg_in;
	}

	void PatGridClass::InitializeGrid(const float* im_ao_in, const float* im_ao_dx_in, const float* im_ao_dy_in)
	{
		im_ao = im_ao_in;
		im_ao_dx = im_ao_dx_in;
		im_ao_dy = im_ao_dy_in;

		new(im_ao_eg) Eigen::Map<const Eigen::MatrixXf>(im_ao, cpt->height, cpt->width);
		new(im_ao_dx_eg) Eigen::Map<const Eigen::MatrixXf>(im_ao_dx, cpt->height, cpt->width);
		new(im_ao_dy_eg) Eigen::Map<const Eigen::MatrixXf>(im_ao_dy, cpt->height, cpt->width);


#pragma omp parallel for schedule(static)
		for (int i = 0; i < nopatches; ++i)
		{
			pat[i]->InitializePatch(im_ao_eg, im_ao_dx_eg, im_ao_dy_eg, pt_ref[i]);
			p_init[i].setZero();
			prob_init[i].setZero();
		}


	}

	void PatGridClass::SetTargetImage(const float* im_bo_in, const float* im_bo_dx_in, const float* im_bo_dy_in)
	{
		im_bo = im_bo_in;
		im_bo_dx = im_bo_dx_in;
		im_bo_dy = im_bo_dy_in;

		new(im_bo_eg) Eigen::Map<const Eigen::MatrixXf>(im_bo, cpt->height, cpt->width);
		new(im_bo_dx_eg) Eigen::Map<const Eigen::MatrixXf>(im_bo_dx, cpt->height, cpt->width);
		new(im_bo_dy_eg) Eigen::Map<const Eigen::MatrixXf>(im_bo_dy, cpt->height, cpt->width);

#pragma omp parallel for schedule(static)
		for (int i = 0; i < nopatches; ++i)
			pat[i]->SetTargetImage(im_bo_eg, im_bo_dx_eg, im_bo_dy_eg);
	}

	void PatGridClass::Optimize(const bool first_level)
	{
#pragma omp parallel for schedule(dynamic, 10)
		for (int i = 0; i < nopatches; ++i)
		{


			if (p_init[i](0) != 0 || first_level == true)
			{
				float probability = prob_init[i](0);
				if (first_level)
					pat[i]->OptimizeIter(p_init[i], true);
				else if (probability > 0)
					pat[i]->OptimizeIter(p_init[i], true);
			}
		}

	}

	void PatGridClass::removeSmallSegments_pat(const float speckle_size, const float speckle_sim_threshold)
	{


		int32_t D_width = nopw;
		int32_t D_height = noph;
		int32_t D_speckle_size = speckle_size;


		int32_t* D_done = (int32_t*)calloc(D_width * D_height, sizeof(int32_t));
		int32_t* seg_list_u = (int32_t*)calloc(D_width * D_height, sizeof(int32_t));
		int32_t* seg_list_v = (int32_t*)calloc(D_width * D_height, sizeof(int32_t));
		int32_t seg_list_count;
		int32_t seg_list_curr;
		int32_t u_neighbor[4];
		int32_t v_neighbor[4];
		int32_t u_seg_curr;
		int32_t v_seg_curr;


		int32_t addr_start, addr_curr, addr_neighbor;
		addr_curr = 1;

		for (int32_t u = 0; u < D_width; u++)
		{
			for (int32_t v = 0; v < D_height; v++)
			{


				addr_start = v * D_width + u;


				if (*(D_done + addr_start) == 0 && (pat[addr_start]->get_bayesian_prob()) != 0)
				{


					*(seg_list_u + 0) = u;
					*(seg_list_v + 0) = v;
					seg_list_count = 1;
					seg_list_curr = 0;


					while (seg_list_curr < seg_list_count)
					{


						u_seg_curr = *(seg_list_u + seg_list_curr);
						v_seg_curr = *(seg_list_v + seg_list_curr);


						addr_curr = v_seg_curr * D_width + u_seg_curr;


						u_neighbor[0] = u_seg_curr - 1;
						v_neighbor[0] = v_seg_curr;
						u_neighbor[1] = u_seg_curr + 1;
						v_neighbor[1] = v_seg_curr;
						u_neighbor[2] = u_seg_curr;
						v_neighbor[2] = v_seg_curr - 1;
						u_neighbor[3] = u_seg_curr;
						v_neighbor[3] = v_seg_curr + 1;


						for (int32_t i = 0; i < 4; i++)
						{


							if (u_neighbor[i] >= 0 && v_neighbor[i] >= 0 && u_neighbor[i] < D_width &&
								v_neighbor[i] < D_height)
							{

								addr_neighbor = v_neighbor[i] * D_width + u_neighbor[i];


								if (*(D_done + addr_neighbor) == 0 && (pat[addr_neighbor]->get_bayesian_prob()) != 0)
								{

									if (abs((*(pat[addr_curr]->GetParam()))[0] -
											(*(pat[addr_neighbor]->GetParam()))[0]) <= speckle_sim_threshold)
									{


										*(seg_list_u + seg_list_count) = u_neighbor[i];
										*(seg_list_v + seg_list_count) = v_neighbor[i];
										seg_list_count++;


										*(D_done + addr_neighbor) = 1;
									}
								}

							}
						}


						seg_list_curr++;


						*(D_done + addr_curr) = 1;

					}


					if (seg_list_count <= D_speckle_size)
					{


						for (int32_t i = 0; i < seg_list_count; i++)
						{
							addr_curr = *(seg_list_v + i) * D_width + *(seg_list_u + i);
							pat[int(addr_curr)]->set_bayesian_prob(0);
						}
					}
				}

			}
		}


		free(D_done);
		free(seg_list_u);
		free(seg_list_v);
	}


	void PatGridClass::removeSmallSegments(const float speckle_size, const float speckle_sim_threshold, float* flowout,
			float* probabilityout)
	{


		int32_t D_width = cpt->width;
		int32_t D_height = cpt->height;
		int32_t D_speckle_size = speckle_size;


		int32_t* D_done = (int32_t*)calloc(D_width * D_height, sizeof(int32_t));
		int32_t* seg_list_u = (int32_t*)calloc(D_width * D_height, sizeof(int32_t));
		int32_t* seg_list_v = (int32_t*)calloc(D_width * D_height, sizeof(int32_t));
		int32_t seg_list_count;
		int32_t seg_list_curr;
		int32_t u_neighbor[4];
		int32_t v_neighbor[4];
		int32_t u_seg_curr;
		int32_t v_seg_curr;


		int32_t addr_start, addr_curr, addr_neighbor;
		addr_curr = 1;


		for (int32_t u = 0; u < D_width; u++)
		{
			for (int32_t v = 0; v < D_height; v++)
			{


				addr_start = v * D_width + u;


				if (*(D_done + addr_start) == 0 && *(flowout + addr_start) != 0)
				{


					*(seg_list_u + 0) = u;
					*(seg_list_v + 0) = v;
					seg_list_count = 1;
					seg_list_curr = 0;


					while (seg_list_curr < seg_list_count)
					{


						u_seg_curr = *(seg_list_u + seg_list_curr);
						v_seg_curr = *(seg_list_v + seg_list_curr);


						addr_curr = v_seg_curr * D_width + u_seg_curr;


						u_neighbor[0] = u_seg_curr - 1;
						v_neighbor[0] = v_seg_curr;
						u_neighbor[1] = u_seg_curr + 1;
						v_neighbor[1] = v_seg_curr;
						u_neighbor[2] = u_seg_curr;
						v_neighbor[2] = v_seg_curr - 1;
						u_neighbor[3] = u_seg_curr;
						v_neighbor[3] = v_seg_curr + 1;


						for (int32_t i = 0; i < 4; i++)
						{


							if (u_neighbor[i] >= 0 && v_neighbor[i] >= 0 && u_neighbor[i] < D_width &&
								v_neighbor[i] < D_height)
							{


								addr_neighbor = v_neighbor[i] * D_width + u_neighbor[i];


								if (*(D_done + addr_neighbor) == 0 && flowout[addr_neighbor] != 0)
								{


									if (abs((flowout[addr_curr] - flowout[addr_neighbor])) <= speckle_sim_threshold)
									{


										*(seg_list_u + seg_list_count) = u_neighbor[i];
										*(seg_list_v + seg_list_count) = v_neighbor[i];
										seg_list_count++;
										*(D_done + addr_neighbor) = 1;
									}
								}

							}
						}


						seg_list_curr++;
						*(D_done + addr_curr) = 1;

					}
					if (seg_list_count <= D_speckle_size)
					{


						for (int32_t i = 0; i < seg_list_count; i++)
						{
							addr_curr = *(seg_list_v + i) * D_width + *(seg_list_u + i);
							flowout[addr_curr] = std::numeric_limits<float>::quiet_NaN();;
							probabilityout[addr_curr] = 0;
						}
					}
				}

			}
		}


		free(D_done);
		free(seg_list_u);
		free(seg_list_v);
	}

	void PatGridClass::InitializeFromCoarserOF(const float* flow_prev, const float* prob_prev)
	{
#pragma omp parallel for schedule(dynamic, 10)
		for (int ip = 0; ip < nopatches; ++ip)
		{
			int x = floor(pt_ref[ip][0] / 2);
			int y = floor(pt_ref[ip][1] / 2);
			int i = y * (cpt->width / 2) + x;
			p_init[ip](0) = flow_prev[i] * 2;
			prob_init[ip](0) = prob_prev[i];
		}

	}


	void PatGridClass::AggregateFlowDense(float* flowout) const
	{
		float* we = new float[cpt->width * cpt->height];

		memset(flowout, 0, sizeof(float) * (op->nop * cpt->width * cpt->height));
		memset(we, 0, sizeof(float) * (cpt->width * cpt->height));

#ifdef USE_PARALLEL_ON_FLOWAGGR
#pragma omp parallel for schedule(static)
#endif
		for (int ip = 0; ip < nopatches; ++ip)
		{

			if (pat[ip]->IsValid())
			{
				const Eigen::Matrix<float, 1, 1>* fl = pat[ip]->GetParam();
				Eigen::Matrix<float, 1, 1> flnew;
				const float* pweight = pat[ip]->GetpWeightPtr();

				int lb = -op->p_samp_s / 2;
				int ub = op->p_samp_s / 2 - 1;

				for (int y = lb; y <= ub; ++y)
				{
					for (int x = lb; x <= ub; ++x, ++pweight)
					{
						int yt = (y + pt_ref[ip][1]);
						int xt = (x + pt_ref[ip][0]);

						if (xt >= 0 && yt >= 0 && xt < cpt->width && yt < cpt->height)
						{

							int i = yt * cpt->width + xt;
							float absw = 1.0f / (float)(std::max(op->minerrval, *pweight));

							flnew = (*fl) * absw;
							we[i] += absw;
							flowout[i] += flnew[0];
						}
					}
				}
			}
		}


		if (cg)
		{
			Eigen::Vector4f wbil;
			Eigen::Vector4i pos;

#ifdef USE_PARALLEL_ON_FLOWAGGR
#pragma omp parallel for schedule(static)
#endif
			for (int ip = 0; ip < cg->nopatches; ++ip)
			{
				if (cg->pat[ip]->IsValid())
				{
					const Eigen::Matrix<float, 1, 1>* fl = (cg->pat[ip]->GetParam());
					Eigen::Matrix<float, 1, 1> flnew;

					const Eigen::Vector2f rppos = cg->pat[ip]->GetPointPos();
					const float* pweight = cg->pat[ip]->GetpWeightPtr();

					Eigen::Vector2f resid;


					pos[0] = ceil(rppos[0] + .00001);
					pos[1] = ceil(rppos[1] + .00001);
					pos[2] = floor(rppos[0]);
					pos[3] = floor(rppos[1]);

					resid[0] = rppos[0] - pos[2];
					resid[1] = rppos[1] - pos[3];
					wbil[0] = resid[0] * resid[1];
					wbil[1] = (1 - resid[0]) * resid[1];
					wbil[2] = resid[0] * (1 - resid[1]);
					wbil[3] = (1 - resid[0]) * (1 - resid[1]);

					int lb = -op->p_samp_s / 2;
					int ub = op->p_samp_s / 2 - 1;


					for (int y = lb; y <= ub; ++y)
					{
						for (int x = lb; x <= ub; ++x, ++pweight)
						{

							int yt = y + pos[1];
							int xt = x + pos[0];
							if (xt >= 1 && yt >= 1 && xt < (cpt->width - 1) && yt < (cpt->height - 1))
							{

								float absw = 1.0f / (float)(std::max(op->minerrval, *pweight));


								flnew = (*fl) * absw;

								int idxcc = xt + yt * cpt->width;
								int idxfc = (xt - 1) + yt * cpt->width;
								int idxcf = xt + (yt - 1) * cpt->width;
								int idxff = (xt - 1) + (yt - 1) * cpt->width;

								we[idxcc] += wbil[0] * absw;
								we[idxfc] += wbil[1] * absw;
								we[idxcf] += wbil[2] * absw;
								we[idxff] += wbil[3] * absw;

								flowout[idxcc] -=
										wbil[0] * flnew[0];
								flowout[idxfc] -= wbil[1] * flnew[0];
								flowout[idxcf] -= wbil[2] * flnew[0];
								flowout[idxff] -= wbil[3] * flnew[0];
							}
						}
					}
				}
			}
		}

#pragma omp parallel for schedule(static, 100)

		for (int yi = 0; yi < cpt->height; ++yi)
		{
			for (int xi = 0; xi < cpt->width; ++xi)
			{
				int i = yi * cpt->width + xi;
				if (we[i] > 0)
				{
					flowout[i] /= we[i];
				}
			}
		}

		delete[] we;
	}


	void
	PatGridClass::AggregateFlowDense(float* flowout, float* probabilityout, const bool bool_first, const bool bool_last,
			const float* spatial_prob, const float likelihood) const
	{

		float* we = new float[cpt->width * cpt->height];
		unsigned int* num_pat = new unsigned int[cpt->width * cpt->height];
		float* var_tmp = new float[cpt->width * cpt->height];

		memset(flowout, 0, sizeof(float) * (op->nop * cpt->width * cpt->height));
		memset(probabilityout, 0, sizeof(float) * (op->nop * cpt->width * cpt->height));
		memset(we, 0, sizeof(float) * (cpt->width * cpt->height));
		memset(num_pat, 0, sizeof(unsigned int) * (cpt->width * cpt->height));

		if (op->bool_var && bool_last)
		{
			memset(var_tmp, 0, sizeof(float) * (op->nop * cpt->width * cpt->height));
		}


#ifdef USE_PARALLEL_ON_FLOWAGGR
#pragma omp parallel for schedule(static)
#endif
		for (int ip = 0; ip < nopatches; ++ip)
		{

			if (pat[ip]->IsValid())
			{
				const Eigen::Matrix<float, 1, 1>* fl = pat[ip]->GetParam();
				Eigen::Matrix<float, 1, 1> flnew;
				const float bayesian_prob = pat[ip]->get_bayesian_prob();

				int lb = -op->p_samp_s / 2;
				int ub = op->p_samp_s / 2 - 1;


				float probability;
				if (bool_first)
					probability = bayesian_prob * likelihood;
				else
					probability = bayesian_prob * likelihood + prob_init[ip](0) * 1;


				double var_pat = 0;
				if (op->bool_var && bool_last)
				{
					pat[ip]->getVarMAP(&var_pat);
				}


				if (bayesian_prob == 0)
					continue;
				for (int y = lb; y <= ub; ++y)
				{
					for (int x = lb; x <= ub; ++x)
					{
						int yt = (y + pt_ref[ip][1]);
						int xt = (x + pt_ref[ip][0]);

						if (xt >= 0 && yt >= 0 && xt < cpt->width && yt < cpt->height)
						{
							int i = yt * cpt->width + xt;
							num_pat[i]++;
							int ind_in_patch = (x - lb) * op->p_samp_s + y - lb;

							flnew = (*fl) * probability * spatial_prob[ind_in_patch];
							we[i] += probability * spatial_prob[ind_in_patch];
/*							if(op->bool_var && bool_last)
								probabilityout[i] += pow(probability * spatial_prob[ind_in_patch],2)*var_pat;
							else*/
							if (op->bool_var && bool_last)
								var_tmp[i] += pow(probability * spatial_prob[ind_in_patch], 2) * var_pat;
							probabilityout[i] += probability;
							flowout[i] += flnew[0];
						}
					}
				}
			}
		}


#pragma omp parallel for schedule(static, 100)

		for (int yi = 0; yi < cpt->height; ++yi)
		{
			for (int xi = 0; xi < cpt->width; ++xi)
			{
				int i = yi * cpt->width + xi;

				if (we[i] > 0 && num_pat[i] > 0)
				{
					flowout[i] /= we[i];
					if (op->bool_var && bool_last)
						var_tmp[i] /= we[i] * we[i];
					probabilityout[i] /= num_pat[i];
				}
				else
				{
					flowout[i] = std::numeric_limits<float>::quiet_NaN();
					probabilityout[i] = 0;
				}
			}
		}

		if (bool_last)
		{
			for (int i = 0; i < op->nop * cpt->width * cpt->height; ++i)
			{
				if (probabilityout[i] <= op->min_prob)
				{
					probabilityout[i] = std::numeric_limits<float>::quiet_NaN();
					flowout[i] = std::numeric_limits<float>::quiet_NaN();
				}
				else
				{
					probabilityout[i] = -var_tmp[i];
				}
			}
		}

		delete[] we;
		delete[] num_pat;
		if (op->bool_var && bool_last)
			delete[] var_tmp;
	}
}


