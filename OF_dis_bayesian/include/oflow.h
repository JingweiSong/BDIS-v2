


#ifndef OFC_HEADER
#define OFC_HEADER

using std::cout;
using std::endl;

namespace OFC
{

	typedef __v4sf v4sf;


	typedef struct
	{
		int width;
		int height;
		int imgpadding;
		float tmp_lb;
		float tmp_ubw;
		float tmp_ubh;
		int tmp_w;
		int tmp_h;
		float sc_fct;
		int curr_lv;
		int camlr;
	} camparam;

	typedef struct
	{

		int sc_f;
		int sc_l;
		int p_samp_s;
		int max_iter;
		int min_iter;
		float dp_thresh;
		float dr_thresh;
		float res_thresh;
		int patnorm;
		int verbosity;
		bool usefbcon;
		int costfct;
		bool usetvref;
		float tv_alpha;
		float tv_gamma;
		float tv_delta;
		int tv_innerit;
		int tv_solverit;
		float tv_sor;


		int nop;
		float patove;
		float outlierthresh;
		int steps;
		int novals;
		int noc;
		int noscales;
		float minerrval = 2.0f;
		float normoutlier = 5.0f;


		v4sf zero = (v4sf){ 0.0f, 0.0f, 0.0f, 0.0f };
		v4sf negzero = (v4sf){ -0.0f, -0.0f, -0.0f, -0.0f };
		v4sf half = (v4sf){ 0.5f, 0.5f, 0.5f, 0.5f };
		v4sf ones = (v4sf){ 1.0f, 1.0f, 1.0f, 1.0f };
		v4sf twos = (v4sf){ 2.0f, 2.0f, 2.0f, 2.0f };
		v4sf fours = (v4sf){ 4.0f, 4.0f, 4.0f, 4.0f };
		v4sf normoutlier_tmpbsq;
		v4sf normoutlier_tmp2bsq;
		v4sf normoutlier_tmp4bsq;


		float ratio_patch_valid = 1.0;
		int num_window = 2;
		float unit_disturb = 0.5;
		bool bool_var = false;
		float min_prob = 0;    //	Minimal probability for depth filtering

	} optparam;


	class OFClass
	{

	public:
		OFClass(const float** im_ao_in, const float** im_ao_dx_in, const float** im_ao_dy_in,
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
				const int padval_in,
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
				const float min_prob);

	public:

		~OFClass();

	private:


		const float** im_ao, ** im_ao_dx, ** im_ao_dy;
		const float** im_bo, ** im_bo_dx, ** im_bo_dy;

		optparam op;
		std::vector<camparam> cpl, cpr;


		inline double fast_exp(double y)
		{
			double d;
			*(reinterpret_cast<int*>(&d) + 0) = 0;
			*(reinterpret_cast<int*>(&d) + 1) = static_cast<int>(1512775 * y + 1072632447);
			return d;
		}

		float* spatial_prob;

	};


}

#endif


