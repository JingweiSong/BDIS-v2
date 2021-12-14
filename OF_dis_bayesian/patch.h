


#ifndef PAT_HEADER
#define PAT_HEADER


#include "oflow.h"

namespace OFC
{


	typedef struct
	{
		bool hasconverged;
		bool hasoptstarted;


		Eigen::Matrix<float, Eigen::Dynamic, 1> pdiff;
		Eigen::Matrix<float, Eigen::Dynamic, 1> pweight;


		Eigen::Matrix<float, 1, 1> Hes;
		Eigen::Matrix<float, 1, 1> p_in, p_iter, delta_p;


		Eigen::Matrix<float, 1, 1> normtmp;
		Eigen::Vector2f pt_iter;
		Eigen::Vector2f pt_st;

		float delta_p_sqnorm = 1e-10;
		float delta_p_sqnorm_init = 1e-10;
		float mares = 1e20;
		float mares_old = 1e20;
		int cnt = 0;
		bool invalid = false;
	} patchstate;


	class PatClass
	{

	public:
		PatClass(const camparam* cpt_in,
				const camparam* cpo_in,
				const optparam* op_in,
				const int patchid_in);

		~PatClass();

		void
		InitializePatch(Eigen::Map<const Eigen::MatrixXf>* im_ao_in, Eigen::Map<const Eigen::MatrixXf>* im_ao_dx_in,
				Eigen::Map<const Eigen::MatrixXf>* im_ao_dy_in, const Eigen::Vector2f pt_ref_in);

		void SetTargetImage(Eigen::Map<const Eigen::MatrixXf>* im_bo_in, Eigen::Map<const Eigen::MatrixXf>* im_bo_dx_in,
				Eigen::Map<const Eigen::MatrixXf>* im_bo_dy_in);

		void OptimizeIter(const Eigen::Matrix<float, 1, 1> p_in_arg, const bool untilconv);

		inline const bool isConverged() const
		{
			return pc->hasconverged;
		}

		inline const bool hasOptStarted() const
		{
			return pc->hasoptstarted;
		}

		inline const Eigen::Vector2f GetPointPos() const
		{
			return pc->pt_iter;
		}

		inline const bool IsValid() const
		{
			return (!pc->invalid);
		}

		inline const bool SetInValid()
		{
			pc->invalid = false;
		}

		inline const float* GetpWeightPtr() const
		{
			return (float*)pc->pweight.data();
		}

		inline const Eigen::Matrix<float, 1, 1>* GetParam() const
		{
			return &(pc->p_iter);
		}

		inline const Eigen::Matrix<float, 1, 1>* GetParamStart() const
		{
			return &(pc->p_in);
		}

		inline const float get_bayesian_prob() const
		{
			return *bayesian_prob;
		}

		inline void set_bayesian_prob(const float prob)
		{
			*bayesian_prob = prob;
		}

		inline void set_last_level_true()
		{
			last_level = true;
		}

		void getVarMAP(double* var);

	private:


		void OptimizeStart(const Eigen::Matrix<float, 1, 1> p_in_arg);

		void OptimizeComputeErrImg();

		void paramtopt();

		void ResetPatch();

		void ComputeHessian();

		void CreateStatusStruct(patchstate* psin);

		void LossComputeErrorImage(Eigen::Matrix<float, Eigen::Dynamic, 1>* patdest,
				Eigen::Matrix<float, Eigen::Dynamic, 1>* wdest, const Eigen::Matrix<float, Eigen::Dynamic, 1>* patin,
				const Eigen::Matrix<float, Eigen::Dynamic, 1>* tmpin);


		bool
		getPatchStaticNNGrad(const float* img, const float* img_dx, const float* img_dy, const Eigen::Vector2f* mid_in,
				Eigen::Matrix<float, Eigen::Dynamic, 1>* tmp_in, Eigen::Matrix<float, Eigen::Dynamic, 1>* tmp_dx_in,
				Eigen::Matrix<float, Eigen::Dynamic, 1>* tmp_dy_in);

		bool getPatchStaticBil(const float* img, const Eigen::Vector2f* mid_in,
				Eigen::Matrix<float, Eigen::Dynamic, 1>* tmp_in_e);

		Eigen::Vector2f pt_ref;
		Eigen::Matrix<float, Eigen::Dynamic, 1> tmp;
		Eigen::Matrix<float, Eigen::Dynamic, 1> dxx_tmp;
		Eigen::Matrix<float, Eigen::Dynamic, 1> dyy_tmp;

		Eigen::Map<const Eigen::MatrixXf>* im_ao, * im_ao_dx, * im_ao_dy;
		Eigen::Map<const Eigen::MatrixXf>* im_bo, * im_bo_dx, * im_bo_dy;

		const camparam* cpt;
		const camparam* cpo;
		const optparam* op;
		const int patchid;

		patchstate* pc = nullptr;

		bool MonteCarloConfidence(const int num_window, const float unit_disturb, float* bayesian_prob);

		void LossComputeErrorImage_noNaN(Eigen::Matrix<float, Eigen::Dynamic, 1>* patdest,
				Eigen::Matrix<float, Eigen::Dynamic, 1>* wdest, const Eigen::Matrix<float, Eigen::Dynamic, 1>* patin,
				const Eigen::Matrix<float, Eigen::Dynamic, 1>* tmpin);

		bool last_level;
		int num_valid;
		Eigen::Matrix<float, Eigen::Dynamic, 1> tmp_nomeannorm;
		float* bayesian_prob;

		inline double fast_exp(double y)
		{
			double d;
			*(reinterpret_cast<int*>(&d) + 0) = 0;
			*(reinterpret_cast<int*>(&d) + 1) = static_cast<int>(1512775 * y + 1072632447);
			return d;
		}

	};


}

#endif


