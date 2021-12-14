
#ifndef PATGRID_HEADER
#define PATGRID_HEADER

#include "patch.h"
#include "oflow.h"


namespace OFC
{

	class PatGridClass
	{

	public:
		PatGridClass(const camparam* cpt_in,
				const camparam* cpo_in,
				const optparam* op_in);

		~PatGridClass();

		void InitializeGrid(const float* im_ao_in, const float* im_ao_dx_in, const float* im_ao_dy_in);

		void SetTargetImage(const float* im_bo_in, const float* im_bo_dx_in, const float* im_bo_dy_in);

		void InitializeFromCoarserOF(const float* flow_prev);

		void AggregateFlowDense(float* flowout) const;


		void Optimize();



		void SetComplGrid(PatGridClass* cg_in);

		inline const int GetNoPatches() const
		{
			return nopatches;
		}

		inline const int GetNoph() const
		{
			return noph;
		}

		inline const int GetNopw() const
		{
			return nopw;
		}

		inline const Eigen::Vector2f GetRefPatchPos(int i) const
		{
			return pt_ref[i];
		}
		inline const Eigen::Vector2f GetQuePatchPos(int i) const
		{
			return pat[i]->GetPointPos();
		}
		inline const Eigen::Vector2f GetQuePatchDis(int i) const
		{
			return pt_ref[i] - pat[i]->GetPointPos();
		}



		void AggregateFlowDense(float* flowout, float* probabilityout, const bool bool_first, const bool bool_last,
				const float* spatial_prob, const float likelihood) const;


		void Optimize(const bool last_level);
		void InitializeFromCoarserOF(const float* flow_prev, const float* prob_prev);
		void removeSmallSegments_pat(const float speckle_size,const float speckle_sim_threshold);
		void removeSmallSegments(const float speckle_size, const float speckle_sim_threshold, float* flowout, float* probabilityout);

	private:

		const float* im_ao, * im_ao_dx, * im_ao_dy;
		const float* im_bo, * im_bo_dx, * im_bo_dy;

		Eigen::Map<const Eigen::MatrixXf>* im_ao_eg, * im_ao_dx_eg, * im_ao_dy_eg;
		Eigen::Map<const Eigen::MatrixXf>* im_bo_eg, * im_bo_dx_eg, * im_bo_dy_eg;

		const camparam* cpt;
		const camparam* cpo;
		const optparam* op;

		int steps;
		int nopw;
		int noph;
		int nopatches;


		std::vector<OFC::PatClass*> pat;
		std::vector<Eigen::Vector2f> pt_ref;
		std::vector<Eigen::Matrix<float, 1, 1>> p_init;
		const PatGridClass* cg = nullptr;


		std::vector<Eigen::Matrix<float, 1, 1>> prob_init;
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


