#ifndef __GETDISP_H__
#define __GETDISP_H__

struct ParamBDIS{
	int lv_f = 5;
	int lv_l = 1;
	int maxiter = 12;
	int miniter = 12;
	float mindprate = 0.05;
	float mindrrate = 0.95;
	float minimgerr = 0;
	int patchsz = 10;
	float poverl = 0.60;
	int usefbcon = 0;
	int patnorm = 1;
	int costfct = 0;
	int verbosity = 2;
	float ratio_patch_valid = 0.25;
	int num_window = 2;
	float unit_disturb = 0.5;
	float min_prob = 0.05;
	bool bool_var = false;
};

void ConstructImgPyramide(const cv::Mat& img_ao_fmat, cv::Mat* img_ao_fmat_pyr, cv::Mat* img_ao_dx_fmat_pyr,
		cv::Mat* img_ao_dy_fmat_pyr, const float** img_ao_pyr, const float** img_ao_dx_pyr, const float** img_ao_dy_pyr,
		const int lv_f, const int lv_l, const int rpyrtype, const bool getgrad, const int imgpadding, const int padw,
		const int padh);

int AutoFirstScaleSelect(int imgwidth, int fratio, int patchsize);

void getDisp(cv::Mat& img_ao_mat, cv::Mat& img_bo_mat, int& selectchannel, int& rpyrtype, ParamBDIS& param, cv::Mat& flowout, cv::Mat& probabilityout);

#endif


