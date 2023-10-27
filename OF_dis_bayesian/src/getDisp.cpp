
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <iostream>
#include <sys/time.h>
#include <fstream>

#include "getDisp.h"
#include "oflow.h"


using namespace std;



void ConstructImgPyramide(const cv::Mat& img_ao_fmat, cv::Mat* img_ao_fmat_pyr, cv::Mat* img_ao_dx_fmat_pyr,
		cv::Mat* img_ao_dy_fmat_pyr, const float** img_ao_pyr, const float** img_ao_dx_pyr, const float** img_ao_dy_pyr,
		const int lv_f, const int lv_l, const int rpyrtype, const bool getgrad, const int imgpadding, const int padw,
		const int padh)
{
	for (int i = 0; i <= lv_f; ++i)
	{
		if (i == 0)
		{
			img_ao_fmat_pyr[i] = img_ao_fmat.clone();
		}
		else
			cv::resize(img_ao_fmat_pyr[i - 1], img_ao_fmat_pyr[i], cv::Size(), .5, .5, cv::INTER_LINEAR);

		img_ao_fmat_pyr[i].convertTo(img_ao_fmat_pyr[i], rpyrtype);

		if (getgrad)
		{
			cv::Sobel(img_ao_fmat_pyr[i], img_ao_dx_fmat_pyr[i], CV_32F, 1, 0, 1, 1, 0, cv::BORDER_DEFAULT);
			cv::Sobel(img_ao_fmat_pyr[i], img_ao_dy_fmat_pyr[i], CV_32F, 0, 1, 1, 1, 0, cv::BORDER_DEFAULT);
			img_ao_dx_fmat_pyr[i].convertTo(img_ao_dx_fmat_pyr[i], CV_32F);
			img_ao_dy_fmat_pyr[i].convertTo(img_ao_dy_fmat_pyr[i], CV_32F);
		}
	}


	for (int i = 0; i <= lv_f; ++i)
	{

		copyMakeBorder(img_ao_fmat_pyr[i], img_ao_fmat_pyr[i], imgpadding, imgpadding, imgpadding, imgpadding,
				cv::BORDER_CONSTANT, 0);
		img_ao_pyr[i] = (float*)img_ao_fmat_pyr[i].data;

		if (getgrad)
		{
			copyMakeBorder(img_ao_dx_fmat_pyr[i], img_ao_dx_fmat_pyr[i], imgpadding, imgpadding, imgpadding, imgpadding,
					cv::BORDER_CONSTANT, 0);
			copyMakeBorder(img_ao_dy_fmat_pyr[i], img_ao_dy_fmat_pyr[i], imgpadding, imgpadding, imgpadding, imgpadding,
					cv::BORDER_CONSTANT, 0);

			img_ao_dx_pyr[i] = (float*)img_ao_dx_fmat_pyr[i].data;
			img_ao_dy_pyr[i] = (float*)img_ao_dy_fmat_pyr[i].data;
		}

		for (int i = 0; i <= lv_f; ++i)
		{
			cv::patchNaNs(img_ao_fmat_pyr[i], 0);
			cv::patchNaNs(img_ao_dx_fmat_pyr[i], 0);
			cv::patchNaNs(img_ao_dy_fmat_pyr[i], 0);
		}
	}
}

int AutoFirstScaleSelect(int imgwidth, int fratio, int patchsize)
{
	return std::max(0, (int)std::floor(log2((2.0f * (float)imgwidth) / ((float)fratio * (float)patchsize))));
}

void getDisp(cv::Mat& img_ao_mat, cv::Mat& img_bo_mat, int& selectchannel, int& rpyrtype, ParamBDIS& param, cv::Mat& flowout, cv::Mat& probabilityout)
{
	struct timeval tv_start_all, tv_end_all;
	gettimeofday(&tv_start_all, NULL);

	cv::Mat img_ao_fmat, img_bo_fmat;
	cv::Size sz = img_ao_mat.size();
	int width_org = sz.width;
	int height_org = sz.height;


	//	Added by Jingwei:
	int nochannels = selectchannel;


	int padw = 0, padh = 0;
	int scfct = pow(2, param.lv_f);

	int div = sz.width % scfct;
	if (div > 0) padw = scfct - div;
	div = sz.height % scfct;
	if (div > 0) padh = scfct - div;
	if (padh > 0 || padw > 0)
	{
		copyMakeBorder(img_ao_mat, img_ao_mat, floor((float)padh / 2.0f), ceil((float)padh / 2.0f),
				floor((float)padw / 2.0f), ceil((float)padw / 2.0f), cv::BORDER_CONSTANT, 0);
		copyMakeBorder(img_bo_mat, img_bo_mat, floor((float)padh / 2.0f), ceil((float)padh / 2.0f),
				floor((float)padw / 2.0f), ceil((float)padw / 2.0f), cv::BORDER_CONSTANT, 0);
	}
	sz = img_ao_mat.size();


	if (param.verbosity > 1)
	{
		gettimeofday(&tv_end_all, NULL);
		double tt = (tv_end_all.tv_sec - tv_start_all.tv_sec) * 1000.0f +
					(tv_end_all.tv_usec - tv_start_all.tv_usec) / 1000.0f;
		printf("TIME (Image loading     ) (ms): %3g\n", tt);
		gettimeofday(&tv_start_all, NULL);
	}


	img_ao_mat.convertTo(img_ao_fmat, CV_32F);
	img_bo_mat.convertTo(img_bo_fmat, CV_32F);

	const float* img_ao_pyr[param.lv_f + 1];
	const float* img_bo_pyr[param.lv_f + 1];
	const float* img_ao_dx_pyr[param.lv_f + 1];
	const float* img_ao_dy_pyr[param.lv_f + 1];
	const float* img_bo_dx_pyr[param.lv_f + 1];
	const float* img_bo_dy_pyr[param.lv_f + 1];

	cv::Mat img_ao_fmat_pyr[param.lv_f + 1];
	cv::Mat img_bo_fmat_pyr[param.lv_f + 1];
	cv::Mat img_ao_dx_fmat_pyr[param.lv_f + 1];
	cv::Mat img_ao_dy_fmat_pyr[param.lv_f + 1];
	cv::Mat img_bo_dx_fmat_pyr[param.lv_f + 1];
	cv::Mat img_bo_dy_fmat_pyr[param.lv_f + 1];


	if (selectchannel == 1)
	{
		img_ao_fmat.setTo(std::numeric_limits<float>::quiet_NaN(), img_ao_fmat == 0);
		img_ao_fmat.setTo(std::numeric_limits<float>::quiet_NaN(), img_ao_fmat == 0);
	}


	ConstructImgPyramide(img_ao_fmat, img_ao_fmat_pyr, img_ao_dx_fmat_pyr, img_ao_dy_fmat_pyr, img_ao_pyr,
			img_ao_dx_pyr, img_ao_dy_pyr, param.lv_f, param.lv_l, rpyrtype, 1, param.patchsz, padw, padh);
	ConstructImgPyramide(img_bo_fmat, img_bo_fmat_pyr, img_bo_dx_fmat_pyr, img_bo_dy_fmat_pyr, img_bo_pyr,
			img_bo_dx_pyr, img_bo_dy_pyr, param.lv_f, param.lv_l, rpyrtype, 1, param.patchsz, padw, padh);


	if (param.verbosity > 1)
	{
		gettimeofday(&tv_end_all, NULL);
		double tt = (tv_end_all.tv_sec - tv_start_all.tv_sec) * 1000.0f +
					(tv_end_all.tv_usec - tv_start_all.tv_usec) / 1000.0f;
		printf("TIME (Pyramide+Gradients) (ms): %3g\n", tt);
	}

	float sc_fct = pow(2, param.lv_l);
	flowout = cv::Mat(sz.height / sc_fct, sz.width / sc_fct, CV_32FC1);
	probabilityout = cv::Mat(sz.height / sc_fct, sz.width / sc_fct, CV_32FC1);

	OFC::OFClass ofc(img_ao_pyr, img_ao_dx_pyr, img_ao_dy_pyr,
			img_bo_pyr, img_bo_dx_pyr, img_bo_dy_pyr,
			param.patchsz,
			(float*)flowout.data,
			(float*)probabilityout.data,
			nullptr,
			sz.width, sz.height,
			param.lv_f, param.lv_l, param.maxiter, param.miniter, param.mindprate, param.mindrrate, param.minimgerr, param.patchsz, param.poverl,
			param.usefbcon, param.costfct, nochannels, param.patnorm,
			param.verbosity,
			param.ratio_patch_valid,
			param.num_window,
			param.unit_disturb,
			param.bool_var,
			param.min_prob);

	if (param.verbosity > 1) gettimeofday(&tv_start_all, NULL);


	if (param.lv_l != 0)
	{
		flowout *= sc_fct;
		cv::resize(flowout, flowout, cv::Size(), sc_fct, sc_fct, cv::INTER_LINEAR);
		cv::resize(probabilityout, probabilityout, cv::Size(), sc_fct, sc_fct, cv::INTER_LINEAR);
	}


	flowout = flowout(cv::Rect((int)floor((float)padw / 2.0f), (int)floor((float)padh / 2.0f), width_org, height_org));
	probabilityout = probabilityout(cv::Rect((int)floor((float)padw / 2.0f), (int)floor((float)padh / 2.0f), width_org, height_org));


	//return flowout;
}


    


