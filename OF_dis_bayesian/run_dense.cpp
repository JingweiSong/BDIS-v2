
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <iostream>
#include <sys/time.h>
#include <fstream>

#include "oflow.h"


using namespace std;


void SavePFMFile_jingwei(const cv::Mat& img, const cv::Mat& probability, const float min_prob, const char* filename)
{
	ofstream myfile(filename);
	cv::Size szt = img.size();
	for (int y = 0; y < szt.height; ++y)
	{
		for (int x = 0; x < szt.width; ++x)
		{
			float tmp = -img.at<float>(y, x);
			float prob = probability.at<float>(y, x);

			if (std::isnan(tmp))
				myfile << 0 << ",";
			else
				myfile << tmp << ",";
		}
		myfile << endl;
	}
	myfile.close();

}

void ReadFlowFile(cv::Mat& img, const char* filename)
{
	FILE* stream = fopen(filename, "rb");
	if (stream == 0)
		cout << "ReadFile: could not open %s" << endl;

	int width, height;
	float tag;
	int nc = img.channels();
	float tmp[nc];

	if ((int)fread(&tag, sizeof(float), 1, stream) != 1 ||
		(int)fread(&width, sizeof(int), 1, stream) != 1 ||
		(int)fread(&height, sizeof(int), 1, stream) != 1)
		cout << "ReadFile: problem reading file %s" << endl;

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			if ((int)fread(tmp, sizeof(float), nc, stream) != nc)
				cout << "ReadFile(%s): file is too short" << endl;

			if (nc == 1)
				img.at<float>(y, x) = tmp[0];
			else if (nc == 2)
			{
				img.at<cv::Vec2f>(y, x)[0] = tmp[0];
				img.at<cv::Vec2f>(y, x)[1] = tmp[1];
			}
			else if (nc == 4)
			{
				img.at<cv::Vec4f>(y, x)[0] = tmp[0];
				img.at<cv::Vec4f>(y, x)[1] = tmp[1];
				img.at<cv::Vec4f>(y, x)[2] = tmp[2];
				img.at<cv::Vec4f>(y, x)[3] = tmp[3];
			}
		}
	}

	if (fgetc(stream) != EOF)
		cout << "ReadFile(%s): file is too long" << endl;

	fclose(stream);
}

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

int main(int argc, char** argv)
{
	struct timeval tv_start_all, tv_end_all;
	gettimeofday(&tv_start_all, NULL);

	std::string folderpath(argv[1]);
	std::string filename(argv[2]);
	std::string outfilename(argv[3]);
	string imgfile_ao = folderpath + "/left/" + filename + ".png";
	string imgfile_bo = folderpath + "/right/" + filename + ".png";
	string outfile = folderpath + "/" + "bayesian_" + outfilename + filename;
	string outfile_prob = folderpath + "/" + "prob_out_" + outfilename + filename;


	cv::Mat img_ao_mat, img_bo_mat, img_tmp;
	int rpyrtype, nochannels, incoltype;
	incoltype = CV_LOAD_IMAGE_GRAYSCALE;
	rpyrtype = CV_32FC1;
	nochannels = 1;
	img_ao_mat = cv::imread(imgfile_ao, incoltype);
	img_bo_mat = cv::imread(imgfile_bo, incoltype);
	cv::Mat img_ao_fmat, img_bo_fmat;
	cv::Size sz = img_ao_mat.size();
	int width_org = sz.width;
	int height_org = sz.height;


	int lv_f, lv_l, maxiter, miniter, patchsz, patnorm, costfct, tv_innerit, tv_solverit, verbosity;
	float mindprate, mindrrate, minimgerr, poverl, tv_alpha, tv_gamma, tv_delta, tv_sor;
	float ratio_patch_valid;
	int num_window;
	float unit_disturb;
	float min_prob;
	bool bool_var;
	bool usefbcon, usetvref;


	if (argc <= 5)
	{
		mindprate = 0.05;
		mindrrate = 0.95;
		minimgerr = 0.0;
		usefbcon = 0;
		patnorm = 1;
		costfct = 0;
		tv_alpha = 10.0;
		tv_gamma = 10.0;
		tv_delta = 5.0;
		tv_innerit = 1;
		tv_solverit = 3;
		tv_sor = 1.6;
		ratio_patch_valid = 1.0;
		num_window = 2;
		unit_disturb = 0.5;
		min_prob = 0;
		bool_var = false;
		verbosity = 2;

		int fratio = 5;

		int sel_oppoint = 2;
		if (argc == 5)
			sel_oppoint = atoi(argv[4]);

		switch (sel_oppoint)
		{
		case 1:
			patchsz = 8;
			poverl = 0.3;
			lv_f = AutoFirstScaleSelect(width_org, fratio, patchsz);
			lv_l = std::max(lv_f - 2, 0);
			maxiter = 16;
			miniter = 16;
			usetvref = 0;
			break;
		case 3:
			patchsz = 12;
			poverl = 0.75;
			lv_f = AutoFirstScaleSelect(width_org, fratio, patchsz);
			lv_l = std::max(lv_f - 4, 0);
			maxiter = 16;
			miniter = 16;
			usetvref = 1;
			break;
		case 4:
			patchsz = 12;
			poverl = 0.75;
			lv_f = AutoFirstScaleSelect(width_org, fratio, patchsz);
			lv_l = std::max(lv_f - 5, 0);
			maxiter = 128;
			miniter = 128;
			usetvref = 1;
			break;
		case 2:
		default:
			patchsz = 8;
			poverl = 0.4;
			lv_f = AutoFirstScaleSelect(width_org, fratio, patchsz);
			lv_l = std::max(lv_f - 2, 0);
			maxiter = 12;
			miniter = 12;
			usetvref = 1;
			break;

		}
	}
	else
	{
		int acnt = 4;
		lv_f = atoi(argv[acnt++]);
		lv_l = atoi(argv[acnt++]);
		maxiter = atoi(argv[acnt++]);
		miniter = atoi(argv[acnt++]);
		mindprate = atof(argv[acnt++]);
		mindrrate = atof(argv[acnt++]);
		minimgerr = atof(argv[acnt++]);
		patchsz = atoi(argv[acnt++]);
		poverl = atof(argv[acnt++]);
		usefbcon = atoi(argv[acnt++]);
		patnorm = atoi(argv[acnt++]);
		costfct = atoi(argv[acnt++]);
		verbosity = atoi(argv[acnt++]);
		ratio_patch_valid = atof(argv[acnt++]);
		num_window = atoi(argv[acnt++]);
		unit_disturb = atof(argv[acnt++]);
		min_prob = atof(argv[acnt++]);
		if (atof(argv[acnt++]) != 0)
			bool_var = true;
	}


	int padw = 0, padh = 0;
	int scfct = pow(2, lv_f);

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


	if (verbosity > 1)
	{
		gettimeofday(&tv_end_all, NULL);
		double tt = (tv_end_all.tv_sec - tv_start_all.tv_sec) * 1000.0f +
					(tv_end_all.tv_usec - tv_start_all.tv_usec) / 1000.0f;
		printf("TIME (Image loading     ) (ms): %3g\n", tt);
		gettimeofday(&tv_start_all, NULL);
	}


	img_ao_mat.convertTo(img_ao_fmat, CV_32F);
	img_bo_mat.convertTo(img_bo_fmat, CV_32F);

	const float* img_ao_pyr[lv_f + 1];
	const float* img_bo_pyr[lv_f + 1];
	const float* img_ao_dx_pyr[lv_f + 1];
	const float* img_ao_dy_pyr[lv_f + 1];
	const float* img_bo_dx_pyr[lv_f + 1];
	const float* img_bo_dy_pyr[lv_f + 1];

	cv::Mat img_ao_fmat_pyr[lv_f + 1];
	cv::Mat img_bo_fmat_pyr[lv_f + 1];
	cv::Mat img_ao_dx_fmat_pyr[lv_f + 1];
	cv::Mat img_ao_dy_fmat_pyr[lv_f + 1];
	cv::Mat img_bo_dx_fmat_pyr[lv_f + 1];
	cv::Mat img_bo_dy_fmat_pyr[lv_f + 1];


	if (SELECTCHANNEL == 1)
	{
		img_ao_fmat.setTo(std::numeric_limits<float>::quiet_NaN(), img_ao_fmat == 0);
		img_ao_fmat.setTo(std::numeric_limits<float>::quiet_NaN(), img_ao_fmat == 0);
	}


	ConstructImgPyramide(img_ao_fmat, img_ao_fmat_pyr, img_ao_dx_fmat_pyr, img_ao_dy_fmat_pyr, img_ao_pyr,
			img_ao_dx_pyr, img_ao_dy_pyr, lv_f, lv_l, rpyrtype, 1, patchsz, padw, padh);
	ConstructImgPyramide(img_bo_fmat, img_bo_fmat_pyr, img_bo_dx_fmat_pyr, img_bo_dy_fmat_pyr, img_bo_pyr,
			img_bo_dx_pyr, img_bo_dy_pyr, lv_f, lv_l, rpyrtype, 1, patchsz, padw, padh);


	if (verbosity > 1)
	{
		gettimeofday(&tv_end_all, NULL);
		double tt = (tv_end_all.tv_sec - tv_start_all.tv_sec) * 1000.0f +
					(tv_end_all.tv_usec - tv_start_all.tv_usec) / 1000.0f;
		printf("TIME (Pyramide+Gradients) (ms): %3g\n", tt);
	}

	float sc_fct = pow(2, lv_l);
	cv::Mat flowout(sz.height / sc_fct, sz.width / sc_fct, CV_32FC1);
	cv::Mat probabilityout(sz.height / sc_fct, sz.width / sc_fct, CV_32FC1);

	OFC::OFClass ofc(img_ao_pyr, img_ao_dx_pyr, img_ao_dy_pyr,
			img_bo_pyr, img_bo_dx_pyr, img_bo_dy_pyr,
			patchsz,
			(float*)flowout.data,
			(float*)probabilityout.data,
			nullptr,
			sz.width, sz.height,
			lv_f, lv_l, maxiter, miniter, mindprate, mindrrate, minimgerr, patchsz, poverl,
			usefbcon, costfct, nochannels, patnorm,
			verbosity,
			ratio_patch_valid,
			num_window,
			unit_disturb,
			bool_var,
			min_prob);

	if (verbosity > 1) gettimeofday(&tv_start_all, NULL);


	if (lv_l != 0)
	{
		flowout *= sc_fct;
		cv::resize(flowout, flowout, cv::Size(), sc_fct, sc_fct, cv::INTER_LINEAR);
		cv::resize(probabilityout, probabilityout, cv::Size(), sc_fct, sc_fct, cv::INTER_LINEAR);
	}


	flowout = flowout(cv::Rect((int)floor((float)padw / 2.0f), (int)floor((float)padh / 2.0f), width_org, height_org));
	probabilityout = probabilityout(
			cv::Rect((int)floor((float)padw / 2.0f), (int)floor((float)padh / 2.0f), width_org, height_org));
	SavePFMFile_jingwei(flowout, probabilityout, min_prob, outfile.c_str());

	if (bool_var)
		SavePFMFile_jingwei(probabilityout, probabilityout, min_prob, outfile_prob.c_str());


	if (verbosity > 1)
	{
		gettimeofday(&tv_end_all, NULL);
		double tt = (tv_end_all.tv_sec - tv_start_all.tv_sec) * 1000.0f +
					(tv_end_all.tv_usec - tv_start_all.tv_usec) / 1000.0f;
		printf("TIME (Saving flow file  ) (ms): %3g\n", tt);
	}

	return 0;
}


    


