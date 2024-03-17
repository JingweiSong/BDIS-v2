
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <iostream>
#include <sys/time.h>
#include <fstream>

#include "getDisp.h"


using namespace std;


/*void SavePFMFile(const cv::Mat& img, const cv::Mat& probability, const float min_prob, const char* filename)
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

}*/
void SavePFMFile(const cv::Mat& img, const char* filename)
{
	cv::Size szt = img.size();
	cv::Mat out(szt.height, szt.width, CV_16UC1);
	float scale = -10;
	for (int y = 0; y < szt.height; ++y)
	{
		for (int x = 0; x < szt.width; ++x)
		{
			float tmp = img.at<float>(y, x);
			if (std::isnan(tmp))
				out.at<uint16_t>(y, x) = 0;
			else
				out.at<uint16_t>(y, x) = (uint16_t)(scale*tmp);
		}
	}
	imwrite(filename, out);
}
void SavePFMFile_prob(const cv::Mat& img, const char* filename)
{
	cv::Size szt = img.size();
	cv::Mat out(szt.height, szt.width, CV_16UC1);
	float scale = -10000;
	for (int y = 0; y < szt.height; ++y)
	{
		for (int x = 0; x < szt.width; ++x)
		{
			float tmp = img.at<float>(y, x);
			if (std::isnan(tmp))
				out.at<uint16_t>(y, x) = 0;
			else if(std::abs(tmp)>1)
				out.at<uint16_t>(y, x) = 1;
			else
				out.at<uint16_t>(y, x) = (uint16_t)(scale*tmp);
		}
	}
	imwrite(filename, out);
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
	string outfile = folderpath + "/" + "bayesian_" + outfilename + filename + ".png";
	string outfile_prob = folderpath + "/" + "prob_out_" + outfilename + filename + ".png";


	cv::Mat img_ao_mat, img_bo_mat, img_tmp;
	int rpyrtype, nochannels, incoltype;
	incoltype = CV_LOAD_IMAGE_GRAYSCALE;
	rpyrtype = CV_32FC1;
	nochannels = 1;
	img_ao_mat = cv::imread(imgfile_ao, incoltype);
	img_bo_mat = cv::imread(imgfile_bo, incoltype);
	
	int selectchannel = SELECTCHANNEL;
	ParamBDIS param;

	int acnt = 4;
	param.lv_f = atoi(argv[acnt++]);
	param.lv_l = atoi(argv[acnt++]);
	param.maxiter = atoi(argv[acnt++]);
	param.miniter = atoi(argv[acnt++]);
	param.mindprate = atof(argv[acnt++]);
	param.mindrrate = atof(argv[acnt++]);
	param.minimgerr = atof(argv[acnt++]);
	param.patchsz = atoi(argv[acnt++]);
	param.poverl = atof(argv[acnt++]);
	param.usefbcon = atoi(argv[acnt++]);
	param.patnorm = atoi(argv[acnt++]);
	param.costfct = atoi(argv[acnt++]);
	param.verbosity = atoi(argv[acnt++]);
	param.ratio_patch_valid = atof(argv[acnt++]);
	param.num_window = atoi(argv[acnt++]);
	param.unit_disturb = atof(argv[acnt++]);
	param.min_prob = atof(argv[acnt++]);
	//if (atof(argv[acnt++]) != 0)
	param.bool_var = true;

	cv::Mat flowout;
	cv::Mat probabilityout; 
	getDisp(img_ao_mat, img_bo_mat, selectchannel, rpyrtype, param, flowout, probabilityout);


	/*probabilityout = probabilityout(
			cv::Rect((int)floor((float)padw / 2.0f), (int)floor((float)padh / 2.0f), width_org, height_org));*/
	SavePFMFile(flowout, outfile.c_str());
	/*SavePFMFile_jingwei(flowout, probabilityout, param.min_prob, outfile.c_str());
	SavePFMFile_jingwei(probabilityout, probabilityout, param.min_prob, outfile_prob.c_str());*/
	if (param.verbosity > 1)
	{
		gettimeofday(&tv_end_all, NULL);
		double tt = (tv_end_all.tv_sec - tv_start_all.tv_sec) * 1000.0f +
				(tv_end_all.tv_usec - tv_start_all.tv_usec) / 1000.0f;
		printf("TIME (Saving flow file  ) (ms): %3g\n", tt);
	}
	if (param.bool_var)
	{
		SavePFMFile_prob(probabilityout, outfile_prob.c_str());
	}

	return 0;
}


    


