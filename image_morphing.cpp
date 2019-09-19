#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <cstdio>
#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;

#define _USE_MATH_DEFINES
#define GLM_FORCE_RADIANS
void onMouse(int Event, int x, int y, int flags, void* param);

float aspect;
Point VertexLeftTop(-1, -1);
Point VertexRightDown(-1, -1);
vector<Point> ptarray_a;
vector<Point> ptarray_b;
vector<Point> ptarray_dest_lp;
vector<Point> ptarraytmp, transformed_a, transformed_b;
vector<float> linelength_a, linelength_b, linelength_dest;
vector<float> parmaa;
vector<float> parmab;
vector<float> parmba;
vector<float> parmbb, parmdb, parmda;

vector<vector<float>>weight_a, weight_b;

//滑鼠事件
void onMouse(int Event, int x, int y, int flags, void* param){
	if (Event == CV_EVENT_LBUTTONDOWN){
		ptarraytmp.push_back(cvPoint(x, y));
	}
}

//圖上畫線
void drawlineonpic(string win, Mat img, vector<Point> *ptarray){
	imshow(win, img);
	setMouseCallback(win, onMouse, NULL);
	while (true){
		for (int i = 0; (i + 1) < ptarraytmp.size(); i += 2){
			line(img, ptarraytmp[i], ptarraytmp[i + 1], Scalar(255, 0, 0), 2);
		}
		imshow(win, img);
		if (cvWaitKey(33) == 27){
			*ptarray = ptarraytmp;
			ptarraytmp.clear();
			break;
		}
	}
}

//計算長度
double getlength(Point a1, Point a2){
	return sqrt(pow(a1.x - a2.x, 2) + pow(a1.y - a2.y, 2));
}

//目標圖座標計算
void get_dest_line_point(vector<Point> pt_a, vector<Point> pt_b, float t){
	for (int i = 0; i < (pt_a).size(); i++){
		cout << "a= " << pt_a << " b= " << pt_b << endl;
		ptarray_dest_lp.push_back(Point(pt_a[i].x * t + pt_b[i].x * (1 - t), pt_a[i].y * t + pt_b[i].y * (1 - t)));
	}
	cout << "mid = " << ptarray_dest_lp << endl;
}
float compare(float a, float b, float c){
	if (a >= b)
		if (b >= c)
			return c;
		else
			return b;
	else
		if (a >= c)
			return c;
		else
			return a;
}


//權重計算
void trans(Point now, Point q, Point q_, Point p, Point p_, float *s_x, float *s_y, float *weight){
	Point x_p((now.x - p.x), (now.y - p.y)), q_p((q.x - p.x), (q.y - p.y));

	float q_p_l = getlength(q, p);
	float d = 0.1, c = 1, k = 1;
	float u = (x_p.x*q_p.x + x_p.y*q_p.y) / pow(q_p_l, 2);
	float v = (x_p.x*q_p.y - x_p.y*q_p.x) / (q_p_l);
	Point s, q_p_(q_.x - p_.x, q_.y - p_.y);
	float q_p_l_ = getlength(q_, p_);

	*s_x = p_.x + u*q_p_.x + v*q_p_.y / q_p_l_ - now.x;
	*s_y = p_.y + u*q_p_.y - v*q_p_.x / q_p_l_ - now.y;
	float dist;
	if (u < 0)
		dist = getlength(p, now);
	else if (u>1)
		dist = getlength(q, now);
	else
		dist = abs(v);
	*weight = pow(pow(q_p_l, k) / (d + dist), c);

}

//顏色分配
Vec3b redecide(float x, float y, Mat img){
	//cout << "a=" << a << " b= " << b << endl;
	if (x < 0)
		x = 0;
	if (x > img.cols - 2)
		x = img.cols - 2;
	if (y < 0)
		y = 0;
	if (y > img.rows - 2)
		y = img.rows - 2;
	Vec3b finalcolor;

	//cout << "a=" << a << " b= " << b << endl;
	Point a(floor(x), floor(y)), b(floor(x + 1), floor(y)), c(floor(x), floor(y + 1)), d(floor(x + 1), floor(y + 1));
	Vec3b color_a = img.at<Vec3b>(a.y, a.x), color_b = img.at<Vec3b>(b.y, b.x), color_c = img.at<Vec3b>(c.y, c.x), color_d = img.at<Vec3b>(d.y, d.x);

	Vec3b tmpa, tmpb;
	tmpa = (b.x - x)*color_a + (x - a.x)*color_b;
	tmpb = (d.x - x)*color_c + (x - c.x)*color_d;
	finalcolor = (d.y - y)*tmpa + (y - a.y)*tmpb;
	//cout << "aaaaaaaaaaaaaaa"<< endl;

	return finalcolor;
	/*
	finalcolor = img.at<Vec3b>(b, a);
	return finalcolor;*/
}
Point get_mid(Point a, Point b){
	Point mid;
	mid.x = (a.x + b.x) / 2;
	mid.y = (a.y + b.y) / 2;
	return mid;
}
int main()
{
	Mat imga = imread("a1.jpg", CV_LOAD_IMAGE_COLOR);
	Mat imgb = imread("a2.jpg", CV_LOAD_IMAGE_COLOR);
	//Mat imga = imread("a.png", CV_LOAD_IMAGE_COLOR);
	//Mat imgb = imread("b.png", CV_LOAD_IMAGE_COLOR);
	Mat showa, showb;
	imga.copyTo(showa);
	imgb.copyTo(showb);
	Mat dest_from_a, dest_from_b;
	imga.copyTo(dest_from_a);
	imgb.copyTo(dest_from_b);
	//imshow("ddest_from_a", dest_from_a);
	//imshow("dest_from_b", dest_from_b);
	int countline;
	drawlineonpic("image a", showa, &ptarray_a);
	drawlineonpic("image b", showb, &ptarray_b);
	float t = 0.5;

	if (ptarray_a.size() == ptarray_b.size() && ptarray_a.size()>0){
		countline = ptarray_a.size() / 2;
		get_dest_line_point(ptarray_a, ptarray_b, t);
		cout << "dest = " << ptarray_dest_lp << endl;

		/*
		for (int i = 0; (i + 1) < ptarray_a.size(); i += 2){
		Get_k_and_b(ptarray_dest_lp[i], ptarray_dest_lp[i + 1], &parmda, &parmdb);
		linelength_dest.push_back(getlength(ptarray_dest_lp[i], ptarray_dest_lp[i + 1]));
		}*/

		for (int i = 0; i < imga.cols; i++){ //cols=x     rows=y
			for (int j = 0; j < imga.rows; j++){
				int countp = 0;
				Point now(i, j), source_a(0, 0), source_b(0, 0);
				float source_ax = 0, source_ay = 0, source_bx = 0, source_by = 0;
				float a_x = 0, a_y = 0, b_x = 0, b_y = 0;
				float wsum_a = 0, wsum_b = 0;
				for (int k = 0; k < ptarray_dest_lp.size(); k += 2){
					float nowweight_a, nowweight_b;
					float sa_x, sa_y, sb_x, sb_y;
					trans(now, ptarray_dest_lp[k + 1], ptarray_a[k + 1], ptarray_dest_lp[k], ptarray_a[k], &sa_x, &sa_y, &nowweight_a);
					trans(now, ptarray_dest_lp[k + 1], ptarray_b[k + 1], ptarray_dest_lp[k], ptarray_b[k], &sb_x, &sb_y, &nowweight_b);

					source_ax += sa_x*nowweight_a;
					source_ay += sa_y*nowweight_a;
					source_bx += sb_x*nowweight_b;
					source_by += sb_y*nowweight_b;

					wsum_a += nowweight_a;
					wsum_b += nowweight_b;
					countp++;
				}
				a_x = now.x + source_ax / wsum_a;
				a_y = now.y + source_ay / wsum_a;
				b_x = now.x + source_bx / wsum_b;
				b_y = now.y + source_by / wsum_b;

				Vec3b intensity_a = redecide(a_x, a_y, imga);
				Vec3b intensity_b = redecide(b_x, b_y, imgb);

				dest_from_a.at<Vec3b>(j, i) = intensity_a;
				dest_from_b.at<Vec3b>(j, i) = intensity_b;

			}
			cout << "i=" << i << endl;
		}

		imshow("ddest_from_a", dest_from_a);

		imshow("dest_from_b", dest_from_b);
		Mat dst;
		addWeighted(dest_from_b, 0.5, dest_from_a, 0.5, 0, dst);
		imshow("dst", dst);
		cvWaitKey(0);
	}
	else{
		cout << "no" << "\n";
		return 0;
	}

	return 0;
}

