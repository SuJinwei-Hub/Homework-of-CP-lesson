/*
	利用Markov链+Metropolis更新方案做随机抽样，用蒙特卡洛方法求六维积分，并求不同抽样频次对应的协方差。
	By 苏锦伟18352035
	2021.05.16
*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include <string.h>
#include<math.h>
#include <windows.h>
#define PI acos(-1)
#define NI 6			// 积分维度


int main() {
	double x0 = 1.0;
	int delta = 1;
	double delta_n;
	double w = 0.0, r = 0.0;
	int N_balance = 10000;
	int N_max = 50000;
	int N0 = 50;
	int n = 0;
	int N_num = (N_max - N_balance) / N0;
	int npx = 0;
	//定义一个二维数组，横坐标为积分变量，每一行为抽出的随机点
	double** px_n = new double* [NI];
	for (int i = 0; i < NI; i++) px_n[i] = new double[N_num];
	double* x_ = new double[NI];
	double* xn = new double[NI];

	//
	double** px_nn = new double* [NI];
	for (int i = 0; i < NI; i++) px_nn[i] = new double[N_max - N_balance];
	void cov(double** px, int row);
	//

	void txt_save(int m, double* xn, const char* name);
	double func_px(double* x);
	double func_gx(double* x);
	//生成满足p(x)分布的取样点
	int m = 0, mMax = 1000;
	double* i_resm = new double[mMax];
	while (m < mMax){
		n = 0;
		srand(clock());		//此处不要用time(NULL)，因为它单位是秒，一秒内随机数种子是相同的
		for (int i = 0; i < NI; i++) xn[i] = x0;		// 初值
		while (n < N_max)			//Metropolis随机游走
		{
			for (int i = 0; i < NI; i++) {
				delta_n = (double(rand() % (2000 * delta))) / 1000.0 - delta;
				x_[i] = xn[i] + delta_n;
			}
			w = func_px(x_) / func_px(xn);
			if (w >= 1) {
				for (int i = 0; i < NI; i++) {
					xn[i] = x_[i];
				}
			}
			else if (w < 1) {
				r = double(rand() % 1000 / 1000.0);
				if (r <= w) {
					for (int i = 0; i < NI; i++) {
						xn[i] = x_[i];
					}
				}
				else if (r > w) {
					for (int i = 0; i < NI; i++) {
						xn[i] = xn[i];
					}
				}
			}
			if (n >= N_balance && (n % N0 == 0)) {
				npx = (n - N_balance) / N0;
				for (int i = 0; i < NI; i++) {
					px_n[i][npx] = xn[i];
				}
			}
			//用于计算sigma和协方差，在单独进行积分时可注释掉
			if (n >= N_balance) {
				for (int i = 0; i < NI; i++) {
					px_nn[i][n - N_balance] = xn[i];
				}
			}
			n += 1;
		}
		//积分
		double* xl = new double[NI];
		double i_res = 0.0, i_res2 = 0.0;
		double sigma = 0.0;
		for (int j = 0; j < N_num; j++) {
			for (int i = 0; i < NI; i++) xl[i] = px_n[i][j];
			i_res += func_gx(xl);
			//i_res2 += pow(func_gx(xl), 2);
		}
		//sigma = sqrt((i_res2 / double(N_num) - pow(i_res / double(N_num), 2)) / (double(N_num) - 1));
		i_res = i_res / double(N_num);
		//printf("The result is %lf\n", i_res);
		i_resm[m] = i_res;
		m++;
	}
	//计算协方差
	cov(px_nn, N_max - N_balance);
	//保存多次积分结果
	txt_save(mMax, i_resm, "result.txt");
	delete[]px_n;
	delete[]px_nn;
	delete[]i_resm;
}


double func_px(double* x) {		//抽样函数
	return pow(2 * PI / sqrt(3), -3) * exp(-(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + x[4] * x[4] + x[5] * x[5]) + (x[0] * x[3] + x[1] * x[4] + x[2] * x[5]));
}


double func_gx(double* x) {		//被积函数
	return pow(2 * PI / sqrt(3), 3) * exp(-(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + x[4] * x[4] + x[5] * x[5]) * 0.5);
}



void cov(double** px, int row) {	//计算协方差
	double xn, xn2, xxl;
	double covl = 0.0;
	int lm = 200, l0 = 10, k = 0, ln = 0;
	double* cov = new double[21];
	void txt_save(int m, double* xn, const char* name);
	for (int l = 0; l <= lm; l += l0) {
		covl = 0.0;
		for (int i = 0; i < NI; i++) {
			xn = xn2 = xxl = 0.0;
			ln = 0;
			for (int j = 0; j < row; j++) {
				xn += px[i][j];
				xn2 += (px[i][j] * px[i][j]);
				if (j + l > row) break;
				xxl += px[i][j] * px[i][j + l];
				ln++;
			}
			xxl = xxl / double(ln);
			xn = xn / double(row);
			xn2 = xn2 / double(row);
			covl += fabs((xxl - xn * xn) / (xn2 - xn * xn));
		}
		cov[k] = covl / NI;
		k++;
	}
	txt_save(21, cov, "cov.txt");
	delete[]cov;
}


void txt_save(int m, double* xn, const char* name) {		//将数组写入txt
	FILE* fp;
	errno_t err;
	err = fopen_s(&fp, name, "w");
	for (int i = 0; i < m; i++) {
		fprintf(fp, "%lf\n", xn[i]);
	}
	//fwrite(xn, sizeof(double), m, fp);
	fclose(fp);
}