// GaussProjection.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <iostream>
#include <string>
#include<iomanip>
const double pi = acos(-1);
using namespace std;
const double rou = 180.00 * 3600.00 / pi;

//首先准备椭球参数，有克氏椭球和1975国际标准椭球两套参数
struct KrasEllipsoid
{
	double a = 6378245;
	double b = 6356863.0187730473;
	double c = 6399698.9017827110;
	double f = 1 / 298.3;
	double e2 = 0.006693421622966;
	double SecondE2 = 0.006738525414683;
}Kras;
struct seventy_five_Ellipsoid
{
	double a = 6378140;
	double b = 6356755.2881575287;
	double c = 6399596.6519880105;
	double f = 1 / 298.257;
	double e2 = 0.006694384999588;
	double SecondE2 = 0.006739501819473;
}seventy_five;

//度分秒化为秒
double DFM2Deg(double du,double fen,double miao)
{
	double result = du + (fen / 60.00) + (miao / 3600.00);
	return result;
}
//秒化为弧度
double Deg2Rad(double deg)
{
	double rad = deg * pi / 180.00;
	return rad;
}
//度分秒化为弧度
double DFM2Rad(double du, double fen, double miao)
{
	double result = DFM2Deg(du, fen, miao);
	double rad = Deg2Rad(result);
	return rad;
}

//高斯投影正算函数
void Gauss0(double B, double L, double k, string name)
{
	    double x=0, y=0;
		double L0,l;
		L0 = 6 * (int(L / 6)) - 3;
		l = L - L0;
		l = Deg2Rad(l);
		double sin2B = sin(B)*sin(B);
		double cos2B = cos(B)*cos(B);
		double t = tan(B);
		//下面计算自赤道起的子午线弧长X和计算用到的各系数
		double N, a0, a3, a4, a5, a6;
		if (k == 0)//根据选择的椭球来进行参数的选择
		{
			N = 6399698.902 - (21562.267 - (108.973 - 0.612*cos2B)*cos2B)*cos2B;
			a0 = 32140.404 - (135.3302 - (0.7092 - 0.004*cos2B)*cos2B)*cos2B;
			a3 = (0.3333333 + 0.001123*cos2B)*cos2B - 0.1666667;
			a4 = (0.25 + 0.00252*cos2B)*cos2B - 0.04166;
			a5 = 0.0083 - (0.1667 - (0.1968 + 0.004*cos2B)*cos2B)*cos2B;
			a6 = (0.166*cos2B - 0.084)*cos2B;
			x = 6367558.4969*B - (a0 - (0.5 + (a4 + a6 * l*l)*l*l)*l*l*N)*sin(B)*cos(B);
			y = (1 + (a3 + a5 * l*l)*l*l)*l*N*cos(B);
		}
		if (k == 1)
		{
			N = 6399596.652 - (21565.045 - (108.996 - 0.603*cos2B)*cos2B)*cos2B;
			a0 = 32144.5189 - (135.3646 - (0.7034 - 0.0041*cos2B)*cos2B)*cos2B;
			a3 = (0.3333333 + 0.001123*cos2B)*cos2B - 0.1666667;
			a4 = (0.25 + 0.00253*cos2B)*cos2B - 0.04167;
			a5 = 0.00878 - (0.1702-0.20382*cos2B)*cos2B;
			a6 = (0.167*cos2B - 0.083)*cos2B;
			x = 6367452.1328*B - (a0 - (0.5 + (a4 + a6 * l*l)*l*l)*l*l*N)*sin(B)*cos(B);
			y = (1 + (a3 + a5 * l*l)*l*l)*l*N*cos(B);
		}
		cout << "下面输出计算结果：" << endl;
		cout << "点名：" << name << endl;
		cout << "x=" << setprecision(10)<< x << "m," << "y=" <<setprecision(10)<< y << "m" << endl;
		system("pause");
}

//高斯反算函数
void Gauss1(double x, double y, double k, string name)
{
	double Bf, Bf1, beta, Z, Nf, b2, b3, b4, b5,B,l;
	if (k == 0)//根据选择的椭球进行参数计算
	{
		beta = x / 6367558.4969;
		Bf1 = beta*rou + (50221746 + (293622 + (2350 + 22 * cos(beta)*cos(beta))*cos(beta)*cos(beta))*cos(beta)*cos(beta))*sin(beta)*cos(beta)*rou*1e-10;
		Bf = Bf1/rou;//由于参与三角函数的运算必须是弧度制，因此Bf1是弧度制变量，Bf则是以秒为单位，便于进行输出
		Nf = 6399698.902 - (21562.267 - (108.973 - 0.612*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf);
		Z = y / (Nf*cos(Bf));
		b2 = (0.5 + 0.003369*cos(Bf)*cos(Bf))*sin(Bf)*cos(Bf);
		b3 = 0.333333 - (0.166667 - 0.001123*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf);
		b4 = 0.25 + (0.16161 + 0.00562*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf);
		b5 = 0.2 - (0.1667 - 0.0088*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf);
		B = Bf1 - (1 - (b4 - 0.12*Z*Z)*Z*Z)*Z*Z*b2*rou;
		l = (1 - (b3 - b5 * Z*Z)*Z*Z)*Z*rou;
	}
	if (k == 1)
	{
		beta = x  / 6367452.133;
		Bf1 = beta*rou + (50228976 + (293697 + (2383 + 22 * cos(beta)*cos(beta))*cos(beta)*cos(beta))*cos(beta)*cos(beta))*sin(beta)*cos(beta)*rou*1e-10;
		Bf = Bf1 / rou;
		Nf = 6399596.652 - (21565.047 - (109.003 - 0.612*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf);
		Z = y / (Nf*cos(Bf));
		b2 = (0.5 + 0.00336975*cos(Bf)*cos(Bf))*sin(Bf)*cos(Bf);
		b3 = 0.333333 - (0.166667 - 0.001123*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf);
		b4 = 0.25 + (0.161612 + 0.005617*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf);
		b5 = 0.2 - (0.1667 - 0.00878*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf);
		B = Bf1 - (1 - (b4 - 0.147*Z*Z)*Z*Z)*Z*Z*b2*rou;
		l = (1 - (b3 - b5 * Z*Z)*Z*Z)*Z*rou;
	}
	cout << "下面输出计算结果：" << endl;
	cout << "点名：" << name << endl;
	cout << "B=" <<setprecision(10)<< B << "秒，" << "l=" <<setprecision(10)<< l << "秒" << endl;
	system("pause");
}

//下面是主函数
int main()
{
	double B, L, x, y;
	double du, fen, miao;
	int i,k;//这里需要两个变量进行分支选择，一个确定换算所用的椭球参数，一个选择是正算还是反算，以便调用对应的函数
	string name;
	cout << "请选择高斯投影使用的椭球，克氏椭球为0,1975国际椭球为1"<<endl;
	cin >> k;
	cout << "请选择高斯投影的计算方式，正算为0，反算为1" << endl;
	cin >> i;
	if (i == 0)
	{
		cout << "请输入点名" << endl;
		cin >> name;
		cout << "请输入大地纬度坐标B，形式为度-分-秒，如：30 22 21" << endl;
		cin >> du >> fen >> miao;
		B = DFM2Rad(du, fen, miao);
		cout << "请输入大地经度坐标L，形式为度-分-秒，如：35 11 23" << endl;
		cin >> du >> fen >> miao;
		L = DFM2Deg(du, fen, miao);
		Gauss0(B, L, k, name);
	}	
	if (i == 1)
	{
		cout << "请输入点名" << endl;
		cin >> name;
		cout << "请输入国家统一坐标x，单位为m" << endl;
		cin >> x;
		cout << "请输入国家统一坐标y，单位为m" << endl;
		cin >> y;
		Gauss1(x, y, k, name);
	}
	return 0;
}

