#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "Library.cpp"
using namespace std;

float funcRK4(float u, float v, float x, int n, int m)
{
	if (m == 1){
		if (n == 1)
			return v;
		if (n == 2)
			return (1 - x - v);
	}
	else{
		if (n == 1)
			return v;
		else
			return (1 + v);
	}
}


void cases(float h)
{
	float u[1], v[1], x, t;
	ofstream outfile1;
	ofstream outfile2;
    if(h == 0.5){
		outfile1.open ("RK4_1a.csv");
		outfile2.open ("RK4_1b.csv");
	}
	else if(h == 0.2){
		outfile1.open ("RK4_2a.csv");
		outfile2.open ("RK4_2b.csv");
	}
    else if(h == 0.02){
		outfile1.open ("RK4_3a.csv");
		outfile2.open ("RK4_3b.csv");
    }
	else if(h == 0.05){
		outfile1.open ("RK4_4a.csv");
		outfile2.open ("RK4_4b.csv");
	}

	//for the range [0,5]
	x = 0;
	u[0] = 2;
	v[0] = 1;
	t = 5;
	while (fabs(x) < fabs(t))
	{
		RK4(u, v, x, t, h, 1);
		outfile1 << t << "," << x << "," << u[0] << endl;
		x += h;
	}
	outfile1.close();

	//for the range [-5,0]
	x = 0;
	u[0] = 2;
	v[0] = 1;
	t = -5;
	float p = -h;
	while (fabs(x) < fabs(t))
	{
		RK4(u, v, x, t, p, 1);
		outfile2 << t << "," << x << "," << u[0] << endl;
		x += p;
	}
	outfile2.close();
	cout << "ODE solved using Runge Kutta method-4th order with h = " << h << endl;
}

int main()
{
	cout << "Given ODE is : y'' + y' = 1 - x" << endl;
	cout << "Initial boundary value: y(0) = 2, y'(0) = 1" << endl;


	cases(0.2);
	cases(0.5);
	cases(0.05);
	cases(0.02);

	cout << "Please find the attached csv file and the plots" << endl;
	return 0;
}
/*
RESULT: PLease find the attached csv file and the plot
NOTE: The plots have been generated for "h = 0.2, 0.5,0.02,0.05" and the csv file also"



