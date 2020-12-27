#include <iostream>  //declaring variables
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
using namespace std;

void RK4(float u[1], float v[1], float x, float t, float h, int p)
{
	float k11, k12, k13, k14;
	float k21, k22, k23, k24;
	k11 = h * funcRK4(u[0], v[0], x, 1, p);
	k21 = h * funcRK4(u[0], v[0], x, 2, p);
	k12 = h * funcRK4(u[0] + 0.5 * k11, v[0] + 0.5 * k21, x + 0.5 * h, 1, p);
	k22 = h * funcRK4(u[0] + 0.5 * k11, v[0] + 0.5 * k21, x + 0.5 * h, 2, p);
	k13 = h * funcRK4(u[0] + 0.5 * k12, v[0] + 0.5 * k22, x + 0.5 * h, 1, p);
	k23 = h * funcRK4(u[0] + 0.5 * k12, v[0] + 0.5 * k22, x + 0.5 * h, 2, p);
	k14 = h * funcRK4(u[0] + k13, v[0] + k23, x + h, 1, p);
	k24 = h * funcRK4(u[0] + k13, v[0] + k23, x + h, 2, p);
	u[0] += (k11 + 2 * k12 + 2 * k13 + k14) / 6;
	v[0] += (k21 + 2 * k22 + 2 * k23 + k24) / 6;
}

void Explicit_Euler(float(*f)(float,float),float xi,float xf,float y, int N){
    float h = (xf - xi)/N;
    float x = xi;
    ofstream outfile;
    outfile.open("explicit1.csv");
    for(int i = 0;i<N;i++){
        outfile << x << "," << y << endl;
        float K1 = h*((*f)(x + h,y));
        y = y + K1;
        x = x + h;
    }
    outfile.close();
}

double Simpsons(double (*fn)(float), float a, float b, int N){
    double h = (b - a)/N;
    double s = 0;
    for(int i = 0;i<N;i++){
        double x = (h/3)*((*fn)(a + i*h));
        if(i == 0 || i == N)s = s + x;
        else if(i%2 == 0)s = s + 2*x;
        else s = s + 4*x;
    }
    return s;
}

double Trapezoidal(double (*fn)(float), float a, float b, int N){
    double h = (b - a)/N;
    double s = 0;
    for(int i = 0;i<N;i++){
        double x = (h/2)*((*fn)(a + i*h));
        if(i == 0 || i == N)s = s + x;
        else s = s+2*x;
    }
    return s;
}

double Midpoint(double (*fn)(float),float a, float b, int N ){
    double h = (b - a)/N;
    double s = 0;
    for(int i = 0;i<N;i++){
        double x = (a+i*h + a + (i+1)*h)/2;
        s = s + h*(*fn)(x);
        //cout << s << endl;
    }
    return s;
}

void Bracketing(float (*fn)(float),double a, double b){
    int Count = 0;
    double beta = 0.75;
    for(int j = 0;j< 200;j++){
        if((*fn)(a)*(*fn)(b) > 0){
            if(abs((*fn)(a)) > abs((*fn)(b))){
                b = b + beta*(b - a);
            }
            else a = a - beta*(b - a);
        }
    }
}
double Bisection(float (*fn)(float), double a, double b){
    Bracketing((*fn),a,b);
    int Count = 0;
    double arrC[200];
    //cout << "c:" << c << endl;
    for(int j = 0;j<200;j++){
        double c = (a+b)/2;
        if(abs(a - b) > 0.000001){
            arrC[j] = abs(a - b);
            if((*fn)(a)*(*fn)(c) < 0) b = c;
            else  a = c;
            //Count ++;
        }
        else{
            cout << "after bisection" << endl ;
            cout << "c:" << c << endl;
            cout << "f(c):" << (*fn)(c) << endl;
            arrC[j] = abs(a - b);
            ofstream outfile;
            outfile.open("bisection.csv");
            for(int i = 0;i<= j; i++){
                outfile<< arrC[i] << endl;
            }
            outfile.close();
            cout << "Please find the attached csv file generated and the plot" << endl;
            return (*fn)(c);
        }
    }
}

double False_position(float (*fn)(float),float a, float b){
    Bracketing((*fn),a,b);
    double arrF[200];
    //cout << "c:" << c << endl;
    for(int j = 0;j<200;j++){
        float c = b - (((b-a)*((*fn)(b)))/((*fn)(b) - (*fn)(a)));
        if(abs((*fn)(c)) > 0.000001){
            arrF[j] = abs((*fn)(c));
            if((*fn)(a)*(*fn)(c) < 0) b = c;
            else a = c;
            //Count ++;
        }
        else{
            cout << "after false-positioning" << endl;
            cout << "c:" << c << endl;
            cout << "f(c):" << (*fn)(c) << endl;
            arrF[j] = abs((*fn)(c));
            ofstream outfile;
            outfile.open("falsi.csv");
            for(int i =0;i<j;i++){
                outfile<< arrF[j] << endl;
            }
            outfile.close();
            cout << "Please find the attached csv file generated and the plot" << endl;
            return (*fn)(c);
        }
    }
}

float Newton_Raphson(float (*fn)(float),float x0){
    float h = 0.001;
    float x = x0;
    float arrX[200];
    for(int j = 0;j<200;j++){
        float delf = ((*fn)(x + h) - (*fn)(x - h))/(2*h);
        float a = (*fn)(x)/delf;
        if(abs(a) > 0.0001){
            x = x - a;
            arrX[j] = abs(a);
        }
        else{
            cout << "after newton-raphson" << endl;
            cout << "root:" << x << endl;
            cout << "f(x): " << (*fn)(x) << endl;
            arrX[j] = abs(a);
            ofstream outfile;
            outfile.open("newton_raphson.csv");
            for(int i = 0;i<= j; i++){
                outfile<< arrX[i] << endl;
            }
            outfile.close();
            cout << "Please find the attached csv file generated and the plot" << endl;
            return (*fn)(x);
        }
    }
}

void Matrix_Product(float arrM[4][4], float arrN[4][4],int n){
    float arrP[n][n], x;
    //int arrQ[3][3];
    for(int i = 0;i<n;i++){
        for(int j = 0;j<n;j++){
            arrP[i][j] = 0;
            //arrQ[i][j] = 0;
        }
    }
    for(int a = 0;a<n;a++){
        for(int b = 0;b<n;b++){
            for(int i = 0;i<n;i++){
                x = float(arrM[a][i])*float(arrN[i][b]);
                //cout<<x<<endl;
                arrP[a][b] = arrP[a][b] + x;
                //cout<< arrP[a][b];
            }
            //cout<< arrP[a][b]<< " ";
            //cout<< endl;
        }
    }
    cout<<"The product A*A^(-1) matrix is:"<< endl;
    for(int i = 0;i<n;i++){
        for(int j = 0;j<n;j++){
            cout<<arrP[i][j]<< " ";
        }
        cout<<endl;
    }
}

void Part_Pivoting(float arrA[4][4]){
    for(int r = 0;r<3;r++){
        float pivot = arrA[r][r];
        //cout << pivot << endl;
        if(pivot == 0){
           for(int i = r+1;i<4;i++){
                if(abs(int(arrA[i][r])) > abs(int(pivot))){
                    float temp[4];
                    int j = 0;
                    for(int j = 0;j<4;j++){
                        temp[j] = arrA[r][j];
                        arrA[r][j] = arrA[i][j];
                        arrA[i][j] = temp[j];
                    }
                }
                else continue;
           }
        }
        else continue;
    }
}

void LU_decomposition(float arrA[4][4]){
    //Partial pivoting at first
    Part_Pivoting(arrA);
    for(int j = 0;j<4;j++){
        if(j == 0){
            //u[0][0] = arrA[0][0];
            //cout<< arrA[0][0] << endl;
            for(int i = 1;i<4;i++){
                arrA[i][j] = arrA[i][j]/arrA[0][0];
                //cout<< arrA[i][j] << endl;
            }
        }
        else{
            for(int i =0;i<4;i++){
                float s = 0;
                float a = 0;
                for(int k = 0;k<i;k++){
                    a = arrA[i][k]*arrA[k][j];
                    s = s + a;
                }
                //cout<< s << endl;
                if(i < j)arrA[i][j] = arrA[i][j] - s;
                else if(i == j){
                    arrA[i][i] = arrA[i][j] - s;
                    //if(i == j)cout<< arrA[i][j] << endl;
                }
                else{
                    float s = 0;
                    float a = 0;
                    for(int k = 0;k<j;k++){
                        a = arrA[i][k]*arrA[k][j];
                        s = s + a;
                    }
                    arrA[i][j] = (arrA[i][j] - s)/arrA[j][j];
                }
            }
        }
    }
}
float Determinant(float arrA[4][4]){
    //Pivoting(arrA);
    LU_decomposition(arrA);
    float p = 1;
    for(int i =0;i<4;i++){
        p = p*arrA[i][i];
    }
    return p;
}
// NOTE: This code assumes that the matrix A is already LU decomposed.
void FB_Substitution(float arrA[4][4], int n){
    //RHS matrix
    cout<<"The RHS matrix b is:"<< endl;
    float arrb[4][n];
    if(n == 1){
        ifstream myfileb;
        myfileb.open("A4_Q1b.txt");
        for(int i = 0;i<4;i++){
            myfileb >> arrb[i][0];
            cout<< arrb[i][0] <<" ";
            cout<<endl;
        }
    }
    else if(n == 4){
        ifstream myfileB;
        myfileB.open("A4_Q2b.txt");
        for(int i = 0;i<4;i++){
            for(int k =0;k<n;k++){
                myfileB >> arrb[i][k];
                cout<< arrb[i][k] <<" ";
            }
            cout<<endl;
        }
    }
    //Making augmented matrix
    int m = 4 + n;
    float arrL[4][m];
    for(int i = 0;i<4;i++){
        for(int j = 0;j<4;j++){
            if(i > j) arrL[i][j] = arrA[i][j];
            if(i == j) arrL[i][j] = 1;
            if(i < j) arrL[i][j] = 0;
        }
    }
    float arry[4][n];
    //solving Y
    if(n == 1){
        for(int i = 0;i<4;i++){
            float s = 0;
            for(int j = 0;j<i;j++){
                s = s + arrL[i][j]*arry[j][0];
            }
            arry[i][0] = (arrb[i][0] - s)/arrL[i][i];
        }
        //for backward substitution
        float arrU[4][4];
        for(int i = 0;i<4;i++){
            for(int j = 0;j<4;j++){
                if(i <= j) arrU[i][j] = arrA[i][j];
                if(i > j) arrU[i][j] = 0;
            }
        }
        //solving X
        float arrX[4][1];
        for(int i = 3;i >= 0;i--){
            float s = 0;
            for(int j = i+1;j<4;j++){
                s = s + arrU[i][j]*arrX[j][0];
            }
            //cout<< s << endl;
            arrX[i][0] = (arry[i][0] - s)/arrU[i][i];
        }
        cout<< "The solution is:" << endl;
        for(int i = 0;i<4;i++){
            cout << "x[" << i << "] = "<< arrX[i][0] << endl;
        }
    }

    else{
    //solving Y
        for(int i = 0;i<4;i++){
            for(int j = 0;j<n;j++){
                float s = 0;
                for(int k = 0;k<i;k++){
                    s = s + arrL[i][k]*arry[k][j];
                }
                arry[i][j] = (arrb[i][j] - s)/arrL[i][i];
                //cout<< arry[i][j] << " ";
            }
            //cout<< endl;
        }
        //for backward substitution
        float arrU[4][4];
        for(int i = 0;i<4;i++){
            for(int j = 0;j<4;j++){
                if(i < j)arrU[i][j] = arrA[i][j];
                else if(i == j)arrU[i][j] = arrA[i][j];
                else arrU[i][j] = 0;
            }
        }
        //solving X[4][4]
        float arrX[4][4];
        for(int i = 3;i >= 0;i--){
           for(int j = 0;j<n;j++){
                float s = 0;
                for(int k = i+1;k<n;k++){
                    s = s + arrU[i][k]*arrX[k][j];
                }
                arrX[i][j] = (arry[i][j] - s)/arrU[i][i];
            }
        }
        cout<< "The INVERSE is:" << endl;
        for(int i = 0;i<4;i++){
            for(int j =0; j< 4;j++){
                cout << arrX[i][j] << " ";
            }
            cout<< endl;
        }
    }
}


