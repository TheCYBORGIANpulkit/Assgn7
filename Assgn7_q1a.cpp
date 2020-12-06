#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "Library.cpp"
using namespace std;

float f(float x, float y){
    float z = y*(log(y))/x;
    return z;
}

int main(){
    int N = 100;
    float xi = 2;
    float xf = 10;
    float y = 2.71828;
    Explicit_Euler(f,xi,10,y,200);
    cout << "Please find the attached file and the plot." << endl;
    return 0;
}
//"RESULT : Please find the attached csv file and the plot."
