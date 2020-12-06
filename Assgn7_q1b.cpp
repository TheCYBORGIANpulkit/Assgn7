#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "Library.cpp"
using namespace std;

float f(float x, float y){
    float z = 6 - 2*(y/x);
    return z;
}

int main(){
    int N = 100;
    float xi = 3;
    float xf = 12;
    float y = 1;
    Explicit_Euler(f,xi,12,y,200);
    cout << "Please find the attached file and the plot." << endl;
    return 0;
}
