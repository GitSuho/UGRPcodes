#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;
int Norm(double _x, double _y){
    return sqrt(_x*_x+_y*_y);
}

////////////////////////////initial values/////////////////
double m_val[3] = {1, 1, 1};
double G = 9.8; double dt = 0.01; 
double initial_val[3][4] = {{1, 0, 2, 1 },
                            {0, 1, 1, 2},
                            {0.5, 0.5, 0.5, 0.5}};
long double E = (0.5)*m_val[0]*(initial_val[0][2]*initial_val[0][2] + initial_val[0][3]*initial_val[0][3])
            - G*m_val[0]*m_val[1]/Norm(initial_val[0][0]-initial_val[1][0], initial_val[0][1]-initial_val[1][1])
              + (0.5)*m_val[1]*(initial_val[1][2]*initial_val[1][2] + initial_val[1][3]*initial_val[1][3])
            - G*m_val[1]*m_val[2]/Norm(initial_val[1][0]-initial_val[2][0], initial_val[1][1]-initial_val[2][1])
              + (0.5)*m_val[2]*(initial_val[2][2]*initial_val[2][2] + initial_val[2][3]*initial_val[2][3])
            - G*m_val[2]*m_val[0]/Norm(initial_val[2][0]-initial_val[0][0], initial_val[2][1]-initial_val[0][1]);
///////////////////////////////////////////////////////////////////


long double T( long double** _xv_array){
    long double T = E;
    for (int i = 0 ; i < 3 ; i++){
        int j = (i+1)%3;
        T += G*m_val[i]*m_val[j]/Norm(_xv_array[i][0]-_xv_array[j][0], 
                                      _xv_array[i][1]-_xv_array[j][1]);
    } 
    return T;
}
long double T_diff(long double** _xv_array){
    long double T_diff = 0;
    for (int i = 0; i < 3; i++){
        int j = (i+1)%3;
        T_diff += G*m_val[i]*m_val[j]
        *((_xv_array[i][0]-_xv_array[j][0])*(_xv_array[i][2]-_xv_array[j][2])
         +(_xv_array[i][2]-_xv_array[j][2])*(_xv_array[i][3]-_xv_array[j][3]))
        /((_xv_array[i][0]-_xv_array[j][0])*(_xv_array[i][0]-_xv_array[j][0])
         +(_xv_array[i][1]-_xv_array[j][1])*(_xv_array[i][1]-_xv_array[j][1]));
    }
    return T_diff;
}
long double diff_T_xv(int i, int j, long double** _xv_array){
    long double result = 0;
    for ( int n = 0 ; n < 2 ; n++ ){
        int k = (i+n)%3; int l = (k+1)%3;
        result += pow(-1, n+2)*G*m_val[k]*m_val[l]*2*(_xv_array[k][j]-_xv_array[l][j])
        /(sqrt((_xv_array[i][0]-_xv_array[j][0])*(_xv_array[i][0]-_xv_array[j][0])
              +(_xv_array[i][1]-_xv_array[j][1])*(_xv_array[i][1]-_xv_array[j][1]))
             *((_xv_array[i][0]-_xv_array[j][0])*(_xv_array[i][0]-_xv_array[j][0])
              +(_xv_array[i][1]-_xv_array[j][1])*(_xv_array[i][1]-_xv_array[j][1])));
    }
    return result;
}
long double** F( long double** _xv_array){
    long double** F;
    long double dT;
    for (int i = 0; i < 3 ; i++){
        for( int j = 0 ; j < 2 ; j ++ ){
            F[i][j] = T_diff(_xv_array)/T(_xv_array)*_xv_array[i][j+2];
            for(int k = 0 ; k < 3; k++){
                for(int l = 0; l < 2; l++){
                    dT = diff_T_xv(i, j, _xv_array);
                    F[i][j] += - 1/T(_xv_array)*dT*_xv_array[k][l+2]
                               + 1/(2*T(_xv_array))*dT*(_xv_array[k][l+2]*_xv_array[k][l+2]);
                }}}}
    return F;
}
long double** geod( long double** _X){
    long double** result;
    for (int i = 0 ; i < 3; i++){
        result[i][0] = _X[i][2]; result[i][1] = _X[i][3];
        result[i][3] = F(_X)[i][0] ; result[i][4] = F(_X)[i][1];
    }
    return result;
}
long double** return_ptr34( long double coeff, long double** ptr){
    long double** result;
    for (int i = 0 ; i < 3 ; i ++){
        for (int j = 0 ; j < 4 ; j++ ){
            result[i][j] = coeff * ptr[i][j];
        }}
    return result;
}
long double** sum_two_ptr(long double** ptr1, long double** ptr2){
    long double** result;
    for (int i = 0 ; i < 3 ; i ++){
        for (int j = 0 ; j < 4 ; j++ ){
            result[i][j] = ptr1[i][j] + ptr2[i][j];
        }}
    return result;
}

int main(){
    string filename = "123";
    ofstream writefile(filename + ".txt");


    long double** X;
    for (int i = 0 ; i < 3 ; i ++){
        for (int j = 0; j < 4 ; j++){
            X[i][j] = initial_val[i][j];
        }
    }

    string write_contents;
    long double** k1;long double** k2;long double** k3;long double** k4;
    writefile << "\n";
    for (double t = 0.0; t < 5; t += dt){
        k1 = geod(X);
        k2 = geod( sum_two_ptr(X, return_ptr34(dt/2, k1)));
        k3 = geod( sum_two_ptr(X, return_ptr34(dt/2, k2)));
        k4 = geod( sum_two_ptr(X, return_ptr34(dt  , k3)));

        
        write_contents = "|";
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 4; j++){
                X[i][j] = X[i][j] + (dt / 6) * (k1[i][j] + 2 * k2[i][j] + 3 * k3[i][j] + 4 * k4[i][j]);
                write_contents += to_string(X[i][j]) + "|";
            }}
        
        writefile << write_contents << endl;
    }
    writefile.close();
    return 0;
}