#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;
long double Norm(long double _x, long double _y){
    return sqrtl(_x*_x+_y*_y);
}
string write_contents(long double _X[][4], double t){
    string result = to_string(t) + ":";
    for(int i = 0 ; i < 3; i ++){for (int j = 0; j < 4; j++){
        result += to_string(_X[i][j]) + ","; 
    }}
    return result;
}


double m_val[3] = {1, 1, 1};
double G = 9.8; double dt = 0.00001; 
double initial_val[3][4] = {{-1.0/2.0, -1.0/3.0, 1.0/3.0, -2.0/3.0},
                            {1.0/2.0, -1.0/3.0, -1.0/3.0, 2.0/3.0},
                            {0, 2.0/3.0, -1, 0}};
long double E = (0.5)*m_val[0]*(initial_val[0][2]*initial_val[0][2] + initial_val[0][3]*initial_val[0][3])
            - G*m_val[0]*m_val[1]/Norm(initial_val[0][0]-initial_val[1][0], initial_val[0][1]-initial_val[1][1])
              + (0.5)*m_val[1]*(initial_val[1][2]*initial_val[1][2] + initial_val[1][3]*initial_val[1][3])
            - G*m_val[1]*m_val[2]/Norm(initial_val[1][0]-initial_val[2][0], initial_val[1][1]-initial_val[2][1])
              + (0.5)*m_val[2]*(initial_val[2][2]*initial_val[2][2] + initial_val[2][3]*initial_val[2][3])
            - G*m_val[2]*m_val[0]/Norm(initial_val[2][0]-initial_val[0][0], initial_val[2][1]-initial_val[0][1]);


long double T( long double _xv_array[][4]){
    long double T = E;
    for (int i = 0 ; i < 3 ; i++){
        int j = (i+1)%3;
        T += G*m_val[i]*m_val[j]/Norm(_xv_array[i][0]-_xv_array[j][0], 
                                      _xv_array[i][1]-_xv_array[j][1]);
// cout << _xv_array[i][0] - _xv_array[j][0] << " "<<
//                                       _xv_array[i][1] - _xv_array[j][1] <<endl;
// cout << "T" << Norm(_xv_array[i][0]-_xv_array[j][0], 
//                                       _xv_array[i][1]-_xv_array[j][1]) << endl;
    } 
    return T;
}
long double T_diff(long double _xv_array[][4]){
    long double T_diff = 0;
    for (int i = 0; i < 3; i++){
        int j = (i+1)%3;
        T_diff += -G*m_val[i]*m_val[j]
        *((_xv_array[i][0]-_xv_array[j][0])*(_xv_array[i][2]-_xv_array[j][2])
         +(_xv_array[i][2]-_xv_array[j][2])*(_xv_array[i][3]-_xv_array[j][3]))
        /(powl((_xv_array[i][0]-_xv_array[j][0]), 2)+powl((_xv_array[i][1]-_xv_array[j][1]), 2));
// cout << "T_diff"<< (powl((_xv_array[i][0]-_xv_array[j][0]), 2)+powl((_xv_array[i][1]-_xv_array[j][1]), 2))<<endl;
    }
    return T_diff;
}
long double diff_T_xv(int i, int j, long double _xv_array[][4]){
    long double result = 0;
    for ( int n = 0 ; n < 2 ; n++ ){
        int k = (i+n)%3; int l = (k+1)%3;
        result += powl(-1, n+2)*G*m_val[k]*m_val[l]*2*(_xv_array[k][j]-_xv_array[l][j])
        /( powl( powl((_xv_array[k][0]-_xv_array[l][0]), 2)+powl((_xv_array[l][k]-_xv_array[l][1]), 2), 1.5));
// cout <<"diff_T_xv" <<( powl( powl((_xv_array[k][0]-_xv_array[l][0]), 2)+powl((_xv_array[l][k]-_xv_array[l][1]), 2), 1.5)) <<endl;
    }
    return result;
}
void F( long double _xv_array[][4], long double F[][4]){
    long double dT;
    for (int i = 0; i < 3 ; i++){for( int j = 0 ; j < 2 ; j ++ ){
        F[i][j] = T_diff(_xv_array)/T(_xv_array)*_xv_array[i][j+2];
        for(int k = 0 ; k < 3; k++){for(int l = 0; l < 2; l++){
            dT = diff_T_xv(i, j, _xv_array);
            F[i][j] += - 1/T(_xv_array)*dT*_xv_array[k][l+2]
                        + 1/(2*T(_xv_array))*dT*(_xv_array[k][l+2]*_xv_array[k][l+2]);
        }}}}
}
void geod( long double _X[][4], long double result_ptr[][4]){
    long double F1[3][4]; long double F2[3][4]; 
    for (int i = 0 ; i < 3; i++){
        result_ptr[i][0] = _X[i][2]; result_ptr[i][1] = _X[i][3];
        F(_X, F1); F(_X, F2);
        result_ptr[i][2] = F1[i][0]; result_ptr[i][3] = F2[i][1];
    }
}
void const_multiply_array( long double coeff, long double ptr[][4], long double result_ptr[][4]){
    for (int i = 0 ; i < 3 ; i ++){for (int j = 0 ; j < 4 ; j++ ){
        result_ptr[i][j] = coeff * ptr[i][j];
    }}
}
void sum_two_ptr(long double ptr1[][4], long double ptr2[][4]){
    for (int i = 0 ; i < 3 ; i ++){for (int j = 0 ; j < 4 ; j++ ){
        ptr2[i][j] += ptr1[i][j];
    }}
}


int main(){
    string filename = "123";
    ofstream writefile(filename + ".txt");
    long double X[3][4];
    for (int i = 0 ; i < 3 ; i ++){
        for (int j = 0; j < 4 ; j++){
            X[i][j] = initial_val[i][j];
        }
    }
    writefile << write_contents(X, 0.0) << endl;

    long double k1[3][4];long double k2[3][4];long double k3[3][4];long double k4[3][4];
    long double extra_k[3][4];
    for (long double t = 0.0; t < 0.030; t += dt){
        geod(X, k1);
        const_multiply_array(dt/2, k1, extra_k);sum_two_ptr(X, extra_k);geod(extra_k, k2);
        const_multiply_array(dt/2, k2, extra_k);sum_two_ptr(X, extra_k);geod(extra_k, k3);        
        const_multiply_array(dt  , k3, extra_k);sum_two_ptr(X, extra_k);geod(extra_k, k4);

        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 4; j++){
                X[i][j] = X[i][j] + (dt / 6) * (k1[i][j] + 2 * k2[i][j] + 3 * k3[i][j] + 4 * k4[i][j]);
            }}
        writefile << write_contents(X, t) << endl;
    }
    writefile.close(); cout << "program end!!!"<<endl;
    return 0;
}