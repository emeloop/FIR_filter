#include<iostream>
#include<math.h>
#include<vector>
#include<tuple>
#include<matplotlib-cpp-starter/matplotlibcpp.h>
#include <algorithm> 
#include <iterator>  
#include <algorithm>
#include<eigen-eigen-323c052e1731/Eigen/Dense>
#include<eigen-eigen-323c052e1731/Eigen/Core>
#include<eigen-eigen-323c052e1731/Eigen/LU>

using std::cout;
using std::endl;
using std::cin;

float lowpass_filter(float _x, float _fc)
{
    if(_x < _fc)return 1.0;
    else return 0.0;
}

float tri_polynomial(float _x, std::vector<float> _list_a)
{
    int a=0,k=0,sum=0;
    for(int i=0; i<_list_a.size();i++)
    {
        sum += _list_a[i]*cos(_x*k);
        k++;
    }
    return sum;
}
float d_tri_polynomial(float _x, std::vector<float> _list_a)
{
    int a=0,k=0,sum=0;
    for(int i=0; i<_list_a.size();i++)
    {
        sum += -_list_a[i]*k*sin(_x*k);
        k++;
    }
    return sum;

}
float dd_tri_polynomial(float _x, std::vector<float> _list_a)
{
    int a=0,k=0,sum=0;
    for(int i=0; i<_list_a.size();i++)
    {
        sum += -_list_a[i]*k*k*cos(_x*k);
        k++;
    }
    return sum;
}
float func(float x, std::vector<float> _list_a)
{
    return d_tri_polynomial(x, _list_a);
}
float dfunc(float x, std::vector<float> _list_a)
{
    return dd_tri_polynomial(x, _list_a);
}
std::vector<float> linspace(float a, float b, int n) 
{
    std::vector<float> array;
    float step = (b-a) / (n-1);

    while(a <= b) {
        array.push_back(a);
        a += step;           // could recode to better handle rounding errors
    }
    return array;
}
float newton(float _x0, std::vector<float> _list_a)
{
    float x,x_next;
    x = _x0;
    for(int i = 0; i < 100; i++)
    {
        x_next = x - func(x,_list_a)/dfunc(x, _list_a);
        if(abs(x_next - x) < 1.0*pow(10.0, -8)) break;
        else
        {
            cout << "[WARNING] Newton method iterations reached the limit (x=" << _x0 << ")" << endl;
        }
        x = x_next;
    }
    return x;
}

std::vector<float> search_extreme_points(std::vector<float> _list_a, int _div)
{
    float func,dfunc,x;
    //    func = d_tri_polynomial(x, _list_a);
    //    dfunc = dd_tri_polynomial(x, _list_a);
    //    newton(x0,_list_a);
    std::vector<float> check_points;
    std::copy(linspace(0.0, M_PI, _div).begin(), linspace(0.0, M_PI, _div).end(), std::back_inserter(check_points) );
    std::vector<float> sign_reverse_section;
    std::tuple<float, float> p;
    for(int i = 0; i < check_points.size();i++)
    {
        p = std::make_tuple(check_points[i], check_points[i+1]);
        //if(func(p[]))
    }
    std::vector<float> ans;
    for(int i=0;i<sign_reverse_section.size();i++)
    {
        int x;
        x = sign_reverse_section[i];
        ans.push_back(newton(x, _list_a));
    }
    return ans;
}
std::vector<float> initialize_extreme_points(int _n, float _fc, float _fs)
{
    int num_point = _n + 2;
    int num_passed_point = (int)(num_point * _fc/M_PI);
    int num_stopped_point = num_point-num_passed_point;
    cout << "num_passed_point:" << num_passed_point <<endl;
    cout << "num_stopped_point:" << num_stopped_point << endl;
    std::vector<float> num;
    std::vector<float> linscape1 = linspace(0.0, _fc-0.5*_fs, num_passed_point);
    std::vector<float> linscape2 = linspace(_fs+0.5*_fs, M_PI, num_stopped_point);
    cout << "linscape1:" << linscape1.size() << endl;
    cout << "linscape2:" << linscape2.size() << endl;
    for(int i=0; i<linscape1.size();i++)
    {
        num.push_back(linscape1[i]);
    }
    for(int i=0; i<linscape2.size();i++)
    {
        num.push_back(linscape2[i]);
    }
    return num;
}
std::tuple<std::vector<float>, float> update_tri_polynomial_coefficients(std::vector<float> _list_x, std::vector<float> _list_a, float _fc)
{
    cout << "test1" << endl;
    int x,k=0;
    int x_size = _list_x.size();
    cout << "x_size:" << _list_x.size() << endl;
    Eigen::MatrixXd matrix_A(x_size,x_size+1);
    /*
    for(x=0; x<x_size-1; x++)
    {
        matrix_A << cos(x*k);
    }
    cout << "test" << endl;
    for(int j=0; j<x_size; j++)
    {
        matrix_J << pow(-1,j);
    }
    cout << "test" << endl;
    matrix_A = matrix_A + matrix_J;
    cout << "test" << endl;
    */
    int j=0;
    int count = 1;
    for(x=0;x<x_size;x++)
    {
        for(k=0;k<x_size-1;k++)
        {
            //cout << "list_x["<< j << "]:" << _list_x[j] << endl;
            //cout << matrix_A << endl;
            matrix_A(x,k) = cos(_list_x[j]*k);
            j++;
            cout << count++ << endl;
        }
        cout << "test" << endl;
        matrix_A(x,x_size) = pow(-1.0,(double)x);
        j=0;
    }
    cout << "test1" <<endl;
    Eigen::VectorXd vector_b(x_size);
    for(int i=0; i<x_size; i++)
    {
        cout << "test2" << endl;
        vector_b(i,0) = lowpass_filter(_list_x[i], _fc);
    }
    cout << "test" << endl;
    cout << "matrix_A:" << matrix_A << endl;
    Eigen::VectorXd vector_x(x_size);
    vector_x = matrix_A.fullPivLu().solve(vector_b);
    cout << "vector_x:" << vector_x << endl;
    cout << "test" << endl;
    for(int i=0; i<vector_x.size();i++)
    {
        _list_a[i] = vector_x(i);
    }
    cout << "test" << endl;
    float d = vector_x(vector_x.size()-1);
    cout << "test" << endl;
}
std::vector<float> update_maximum_error_points(std::vector<float> _list_a, float _fc, float _fs)
{
    int n = _list_a.size()-1;
    std::vector<float> extreme_points;
    std::copy(search_extreme_points(_list_a,(n+2)*10).begin(), search_extreme_points(_list_a,(n+2)*10).end(), std::back_inserter(extreme_points));
    extreme_points.push_back(_fc-_fs*0.5);
    extreme_points.push_back(_fc+_fs*0.5);
    std::sort(extreme_points.begin(), extreme_points.end());

    if(extreme_points.size() == n+1)
    {
        extreme_points.push_back(M_PI);
        return extreme_points;
    }
    else if(extreme_points.size() == n+2)
    {
        return extreme_points;
    }    
    else if(extreme_points.size() == n+3)
    {
        assert(!extreme_points.empty());
        extreme_points.erase(extreme_points.begin());
        return extreme_points;
    }
    else
    {
        cout << "[ERROR] number of extreme point. " << endl;
    }
}
float efunc(float _x, std::vector<float> _list_a, float _fc)
{
    return tri_polynomial(_x, _list_a) - lowpass_filter(_x, _fc);
}
float check_convergence(std::vector<float> _list_a, std::vector<float> _list_x, float _fc)
{
    int k=0,x;
    std::vector<float> efunc_array;
    for(int i = 0; i<_list_x.size(); i++)
    {
        x = _list_x[i];
        efunc_array.push_back(efunc(x,_list_a,_fc)*pow((-1),k));
    }
    float x_average,x_sum,x_dispersion;
    for(int i=0;i<efunc_array.size();i++)
    {
        x_sum+=efunc_array[i];
    }
    x_average = x_sum/efunc_array.size();
    x_sum=0;
    for(int i=0;i<efunc_array.size();i++)
    {
        x = _list_x[i];
        x_sum += pow(x - x_average, 2);
    }
    x_dispersion = x_sum/efunc_array.size();
    return x_dispersion < 1.0*pow(10,-12);
}

std::tuple<std::vector<float>, float, std::vector<float>,int> remez(int _order, float _fc, float _fs, int max_iter=100)
{
    int n = (_order - 1) /2;
    std::vector<float> list_x = initialize_extreme_points(n, _fc, _fs);
    cout << "list_x size:" << list_x.size() << endl;
    std::vector<float> list_a,list_h;
    float d;
    int i;
    for(i=1; i<max_iter+1;i++)
    {
        std::tie(list_a, d) = update_tri_polynomial_coefficients(list_x, list_a,_fs);
        int j = list_a.size();
        list_x = update_maximum_error_points(list_a, _fc, _fs);
        if(check_convergence(list_a, list_x, _fc))
        {
            for(int k=0; k<list_a.size();k++)
            {
                list_h.push_back(list_a[j-1]*0.5 + list_a[0] + list_a[k]*0.5);
                j--;
            }
            return forward_as_tuple(list_h, d, list_x, i);
        }
        else
        {
            cout << "[ERROR] Remez algorithm failed." << endl;
        }
    }
}
int main(void)
{
    int filter_order;
    cout << "filter order (odd number): " << endl;
    cin >> filter_order;
    if(filter_order < 3 || filter_order % 2 == 0)
    {
        cout << "inuput add number" << endl;
    }
    float fc;
    cout << "cutoff frequency fc:" << endl;
    cin >> fc;
    if(fc <= 0.0 || M_PI <= fc)
    {
        cout << "[ERROR] Please input range(0vector<tuple>.0,pi)." << endl; 
    }
    float fs;
    cout << "tarnsition width fs:" << endl;
    cin >> fs;
    if(fc-fs*0.5 <= 0.0 || M_PI <= fc-fs*0.5)
    {
        cout << "[ERROR] Transition width is too large." << endl; 
    }

    std::vector<float> list_h;
    float d;
    std::vector<float> list_x;
    int count;
    std::tie(list_h, d, list_x, count) = remez(filter_order, fc, fs);

    for(int i=0;i<list_h.size();i++)
    {
        cout << list_h[i] << endl;
    }
    
}
