// CHMlab2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
const double PI = 3.141592653589793;
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

double func(double x) {
    return sin(x) + x * cos(x);
}

double Rect_Sum(int n) {
    double a = 0, b = PI/2;
    double h = (b - a) / n, ans;
    vector <double> x, y;
    double delta = 0.0001;
    double accurateI = 1.570796327;
    for (int i = 0; i <= n; i++) {
        double tempx = a + i*h;
        x.push_back(tempx);
        y.push_back(func(tempx));
    }
    double rectsum=0;
    for (int i = 1; i < n; i++) {
        rectsum += h * func((x[i - 1] + x[i]) / 2);
    }
    return rectsum;
}

double Rect_Raz(int n) {
    double a = 0, b = PI / 2;
    double h = (b - a) / n, ans;
    vector <double> x, y;
    double delta = 0.0001;
    double accurateI = 1.570796327;
    for (int i = 0; i <= n; i++) {
        double tempx = a + i * h;
        x.push_back(tempx);
        y.push_back(func(tempx));
    }
    double rectsum = 0;
    for (int i = 1; i < n; i++) {
        rectsum += h * func((x[i - 1] + x[i]) / 2);
    }
    ans = abs(accurateI - rectsum);
    return ans;
}

double Trap_Sum(int n) {
    double a = 0, b = PI / 2;
    double h = (b - a) / n, ans;
    vector <double> x, y;
    double delta = 0.0001;
    double accurateI = 1.570796327;
    for (int i = 0; i <= n; i++) {
        double tempx = a + i * h;
        x.push_back(tempx);
        y.push_back(func(tempx));
    }

    double trapsum = h * (func(a) + func(b)) / 2;
    for (int i = 1; i < n; i++) {
        trapsum += h * func(x[i]);
    }

    return trapsum;
}

double Trap_Raz(int n) {
    double a = 0, b = PI / 2;
    double h = (b - a) / n, ans;
    vector <double> x, y;
    double delta = 0.0001;
    double accurateI = 1.570796327;
    for (int i = 0; i <= n; i++) {
        double tempx = a + i * h;
        x.push_back(tempx);
        y.push_back(func(tempx));
    }

    double trapsum = h * (func(a) + func(b)) / 2;
    for (int i = 1; i < n; i++) {
        trapsum += h * func(x[i]);
    }

    ans = abs(accurateI - trapsum);
    return ans;
}

void Ex1_Point2()
{
    ofstream graphout;
    graphout.open("Ex1_Point2_Graphic.txt");

    for (int i = 1; i < 10000; i++)
    {
        graphout <<i << "," << Rect_Raz(i) << "," << Trap_Raz(i) << "\n";
    }
}
void Ex1_Point3()
{
    double needed_delta = 0.0001;
    int flag = 0;
    int temp_n_Rec = 1;
    int temp_n_Trap = 1;
    double temp_delta;

    while (flag == 0)
    {
        temp_delta = Rect_Raz(temp_n_Rec);
        if (temp_delta < needed_delta)
        {
            flag++;
        }
        else
        {
            temp_n_Rec ++;
        }
    }

    flag = 0;

    while (flag == 0)
    {
        temp_delta = Trap_Raz(temp_n_Trap);
        if (temp_delta < needed_delta)
        {
            flag++;
        }
        else
        {
            temp_n_Trap++;
        }

    }
    cout << "n for rect = " << temp_n_Rec << " and for trap n is = " << temp_n_Trap << "\n";
}

double Sim_sum(int n) {
    double a = 0, b = PI / 2;
    double h = (b - a) / double(n), ans;
    vector <double> x, y;
    double accurateI = 1.570796327;
    for (int i = 0; i < n; i++) {
        double tempx = a + i * h;
        x.push_back(tempx);
        y.push_back(func(tempx));
    }

    double tempx = b;
    x.push_back(tempx);
    y.push_back(func(tempx));
  
 
    double simsum = 0;

    if (n % 2 == 0)
    {
        /*for (int i = 1; i < n; i++)
        {
            simsum += h  * (y[i - 1] + 4 * y[i] + y[i + 1])  / 3;
        }*/
        simsum += y[0];
        for (int i = 0; i < n; i++)
        {
            if (i % 2 == 0)
            {
                simsum += 2 * y[i];
            }
            else
            {
                simsum += 4 * y[i];
            }
        }
        simsum += y[n];
        simsum = h * simsum / (double)3;
    }
    else
    {
        cout << "n!=2m error";
    }
    
    ans = abs(accurateI - simsum);
    //cout << ans << endl;
    return ans;
}

void Ex2_Point2()
{
    ofstream graphout;
    graphout.open("Ex2_Point2_Graphic.txt");

    for (int i = 2; i < 1000; i+=2)
    {
        graphout << i << "," << Sim_sum(i) <<  "\n";
    }
}
void Ex2_Point3()
{
    double needed_delta = 0.0001;
    int flag = 0;
    int temp_n_Sim = 2;
    double temp_delta;

    while (flag == 0)
    {
        temp_delta = Sim_sum(temp_n_Sim);
        if (temp_delta < needed_delta)
        {
            flag++;
        }
        else
        {
            temp_n_Sim+=2;
        }
    }

    cout << "n for Simpson = " << temp_n_Sim << "\n";
}

double Ex3_Point1()
{      
    int n = 1000;
    double a = 0, b = PI / 2;
    double h = (b - a) / double(n);
    double n1 = 0, n2 = 0, n3 = 0;
    double q = 0.5;

    double h1 = h;
    for (double tempx = a; tempx <= b; tempx += h1) {
        n1++;
    }
    double h2 = q * h;
    for (double tempx = a; tempx <= b; tempx += h2) {
        n2++;
    }
    double h3 = q * q * h;
    for (double tempx = a; tempx <= b; tempx += h3) {
        n3++;
    }

    double I1 = Rect_Sum(n1), I2 = Rect_Sum(n2), I3 = Rect_Sum(n3);
    //cout << I1 << " " << I2 << " " << I3 << " \n";
    double B = (I1 - I2) / (I2 - I3);
    double p = - log(B) / log(q);
    return p;
}
void Ex3_Point2()
{
    int n = 1000;
    double p = Ex3_Point1();

    cout << " p = " << " " <<  p << "\n";
    // метод рунге 
    double I2n = Rect_Sum(2 * n), In = Rect_Sum(n);
    double I = I2n + (I2n - In) / (pow(2, p) - 1);
    double accurateI = 1.570796327;
    cout.precision(7);
    cout << " my integral = " << I << " actual integral is  " << accurateI << "\n";
}

void Ex3_Point3()
{
    double p = Ex3_Point1();
   // double q = -6 * log(10) / log()

}

double h(int n, double a, double b) {
    return(b - a) / pow(2.0, n);
}
double R(int n, int m, double a, double b) {
    if (m == 0) {
        if (n == 0) {
            return h(1, a, b) * (func(a) + func(b));
        }
        else
        {
            double sum = 0.5 * R(n - 1, 0, a, b);
            for (int k = 1; k <= pow(2, n - 1); k++) {
                sum += h(n, a, b) * func(a + (2 * k - 1) * h(n, a, b));
            }
            return sum;
        }
    }
    else
    {
        return (pow(4.0, m) * R(n, m - 1, a, b) - R(n - 1, m - 1, a, b)) / (pow(4.0, m) - 1);
    }
}

void Ex3_Point3Richardson() {
    double r = 0.0, a = 0.0, b = PI / 2;
    int n=0;
    double accurateI = 1.570796327;
    while (abs(r - accurateI)>0.000001) {
        r = R(n, n,a,b);
        n++;
    }
    cout << "N:" << n - 1 << " " << fixed<<setprecision(10) << r <<" " << abs(r - accurateI);
}
void Table_Data(int n, vector<double>& x, vector<double>& W) {
    if (n == 1)
    {
        x.push_back(0);
        W.push_back(2);
    }
    else if (n == 2)
    {
        x.push_back(-0.5773503);
        x.push_back(0.5773503);
        W.push_back(1);
        W.push_back(1);
    }
    else if (n == 3)
    {
        x.push_back(-0.7745967);
        x.push_back(0);
        x.push_back(0.7745967);
        W.push_back(0.5555556);
        W.push_back(0.8888889);
        W.push_back(0.5555556);
    }
    else if (n == 4)
    {
        x.push_back(-0.8611363);
        x.push_back(-0.3399810);
        x.push_back(0.3399810);
        x.push_back(0.8611363);
        W.push_back(0.3478548);
        W.push_back(0.6521451);
        W.push_back(0.6521451);
        W.push_back(0.3478548);
    }
    else if (n == 5)
    {
        x.push_back(-0.9061798);
        x.push_back(-0.5384693);
        x.push_back(0);
        x.push_back(0.5384693);
        x.push_back(0.9061798);
        W.push_back(0.2369269);
        W.push_back(0.4786287);
        W.push_back(0.5688888);
        W.push_back(0.4786287);
        W.push_back(0.2369269);
    }
    else if (n == 6)
    {
        x.push_back(-0.9324700);
        x.push_back(-0.6612094);
        x.push_back(-0.2386142);
        x.push_back(0.2386142);
        x.push_back(0.6612094);
        x.push_back(0.9324700);
        W.push_back(0.1713245);
        W.push_back(0.3607616);
        W.push_back(0.4679140);
        W.push_back(0.4679140);
        W.push_back(0.3607616);
        W.push_back(0.1713245);
    }
    else if (n == 7)
    {
        x.push_back(-0.9491079);
        x.push_back(-0.7415312);
        x.push_back(-0.4058452);
        x.push_back(0);
        x.push_back(0.4058452);
        x.push_back(0.7415312);
        x.push_back(0.9491079);
        W.push_back(0.1294850);
        W.push_back(0.2797054);
        W.push_back(0.3818301);
        W.push_back(0.4179592);
        W.push_back(0.3818301);
        W.push_back(0.2797054);
        W.push_back(0.1294850);
    }
}
void Ex4() {
    //t=4x/pi-1
    int minstep = 0;
    // n - количество узлов сетки, x- вектор узлов, W- вектор весов 
    double a = 0.0,b=PI/2,delta=0.000001,accurateI= 1.570796327;
    for (int n = 2; n < 8; n++)
    {
        float Integral = 0;
        vector<double> x, W;
        Table_Data(n, x, W);
        for (int i = 0; i < n; i++)
        {
            Integral = Integral + func(PI/4 * x[i] + PI/4) * W[i];
        }
        Integral = Integral * (b - a) / 2;
        if (abs(accurateI - Integral) < delta)
        {
            minstep = n;
            cout << "Minimum n for Gauss: " << n << endl;
            cout << abs(accurateI - Integral) << endl;
            //cout << "Minimum n for Gauss: " << n << endl;
            break;
        }
    }
}

void Ex4_Graphic() {
    //t=4x/pi-1
    // n - количество узлов сетки, x- вектор узлов, W- вектор весов 
    ofstream out;
    out.open("Ex4_Graphic.txt");
    double a = 0.0, b = PI / 2, delta = 0.000001, accurateI = 1.570796327;
    for (int n = 2; n < 8; n++)
    {
        float Integral = 0;
        vector<double> x, W;
        Table_Data(n, x, W);
        for (int i = 0; i < n; i++)
        {
            Integral = Integral + func(PI / 4 * x[i] + PI / 4) * W[i];
        }
        Integral = Integral * (b - a) / 2;

        out << n << " " << abs(accurateI - Integral) << endl;
    }

}
double Int_Trap( double x1, double x2)
{
    return (func(x1) + func(x2)) / 2 * (x2 - x1);
}
void Ex5()
{
    double a = 0.0, b = PI / 2, delta = 0.0001, accurateI = 1.570796327;
    double S = 0;
    int n = 1000;
    vector<double> x;
    vector<double> h;
    x.push_back(a);
    h.push_back( 0.5 );
    x.push_back( a + h[0] );
    int i = 1;
    while( 1 )
    {
        
        double Int_Trap = h[i - 1] / 2 * ( func( x[i - 1] ) + func( x[i] )  );
        double Int_h_i = h[i - 1] / 4 *  ( func( x[i - 1] ) + 2 * func( (x[i - 1] + x[i]) / 2 ) + func( x[i]) );
        //cout <<  "int" << " " << Int_h_i << endl;
        double Ei = (Int_h_i - Int_Trap) / 3;
        //cout << "ei " <<  Ei << endl;
        //cout << h[i - 1] * delta / (b - a) << endl;
        //double proverka = h[i - 1] * delta / (b - a);
        if (abs(Ei) > h[i-1] * delta / (b - a))
        {
            h[i - 1] = h[i - 1] / 2;
            x[i] = x[i - 1] + h[i - 1];
        }
        else
        { 

            x[i] = x[i - 1] + h[i - 1];
            S += Int_h_i;
            if (x[i] + h[i - 1] > b)
            {
                h.push_back(b - x[i]);
            }
            else
            {
                h.push_back( h[i-1] );
                //cout << h[i] << endl;
            }
            x.push_back(x[i] + h[i]);
            i++;
            //cout << x[i] << " " << S << endl;
            cout << S << endl;
            if (abs(S - accurateI) < delta) break;
        }
    }
    cout << x.size();
    ofstream out;
    out.open("Ex5_Graphic.txt");
    for (int i = 0; i < x.size(); i++)
    {
        out << x[i] << " " << func(x[i]) << endl;
    }
    out.close();

}

void Ex6()
{
    ofstream out;
    out.open("Ex6_Graphic.txt");
    double a = 0.0, b = PI / 2, accurateI = 1.570796327;
    double xk = RAND_MAX / (b - a);
    srand(time(NULL));
    for (int n = 1; n < 2000; n++)
    {
        double S = 0.0;
        for (int k = 0; k < 100; k++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                double x = rand() / xk + a;
                sum += func(x);
            }
            sum = sum * (b - a) / n;
            S += sum;
        }
        S /= 100;
       
        //cout << S<<"\t"<< "Raz:"<<abs(accurateI-S)<< endl;
        out << n << " " << abs(accurateI - S) << endl;
    }
}

int main()
{
    setlocale(LC_CTYPE, "Russian");
    Ex1_Point2();
    //Ex1_Point3();

    //Ex2_Point2();
    //Ex2_Point3();
    //Ex3_Point2();
    
    //Ex3_Point3Richardson();
    //Ex4();
    //Ex4_Graphic();
    //Ex5();
    //Ex6();
    return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
