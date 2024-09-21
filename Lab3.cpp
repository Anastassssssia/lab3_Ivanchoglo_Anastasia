#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>  
#include <chrono>
#include <stdlib.h>
#include <Windows.h>
#include <math.h>
using namespace std;

class Cartesian2D;
class Cartesian3D;
class Polar;
class Spherical;

//===============================================================================
class Cartesian2D {
public:

    double x, y;

    Cartesian2D() {
        x = 0;
        y = 0;
    }

    Cartesian2D(double _x, double _y) {
        x = _x;
        y = _y;
    }

    double distance(Cartesian2D point) {
        return sqrt(pow(point.x - x, 2) + pow(point.y - y, 2));
    }

    Polar ToPolar();
};
//===============================================================================
class Cartesian3D {
public:

    double x, y, z;

    Cartesian3D() {
        x = 0;
        y = 0;
        z = 0;
    }

    Cartesian3D(double _x, double _y, double _z) {
        x = _x;
        y = _y;
        z = _z;
    }

    double distance(Cartesian3D point) const {
        return sqrt(pow(point.x - x, 2) + pow(point.y - y, 2) + pow(point.z - z, 2));
    }

    Spherical ToSpherical() const;
};
//===============================================================================
class Polar {
public:

    double r, theta;

    Polar() {
        r = 0;
        theta = 0;
    }

    Polar(double _r, double _theta) {
        r = _r;
        theta = _theta;
    }

    double distance(Polar circle) const {
        return sqrt(r * r + circle.r * circle.r - 2 * r * circle.r * cos(circle.theta - theta));
    }

    Cartesian2D ToCartesian2D() const;
};
//===============================================================================
class Spherical {
public:

    double  r, theta, phi;

    Spherical() {
        r = 0;
        theta = 0;
        phi = 0;
    }

    Spherical(double _r, double _theta, double _phi) {
        r = _r;
        theta = _theta;
        phi = _phi;
    }

    double distanceV(Spherical spherical) const {
        return sqrt(r * r + spherical.r * spherical.r - 2 * r * spherical.r * (sin(theta) * sin(spherical.theta) * cos(phi - spherical.phi) + cos(theta) * cos(spherical.theta)));
    }

    double distanceC(Spherical spherical) const {
        return r * acos(sin(phi) * sin(spherical.phi) + cos(phi) * cos(spherical.phi) * cos(theta - spherical.theta));
    }

    Cartesian3D ToCartesian3D() const;
};
//===============================================================================

double radToDeg(double radians) {
    return radians * 180.0 / M_PI;
}
//===============================================================================
double degToRad(double degrees) {
    return degrees * M_PI / 180.0;
}
//===============================================================================


Polar Cartesian2D::ToPolar() {
    double r = sqrt(x * x + y * y);
    double theta = radToDeg(atan2(y, x));
    return Polar{ r , theta };
}
//===============================================================================
Spherical Cartesian3D::ToSpherical() const {
    double r = sqrt(x * x + y * y + z * z);
    double theta = radToDeg(atan2(y, x));
    double phi = radToDeg(acos(z / r));
    return Spherical{ r, theta, phi };
}
//===============================================================================
Cartesian2D Polar::ToCartesian2D() const {
    double x = r * cos(degToRad(theta));
    double y = r * sin(degToRad(theta));
    return Cartesian2D{ x,y };
}
//===============================================================================
Cartesian3D Spherical::ToCartesian3D() const {
    double x = r * sin(degToRad(phi)) * cos(degToRad(theta));
    double y = r * sin(degToRad(phi)) * sin(degToRad(theta));
    double z = r * cos(degToRad(phi));
    return Cartesian3D{ x,y,z };
}
//===============================================================================


int main() {
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);

    Polar A{ 5, 30 };
    Polar B{ 3, 120 };

    cout << "Задані точки в полярній системі:" << endl;
    cout << "A (r = " << A.r << ",  0 = " << A.theta << "°)\n";
    cout << "B (r = " << B.r << ", 0 = " << B.theta << "°)\n\n";

    //---------------------------------------------------------------------------------

    Cartesian2D A_2D = A.ToCartesian2D();
    Cartesian2D B_2D = B.ToCartesian2D();

    cout << "Результати перетворення в декартову систему:" << endl;
    cout << "A (x = " << A_2D.x << ", y = " << A_2D.y << ")\n";
    cout << "B (x = " << B_2D.x << ", y = " << B_2D.y << ")\n\n";

    //---------------------------------------------------------------------------------

    Polar A_back = A_2D.ToPolar();
    Polar B_back = B_2D.ToPolar();

    cout << "Координати точок після зворотного перетворення в полярну систему:" << endl;
    cout << "A (r = " << A_back.r << ", 0 = " << A_back.theta << "°)\n";
    cout << "B (r = " << B_back.r << ", 0 = " << B_back.theta << "°)\n\n";

    //---------------------------------------------------------------------------------

    cout << "Перевірка коректності: " << endl;
    cout << "A: Радіус: " << (abs(A.r - A_back.r) <= 0.001 ? "співпадає" : "не співпадає")
        << ", Кут: " << (abs(A.theta - A_back.theta) <= 0.001 ? "співпадає" : "не співпадає") << "\n";
    cout << "B: Радіус: " << (abs(B.r - B_back.r) <= 0.001 ? "співпадає" : "не співпадає")
        << ", Кут: " << (abs(B.theta - B_back.theta) <= 0.001 ? "співпадає" : "не співпадає") << "\n\n";

    //---------------------------------------------------------------------------------

    Spherical C{ 5 ,45 ,30 };
    Spherical D{ 3 ,60 ,60 };

    cout << "Задані точки в сферичній системі:" << endl;
    cout << "C (r = " << C.r << ",  0 = " << C.theta << "°" << ",  ф = " << C.phi << "°)\n";
    cout << "D (r = " << D.r << ", 0 = " << D.theta << "°" << ",  ф = " << D.phi << "°)\n\n";

    //---------------------------------------------------------------------------------

    Cartesian3D C_3D = C.ToCartesian3D();
    Cartesian3D D_3D = D.ToCartesian3D();

    cout << "Результати перетворення в декартову систему:" << endl;
    cout << "C (x = " << C_3D.x << ", y = " << C_3D.y << ", z = " << C_3D.z << ")\n";
    cout << "D (x = " << D_3D.x << ", y = " << D_3D.y << ", z = " << D_3D.z << ")\n\n";

    //---------------------------------------------------------------------------------

    Spherical C_back = C_3D.ToSpherical();
    Spherical D_back = D_3D.ToSpherical();

    cout << "Координати точок після зворотного перетворення в сферичну систему:" << endl;
    cout << "C (r = " << C_back.r << ", 0 = " << C_back.theta << "°" << ",  ф = " << C_back.phi << "°)\n";
    cout << "D (r = " << D_back.r << ", 0 = " << D_back.theta << "°" << ",  ф = " << D_back.phi << "°)\n\n";

    //---------------------------------------------------------------------------------

    cout << "Перевірка коректності: " << endl;
    cout << "C: Радіус: " << (abs(C.r - C_back.r) <= 0.001 ? "співпадає" : "не співпадає")
        << ", Кут 0: " << (abs(C.theta - C_back.theta) <= 0.001 ? "співпадає" : "не співпадає")
        << ", Кут ф: " << (abs(C.phi - C_back.phi) <= 0.001 ? "співпадає" : "не співпадає") << "\n";
    cout << "D: Радіус: " << (abs(D.r - D_back.r) <= 0.001 ? "співпадає" : "не співпадає")
        << ", Кут 0: " << (abs(D.theta - D_back.theta) <= 0.001 ? "співпадає" : "не співпадає")
        << ", Кут ф: " << (abs(D.phi - D_back.phi) <= 0.001 ? "співпадає" : "не співпадає") << "\n\n";

    //---------------------------------------------------------------------------------

    double dist_cartesian_2D = A_2D.distance(B_2D);
    cout << "Відстань між A і B в декартовій системі (2D): " << dist_cartesian_2D << "\n";

    double dist_polar = A.distance(B);
    cout << "Відстань між A і B в полярній системі: " << dist_polar << "\n";

    double dist_cartesian_3D = C_3D.distance(D_3D);
    cout << "Відстань між C і D в декартовій системі (3D): " << dist_cartesian_3D << "\n";

    double dist_spherical_v = C.distanceV(D);
    cout << "Відстань між C і D в сферичній системі через об'єм сфери: " << dist_spherical_v << "\n";

    double dist_spherical_c = C.distanceC(D);
    cout << "Велика колова відстань між C і D на сфері: " << dist_spherical_c << "\n\n";

    //---------------------------------------------------------------------------------

    const int size = 10000;
    double x, y, z, distance;

    Cartesian2D Arr_Cartesian2D[size];
    Cartesian3D Arr_Cartesian3D[size];
    Polar Arr_Polar[size];
    Spherical Arr_Spherical[size];


    //===============================================================================
    int start = 1, end = 999;

    for (int i = 0; i < size; i++) {
        x = rand() % (end - start + 1) + start;
        y = rand() % (end - start + 1) + start;
        z = rand() % (end - start + 1) + start;
        Arr_Cartesian2D[i] = Cartesian2D{ x, y };
        Arr_Cartesian3D[i] = Cartesian3D{ x, y, z };
        Arr_Polar[i] = Polar{ x, y };
        Arr_Spherical[i] = Spherical{ x, y, z };
    };

    //===============================================================================
    auto tStart = chrono::high_resolution_clock::now();
    for (int i = 0; i < size - 1; i++) {
        distance = Arr_Cartesian2D[i].distance(Arr_Cartesian2D[i + 1]);
    };
    auto tEnd = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = tEnd - tStart;
    cout << "Час для перетворення з полярної в декартову (2D): " << duration.count() << "\n";

    //===============================================================================
    tStart = chrono::high_resolution_clock::now();
    for (int i = 0; i < size - 1; i++) {
        distance = Arr_Cartesian3D[i].distance(Arr_Cartesian3D[i + 1]);
    };
    tEnd = chrono::high_resolution_clock::now();
    duration = tEnd - tStart;
    cout << "Час для перетворення з сферичної в декартову (3D): " << duration.count() << "\n";

    //===============================================================================
    tStart = chrono::high_resolution_clock::now();
    for (int i = 0; i < size - 1; i++) {
        distance = Arr_Polar[i].distance(Arr_Polar[i + 1]);
    };
    tEnd = chrono::high_resolution_clock::now();
    duration = tEnd - tStart;
    cout << "Час для розрахунку відстаней у полярній системі: " << duration.count() << "\n";

    //===============================================================================
    tStart = chrono::high_resolution_clock::now();
    for (int i = 0; i < size - 1; i++) {
        distance = Arr_Spherical[i].distanceV(Arr_Spherical[i + 1]);
    };
    tEnd = chrono::high_resolution_clock::now();
    duration = tEnd - tStart;
    cout << "Час для розрахунку відстаней у сферичній системі (об'єм): " << duration.count() << "\n";

    //===============================================================================
    tStart = chrono::high_resolution_clock::now();
    for (int i = 0; i < size - 1; i++) {
        distance = Arr_Spherical[i].distanceC(Arr_Spherical[i + 1]);
    };
    tEnd = chrono::high_resolution_clock::now();
    duration = tEnd - tStart;
    cout << "Час для розрахунку відстаней у сферичній системі (велика колова відстань): " << duration.count() << "\n";


    return 0;
}



