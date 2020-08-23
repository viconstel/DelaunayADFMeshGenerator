// Prob1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include "IGESLoader.h"
#include "ParametricGeometryModel.h"
#include "IGESData.h"
#include "IGESDataTypes.h"
#include "DelaunayTriangulation.h"
#include <Eigen/Dense>
#include <vector>

template<typename _T>
void Print(std::map<unsigned, _T> a) {
    for (auto i = a.begin(); i = a.end(); ++i) {
        std::cout << *i.second() << std::endl;
    }
}

double ScalarMulV3(double* v1, double* v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

double* VecMulV3(double* v1, double* v2) {
    double* res = new double[3];
    res[0] = v1[1] * v2[2] - v1[2] * v2[1];
    res[1] = v1[2] * v2[0] - v1[0] * v2[2];
    res[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return res;
}

void NormV3(double* v) {
    double vv = pow(v[0] * v[0] + v[1] * v[1] + v[2] * v[2], .5);
    if (vv < 0. || vv > 0.) {
        v[0] = v[0] / vv;
        v[1] = v[1] / vv;
        v[2] = v[2] / vv;
    }
}

void PrincipalCurvatures(IGESData* data, double coord1, double coord2) {
    const auto& s = data->getNURBSSurfaces();
    auto sur = *(s.begin());
    double res[18];

    //double l2[3] = { 0, 0, -1 };
    //double l3[3] = { 0, 0, 0 };

    //sur.second.getDer1Point(0, 0, res);
    //sur.second.getDer2Point(0, 0, res);
    sur.second.getPointAndDerivs(2, coord1, coord2, res);

    //double l1[3] = { res[9], res[10], res[11] };
    //double l2[3] = { res[3], res[4], res[5] };
    double* l3 = new double[3];
    double* l2 = new double[3];
    double* l1 = new double[3];
    sur.second.getDer1Point(coord1, coord2, l1);
    sur.second.getDer2Point(coord1, coord2, l2);
    l3 = VecMulV3(l1, l2);

    NormV3(l3);
    //NormV3(l2);
    //NormV3(l1);

    std::cout << "coord l1 " << *l1 << " " << *(l1 + 1) << " " << *(l1 + 2) << std::endl;
    std::cout << "coord l2 " << *l2 << " " << *(l2 + 1) << " " << *(l2 + 2) << std::endl;
    std::cout << "coord l3 " << *l3 << " " << *(l3 + 1) << " " << *(l3 + 2) << std::endl;

    double Suu[3] = { res[15], res[16], res[17] };
    double Svv[3] = { res[6], res[7], res[8] };
    double Suv[3] = { res[12], res[13], res[14] };

    //NormV3(Suu);
    //NormV3(Suv);
    //NormV3(Svv);

    double L[4] = { ScalarMulV3(Suu, l3), ScalarMulV3(Suv, l3),ScalarMulV3(Suv, l3),ScalarMulV3(Svv, l3) };
    double G[4] = { ScalarMulV3(l1, l1), ScalarMulV3(l1, l2),ScalarMulV3(l1, l2),ScalarMulV3(l2, l2) };

    double K = (L[0] * L[3] - L[2] * L[1]) / (G[0] * G[3] - G[2] * G[1]);

    double point[3];
    sur.second.getPoint(coord1, coord2, point);
    std::cout << "Point Coord: " << point[0] << " " << point[1] << " " << point[2] << std::endl;
    std::cout << "radius: " << pow(point[0] * point[0] + point[1] * point[1] + point[2] * point[2], 0.5) << std::endl;
    std::cout <<"Gauss Curvature: " << K << std::endl;

    //solving part
    // вадратное уравнение на  ривизны

    double A = G[0] * G[3] - G[2] * G[1];
    double B = -(L[0] * G[3]) - (G[0] * L[3]); //+ (G[1] * L[2] + G[2] * L[1]);       //вторую скобку можно опустить
    double C = L[0] * L[3] - L[2] * L[1];

    double lam_1;
    double lam_2;

    std::cout << "Principal Curatives eq: " << A << "x^2 + " << B << "x + " << C << " = 0" << std::endl;
    std::cout << pow(B, 2) - 4 * A * C;
    if (abs(pow(B, 2) - 4 * A * C) < 0.0000001) {
        std::cout << "Correct equation, rational curvatives" << std::endl;
        lam_1 = (-B + pow(pow(B, 2) - 4 * A * C, 2)) / (2 * A);
        lam_2 = (-B - pow(pow(B, 2) - 4 * A * C, 2)) / (2 * A);
        std::cout << "Roots: lam_1 = " << lam_1 << " lam_2 = " << lam_2 << std::endl;
    } else {
        std::cout << "Incorrect equation" << std::endl;
    }

    Eigen::Matrix2d Lw1, Lw2;
    Eigen::Vector2d b;

    b << 0., 0.;
    Lw1 << L[0] - lam_1 * G[0], L[1] - lam_1 * G[1], L[2] - lam_1 * G[2], L[3] - lam_1 * G[3];
    Eigen::Vector2d eigen_vec1 = Lw1.colPivHouseholderQr().solve(b);
    Lw2 << L[0] - lam_2 * G[0], L[1] - lam_2 * G[1], L[2] - lam_2 * G[2], L[3] - lam_2 * G[3];
    Eigen::Vector2d eigen_vec2 = Lw2.colPivHouseholderQr().solve(b);

    std::cout << "Principal curvature lam_1: " << std::endl << eigen_vec1 << " Principal curvature lam_1: " << std::endl << eigen_vec2 << std::endl;

    delete[] l3;
    std::cout << std::endl;
}

double CurveLength(double eps, NURBSSurface surface, double A[], double B[]) {
    double l1[3];
    double l2[3];
    double r = 0;
    double dt = 0.1;
    double AB[2] = { B[0] - A[0], B[1] - A[1] };
    for (double t = 0; t < 1 - eps; t += dt) {
        surface.getDer1Point(A[0] + t * AB[0], A[1] + t * AB[1], l1);
        surface.getDer2Point(A[0] + t * AB[0], A[1] + t * AB[1], l2);
        double metric_A[2][2] = { {ScalarMulV3(l1,l1), ScalarMulV3(l1, l2)}, {ScalarMulV3(l2, l1), ScalarMulV3(l2, l2) } };
        r += pow((AB[0] * metric_A[0][0] + AB[1] * metric_A[1][0]) * AB[0] + (AB[0] * metric_A[0][1] + AB[1] * metric_A[1][1]) * AB[1], 0.5) * dt;
    }
    return r;
}

void SubDivideSegmentPositive(std::vector<double> &points, double metr, double b[], double A_old[], double B_old[], double unit_vector[], NURBSSurface surface, bool flag, double eps) {
    double A[2] = { A_old[0], A_old[1] };
    double B[2] = { B_old[0], B_old[1] };
    double sigma = 0;
    double r = 0;
    double n = 0;
    double l1[3];
    double l2[3];
    while ((B[0] <= b[0]) && (B[1] <= b[1])) {
        n += 1;
        points.push_back(A[0]);
        points.push_back(A[1]);
        surface.getDer1Point(A[0], A[1], l1);
        surface.getDer2Point(A[0], A[1], l2);
        double metric_A[2][2] = { {ScalarMulV3(l1,l1), ScalarMulV3(l1, l2)}, {ScalarMulV3(l2, l1), ScalarMulV3(l2, l2) } };
        double s = metr / pow((unit_vector[0] * metric_A[0][0] + unit_vector[1] * metric_A[1][0]) * unit_vector[0] + (unit_vector[0] * metric_A[0][1] + unit_vector[1] * metric_A[1][1]) * unit_vector[1], 0.5);
        while (1) {
            B[0] = A[0] + s * unit_vector[0];
            B[1] = A[1] + s * unit_vector[1];
            r = CurveLength(eps, surface, A, B);
            s = 2 * s / (metr + r);
            if ((metr - r) < eps) { break; }
        }
        A[0] = B[0];
        A[1] = B[1];
        if (B[0] > b[0]-eps && B[0] < b[0] + eps && B[1] > b[1] - eps && B[1] < b[1] + eps) { break; }
        if (flag && ((B[0] > b[0]) || (B[1] > b[1]))) {
            sigma = CurveLength(eps, surface, b, B) / n;
            if (sigma > 0) {
                n = n * 2;
                while (n > 0) {
                    points.pop_back();
                    n--;
                }
                SubDivideSegmentPositive(points, metr - sigma, b, A_old, B_old, unit_vector, surface, false, eps); 
            }
        }
    }
}

void SubDivideSegmentNegative(std::vector<double> &points, double metr, double b[], double A_old[], double B_old[], double unit_vector[], NURBSSurface surface, bool flag, double eps) {
    double A[2] = { A_old[0], A_old[1] };
    double B[2] = { B_old[0], B_old[1] };
    double sigma = 0;
    double r = 0;
    double n = 0;
    double l1[3];
    double l2[3];
    while ((B[0] >= b[0]) && (B[1] >= b[1])) {
        n += 1;
        points.push_back(A[0]);
        points.push_back(A[1]);
        surface.getDer1Point(A[0], A[1], l1);
        surface.getDer2Point(A[0], A[1], l2);
        double metric_A[2][2] = { {ScalarMulV3(l1,l1), ScalarMulV3(l1, l2)}, {ScalarMulV3(l2, l1), ScalarMulV3(l2, l2) } };
        double s = metr / pow((unit_vector[0] * metric_A[0][0] + unit_vector[1] * metric_A[1][0]) * unit_vector[0] + (unit_vector[0] * metric_A[0][1] + unit_vector[1] * metric_A[1][1]) * unit_vector[1], 0.5);
        while (1) {
            B[0] = A[0] + s * unit_vector[0];
            B[1] = A[1] + s * unit_vector[1];
            r = CurveLength(eps, surface, A, B);
            s = 2 * s / (metr + r);
            if ((metr - r) < eps) { break; }
        }
        A[0] = B[0];
        A[1] = B[1];
        if (B[0] > b[0] - eps && B[0] < b[0] + eps && B[1] > b[1] - eps && B[1] < b[1] + eps) { break; }
        if (flag && ((B[0] < b[0]) || (B[1] < b[1]))) {
            sigma = CurveLength(eps, surface, B, b) / n;
            if (sigma > 0) { 
                n = n * 2;
                while (n > 0) {
                    points.pop_back();
                    n--;
                }
                SubDivideSegmentNegative(points, metr - sigma, b, A_old, B_old, unit_vector, surface, false, eps); 
            }
        }
    }
}

bool are_crossing(double v11[], double v12[], double v21[], double v22[]) {
    if ((v11[0] == v22[0] && v11[1] == v22[1]) || (v11[0] == v21[0] && v11[1] == v21[1]))
        return false;
    double cut1[3] = { v12[0] - v11[0], v12[1] - v11[1], 0.0 };
    double cut2[3] = { v22[0] - v21[0], v22[1] - v21[1], 0.0 };
    double* prod1 = new double[3];
    double* prod2 = new double[3];
    double eps = 0.00000001;
    double tmp_vec1[3] = { v21[0] - v11[0], v21[1] - v11[1], 0.0 };
    double tmp_vec2[3] = { v22[0] - v11[0], v22[1] - v11[1], 0.0 };

    prod1 = VecMulV3(cut1, tmp_vec1);
    prod2 = VecMulV3(cut1, tmp_vec2);

    double res = prod1[2] * prod2[2];
   // if (res > 0 - eps && res < 0 + eps)
       // res = 0;

    //std::cout << prod1[0] << ' ' << prod1[1] << ' ' << prod1[2] << std::endl;
    //std::cout << prod2[0] << ' ' << prod2[1] << ' ' << prod2[2] << std::endl;
    //(prod1[2] > 0 - eps && prod1[2] < 0 + eps) && (prod2[2] > 0 - eps && prod2[2] < 0 + eps)
    if (res > 0 || (prod1[2] == 0 && prod2[2] == 0))
        return false;

    tmp_vec1[0] = v11[0] - v21[0];
    tmp_vec1[1] = v11[1] - v21[1];
    tmp_vec1[2] = 0.0;
    tmp_vec2[0] = v12[0] - v21[0];
    tmp_vec2[1] = v12[1] - v21[1];
    tmp_vec2[2] = 0.0;

    prod1 = VecMulV3(cut2, tmp_vec1);
    prod2 = VecMulV3(cut2, tmp_vec2);
    //std::cout << prod1[0] << ' ' << prod1[1] << ' ' << prod1[2] << std::endl;
    //std::cout << prod2[0] << ' ' << prod2[1] << ' ' << prod2[2] << std::endl;
    
    res = prod1[2] * prod2[2];
    //if (res > 0 - eps && res < 0 + eps)
     //   res = 0;

    if (res > 0)
        return false;

    return true;
}

std::array<double, 2> InsertionPoint(NURBSSurface surface, std::vector<std::array<double, 2>> front, std::array<double, 2> a, std::array<double, 2> b, double metr, double eps) {
    std::array<double, 2> answer;
    double A[2] = { a[0], a[1] };
    double B[2] = { b[0], b[1] };
    double AB = pow((B[0] - A[0]) * (B[0] - A[0]) + (B[1] - A[1]) * (B[1] - A[1]), 0.5);
    double C[2] = { (B[0] - A[0]) * 0.5 - (B[1] - A[1]) * pow(3, 0.5) / 2 + A[0], (B[0] - A[0]) * pow(3, 0.5) / 2 + (B[1] - A[1]) * 0.5 + A[1] };
    double AB_vec[3] = { B[0] - A[0], B[1] - A[1], 0.0 };
    double AC_vec[3] = { C[0] - A[0], C[1] - A[1], 0.0 };
    double* res = new double[3];
    res = VecMulV3(AB_vec, AC_vec);
    if (res[2] < 0) {
        C[0] = (B[0] - A[0]) * 0.5 + (B[1] - A[1]) * pow(3, 0.5) / 2 + A[0];
        C[1] = -(B[0] - A[0]) * pow(3, 0.5) / 2 + (B[1] - A[1]) * 0.5 + A[1];
    }
    while ((CurveLength(eps, surface, A, C) < metr - eps) || (CurveLength(eps, surface, A, C) > metr + eps) || (CurveLength(eps, surface, B, C) < metr - eps) || (CurveLength(eps, surface, B, C) > metr + eps)) {
        double AC1[2] = { metr*(C[0] - A[0]) / CurveLength(eps, surface, A, C), metr*(C[1] - A[1]) / CurveLength(eps, surface, A, C) };
        double BC1[2] = { metr*(C[0] - B[0]) / CurveLength(eps, surface, B, C), metr*(C[1] - B[1]) / CurveLength(eps, surface, B, C) };
        double C1[2] = { A[0] + AC1[0], A[1] + AC1[1] };
        double C2[2] = { B[0] + BC1[0], B[1] + BC1[1] };
        C[0] = (C1[0] + C2[0]) / 2;
        C[1] = (C1[1] + C2[1]) / 2;
    }

    if (C[0] > 0 - eps && C[0] < 0 + eps)
        C[0] = 0;
    if (C[1] > 0 - eps && C[1] < 0 + eps)
        C[1] = 0;
    answer[0] = C[0];
    answer[1] = C[1];

    for (int i=0; i < front.size() - 1; i++) {
        double left[2] = { front[i][0], front[i][1] };
        double right[2] = { front[i + 1][0], front[i + 1][1] };
        if ((answer == front[i] && (a == front[i + 1] || b == front[i + 1])) ||
            (answer == front[i + 1] && (a == front[i] || b == front[i]))) {
            break;
        }
        if (are_crossing(A, C, left, right) || are_crossing(B, C, left, right)) {
            answer[0] = -1.0;
            answer[1] = -1.0;
            return answer;
        }
        if (CurveLength(eps, surface, C, left) < metr / 1.85 || CurveLength(eps, surface, C, right) < metr / 1.85) {
            answer[0] = -1.0;
            answer[1] = -1.0;
            return answer;
        }
    }

    return answer;
}

bool oneTriangleCheck(std::vector<std::array<double, 2>> front, std::array<double, 2> a, std::array<double, 2> b, std::array<double, 2> c) {
    int count = 0;
    std::pair<std::array<double, 2>, std::array<double, 2>> side1(b, c);
    std::pair<std::array<double, 2>, std::array<double, 2>> side2(c, a);
    for (int i = 0; i < front.size() - 1; i++) {
        std::pair<std::array<double, 2>, std::array<double, 2>> temp(front[i], front[i + 1]);
        if (temp == side1) { count++; }
        if (temp == side2) { count++; }
        if (count > 1) {
            return true;
        }
    }
    return false;
}

std::array<double, 2> gradientDescent(NURBSSurface surface, double* p, std::array<double, 2> start_p) {
    double step = 0.00001;
    double eps = 0.01;
    double r = 1;
    std::cout << std::endl << p[0] << " " << p[1] << " " << p[2];
    while (r > eps || r < -eps) {
        double l1[3], l2[3], cur_p[3];
        //double l2[3];
        //double cur_p[3];
        surface.getPoint(start_p[0], start_p[1], cur_p);
        if (cur_p[0] == p[0] && cur_p[1] == p[1] && cur_p[2] == p[2]) { break; }
        surface.getDer1Point(start_p[0], start_p[1], l1);
        surface.getDer2Point(start_p[0], start_p[1], l2);
        double der1 = ((cur_p[0] - p[0]) * l1[0] + (cur_p[1] - p[1]) * l1[1] + (cur_p[2] - p[2]) * l1[2]) / pow(pow(cur_p[0] - p[0], 2) + pow(cur_p[1] - p[1], 2) + pow(cur_p[2] - p[2], 2), 0.5);
        double der2 = ((cur_p[0] - p[0]) * l2[0] + (cur_p[1] - p[1]) * l2[1] + (cur_p[2] - p[2]) * l2[2]) / pow(pow(cur_p[0] - p[0], 2) + pow(cur_p[1] - p[1], 2) + pow(cur_p[2] - p[2], 2), 0.5);
        start_p[0] = start_p[0] - step * der1;
        start_p[1] = start_p[1] - step * der2;
        surface.getPoint(start_p[0], start_p[1], cur_p);
        r = pow(pow(cur_p[0] - p[0], 2) + pow(cur_p[1] - p[1], 2) + pow(cur_p[2] - p[2], 2), 0.5);
    }
    return start_p;
}

std::vector<std::array<double, 2>> discreteBoundary(NURBSSurface surface, NURBSCurve curve, bool flag, std::array<double, 2> prev) {
    double start;
    double end;
    double step;
    std::vector<std::array<double, 2>> answer;
    if (!flag) {
        start = curve.getBegin();
        end = curve.getEnd();
        step = 0.05;
    }
    else {
        start = curve.getEnd();
        end = curve.getBegin();
        step = -0.05;
    }
    for (double i = start; ((i > end + step) && (i < start - step)) || ((i > start - step) && (i < end + step)); i += step) {
        double* a = new double[3];
        curve.getPoint(i, a);
        prev = gradientDescent(surface, a, prev);
        answer.push_back(prev); 
    }
    answer.pop_back();
    return answer;
}

bool sharpAngle(std::array<double, 2> a, std::array<double, 2> b, std::array<double, 2> c) {
    std::array<double, 2> v1 = { a[0] - b[0], a[1] - b[1] };
    std::array<double, 2> v2 = { c[0] - b[0], c[1] - b[1] };
    double len1 = pow(v1[0] * v1[0] + v1[1] * v1[1], 0.5);
    double len2 = pow(v2[0] * v2[0] + v2[1] * v2[1], 0.5);
    double cosinus = (v1[0] * v2[0] + v1[1] * v2[1]) / (len1 * len2);
    double angle = acos(cosinus) * 180 / 3.14159265;
    if (angle <= 145) { return true; }
    return false;
}

void unitDiscretization(std::vector<std::array<double, 2>>& points, NURBSSurface surface, std::vector<std::vector<std::array<double, 2>>> discr_bound, double metr, double eps, bool flag) {
    std::vector<std::array<double, 2>> boundary;
    for (auto i : discr_bound) {
        for (auto j : i) {
            boundary.push_back(j);
        }
    }
    std::array<double, 2> A = boundary[0];
    std::array<double, 2> B = boundary[1];
    std::array<double, 2> var;
    double sigma = 0;
    double r = 0;
    double n = 0;
    double l1[3];
    double l2[3];
    auto pointer = boundary.begin() + 2;
    points.push_back(A);
    while (pointer - 1 != boundary.end()) {
        double ta[2] = { A[0], A[1] };
        double tb[2] = { B[0], B[1] };
        while (CurveLength(eps, surface, ta, tb) < metr && pointer != boundary.end()) {
            if (pointer == boundary.end() || sharpAngle(A, B, *pointer)) { break; }
            B = *pointer;
            pointer++;
            ta[0] = A[0];
            ta[1] = A[1];
            tb[0] = B[0];
            tb[1] = B[1];
        }
        n += 1;
        double unit_vector[2] = { (B[0] - A[0]) / pow((B[0] - A[0]) * (B[0] - A[0]) + (B[1] - A[1]) * (B[1] - A[1]), 0.5), (B[1] - A[1]) / pow((B[0] - A[0]) * (B[0] - A[0]) + (B[1] - A[1]) * (B[1] - A[1]), 0.5) };
        surface.getDer1Point(A[0], A[1], l1);
        surface.getDer2Point(A[0], A[1], l2);
        double metric_A[2][2] = { {ScalarMulV3(l1,l1), ScalarMulV3(l1, l2)}, {ScalarMulV3(l2, l1), ScalarMulV3(l2, l2) } };
        double s = metr / pow((unit_vector[0] * metric_A[0][0] + unit_vector[1] * metric_A[1][0]) * unit_vector[0] + (unit_vector[0] * metric_A[0][1] + unit_vector[1] * metric_A[1][1]) * unit_vector[1], 0.5);
        while (1) {
            var[0] = A[0] + s * unit_vector[0];
            var[1] = A[1] + s * unit_vector[1];
            double ta[2] = { A[0], A[1] };
            double tb[2] = { var[0], var[1] };
            r = CurveLength(eps, surface, ta, tb);
            s = 2 * s / (metr + r);
            if ((metr - r) < eps) { break; }
        }
        ta[0] = A[0];
        ta[1] = A[1];
        tb[0] = B[0];
        tb[1] = B[1];
        if (CurveLength(eps, surface, ta, tb) < r) {
            A[0] = B[0];
            A[1] = B[1];
            points.push_back(A);
            if (pointer == boundary.end()) { break; }
            B = *pointer;
            pointer++;
        }
        else {
            A[0] = var[0];
            A[1] = var[1];
            points.push_back(A);
        }
    }
    if (flag && ((points[points.size() - 1][0] > boundary[boundary.size() - 1][0]) || (points[points.size() - 1][1] > boundary[boundary.size() - 1][1]))) {
        double ta[2] = { boundary[boundary.size() - 1][0], boundary[boundary.size() - 1][1] };
        double tb[2] = { points[points.size() - 1][0], points[points.size() - 1][1] };
        sigma = CurveLength(eps, surface, ta, tb) / n;
        if (sigma > 0) {
            while (n > 0) {
                points.pop_back();
                n--;
            }
            unitDiscretization(points, surface, discr_bound, metr, eps, false);
        }
    }
}

std::array<double, 2> averageAngle(double* A, double* B, double* C) {
    double min;
    std::array<double, 3> v1 = { B[0] - A[0], B[1] - A[1], B[2] - A[2] };
    std::array<double, 3> v2 = { C[0] - A[0], C[1] - A[1], C[2] - A[2] };
    double len1 = pow(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2], 0.5);
    double len2 = pow(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2], 0.5);
    double cosinus = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (len1 * len2);
    double angle1 = acos(cosinus) * 180 / 3.14159265;

    v1[0] = A[0] - B[0];
    v1[1] = A[1] - B[1];
    v1[2] = A[2] - B[2];
    v2[0] = C[0] - B[0];
    v2[1] = C[1] - B[1];
    v2[2] = C[2] - B[2];
    len1 = pow(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2], 0.5);
    len2 = pow(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2], 0.5);
    cosinus = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (len1 * len2);
    double angle2 = acos(cosinus) * 180 / 3.14159265;

    v1[0] = A[0] - C[0];
    v1[1] = A[1] - C[1];
    v1[2] = A[2] - C[2];
    v2[0] = B[0] - C[0];
    v2[1] = B[1] - C[1];
    v2[2] = B[2] - C[2];
    len1 = pow(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2], 0.5);
    len2 = pow(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2], 0.5);
    cosinus = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (len1 * len2);
    double angle3 = acos(cosinus) * 180 / 3.14159265;

    min = angle1;
    if (angle2 < min) { min = angle2; }
    if (angle3 < min) { min = angle3; }
    std::array<double, 2> answer = { (angle1 + angle2 + angle3) / 3, min };
    return answer;
}

std::array<double, 6> qualityEstimate(NURBSSurface surface, std::vector<std::array<double, 2>> coords, std::map<int, std::array<int, 3>> mesh, double metr, std::vector<double>& stat_alpha, std::vector<double>& stat_sigma, std::vector<double>& stat_angle) {
    double min_sigma = 2;
    double min_alpha = 2;
    double avg_sigma = 0;
    double avg_alpha = 0;
    double min_angle = 180;
    double avg_angle = 0;
    for (auto triangle : mesh) {
        double A[2] = { coords[triangle.second[0]][0], coords[triangle.second[0]][1] };
        double B[2] = { coords[triangle.second[1]][0], coords[triangle.second[1]][1] };
        double C[2] = { coords[triangle.second[2]][0], coords[triangle.second[2]][1] };
        double p1[3];
        double p2[3];
        double p3[3];
        surface.getPoint(A[0], A[1], p1);
        surface.getPoint(B[0], B[1], p2);
        surface.getPoint(C[0], C[1], p3);
        auto angles = averageAngle(p1, p2, p3);
        if (angles[1] < min_angle) { min_angle = angles[1]; }
        avg_angle += angles[0];
        stat_angle.push_back(angles[1]);
        double AB[3] = { p2[0] - p1[0],p2[1] - p1[1],p2[2] - p1[2] };
        double AC[3] = { p3[0] - p1[0],p3[1] - p1[1],p3[2] - p1[2] };
        double BC[3] = { p3[0] - p2[0],p3[1] - p2[1],p3[2] - p2[2] };
        double* vec = VecMulV3(AB, AC);
        double len_vec = pow(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2], 0.5);
        double len_AB = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
        double len_AC = AC[0] * AC[0] + AC[1] * AC[1] + AC[2] * AC[2];
        double len_BC = BC[0] * BC[0] + BC[1] * BC[1] + BC[2] * BC[2];
        double alpha = 2 * pow(3, 0.5) * len_vec / (len_AB + len_AC + len_BC);
        if (alpha < min_alpha) { min_alpha = alpha; }
        avg_alpha += alpha;
        stat_alpha.push_back(alpha);
        auto r1 = CurveLength(0.0001, surface, A, B) / 1;
        auto r2 = CurveLength(0.0001, surface, A, C) / 1;
        auto r3 = CurveLength(0.0001, surface, B, C) / 1;
        if (1 / r1 < r1) { r1 = 1 / r1; }
        if (1 / r2 < r2) { r2 = 1 / r2; }
        if (1 / r3 < r3) { r3 = 1 / r3; }
        r1 = r1 * r2 * r3;
        if (r1 < min_sigma) { min_sigma = r1; }
        avg_sigma = avg_sigma + r1;
        stat_sigma.push_back(r1);
    }
    std::array<double, 6> answer = { avg_alpha / mesh.size(), min_alpha, avg_sigma / mesh.size(), min_sigma, avg_angle / mesh.size(), min_angle };
    return answer;
}

void SubDivideBooudary(IGESData* data) {
    auto sbd = data->getSurfacesBoundaryData();
    auto surfaces = data->getNURBSSurfaces();
    auto curves = data->getNURBSCurves();
    std::vector<std::array<double, 9>> res_coords;
    std::vector<std::array<int, 4>> res_mesh;
    int amount_coords = 0;
    int amount_mesh = 0;
    int iteration_step = 0;
    int porog = 0;
    double abs_sigma = 0;
    double abs_alpha = 0;
    double abs_angle = 0;
    double min_sigma = 2;
    double min_alpha = 2;
    double min_angle = 180;
    std::vector<double> stat_alpha;
    std::vector<double> stat_sigma;
    std::vector<double> stat_angle;
    for (auto i = sbd.begin(); i != sbd.end(); i++) {
        std::cout << "---------------------------------------------------------------------" << std::endl;
        auto surface_id = (*i).first;
        auto surface = surfaces.at(surface_id);
        std::vector<std::pair<int, bool>> bound_curves;
        std::vector<std::array<double, 2>> points;
        std::vector<std::vector<std::array<double, 2>>> discr_bound;
        double u_0 = 0;
        double u_1 = 1;
        double v_0 = 0;
        double v_1 = 1;
        double metr = 1;
        double eps = 0.001;
        std::array<double, 2> start_point = { 0.5, 0.5 };
        auto bound = sbd.at(surface_id);
        auto p = *bound.begin();
        for (auto item : bound.at(p.first)) {
            std::pair<int, bool> a(item.first, item.second);
            bound_curves.push_back(a);
        }
        for (auto item : bound_curves) {
            auto temp = discreteBoundary(surface, curves.at(item.first), item.second, start_point);
            discr_bound.push_back(temp);
            start_point = temp[temp.size() - 1];
        }
        std::vector<std::array<double, 2>> front;
        DelaunayTriangulation D_triang = DelaunayTriangulation(u_0, u_1, v_0, v_1, surface, eps, metr);
        unitDiscretization(points, surface, discr_bound, metr, eps, true);
        for (int i = 0; i < points.size(); i++) {
            if (points[i][0] < 0.000000001 && points[i][0] > -0.000000001) { points[i][0] = 0; }
            if (points[i][1] < 0.000000001 && points[i][1] > -0.000000001) { points[i][1] = 0; }
            if (points[i][0] - 1 < 0.000000001 && points[i][0] - 1 > -0.000000001) { points[i][0] = 1; }
            if (points[i][1] - 1 < 0.000000001 && points[i][1] - 1 > -0.000000001) { points[i][1] = 1; }
            std::cout << "(" << points[i][0] << ", " << points[i][1] << ")," << std::endl;
            //std::array<double, 2> p = { points[i], points[i + 1] };
            D_triang.addPoint(points[i], true);
            front.push_back(points[i]);
        }
        std::cout << std::endl << std::endl << std::endl;
        D_triang.boundaryRecovery(front, true);
        auto k = D_triang.exportTriangles();
        for (auto t : k)
            std::cout << "(" << t.second[0] << ", " << t.second[1] << ", " << t.second[2] << ")," << std::endl;
        std::cout << std::endl << std::endl << std::endl;
        //front.push_back(front[0]);
        /*std::array<std::vector<std::array<double, 2>>, 4> bound_lines;
        for (auto point : front) {
            if (point[1] - 1 > -eps && point[1] - 1 < eps) { bound_lines[0].push_back(point); }
            if (point[0] > -eps && point[0] < eps) { bound_lines[1].push_back(point); }
            if (point[1] > -eps && point[1] < eps) { bound_lines[2].push_back(point); }
            if (point[0] - 1 > -eps && point[0] - 1 < eps) { bound_lines[3].push_back(point); }
        }
        auto temporary = bound_lines[3][0];
        bound_lines[3].erase(bound_lines[3].begin());
        bound_lines[3].push_back(temporary);*/
        std::array<double, 2> insertion_point = { -1.0, -1.0 };
        int index = front.size();
        bool flag = false;
        bool delition_flag = false;
        D_triang.recalculateCircumcentre();
        int counter = 0;
        while (front.size() > 1) {
            //index = front.size();
            counter++;
            std::cout << counter << std::endl;
            if (counter == 53) {
                std::cout << std::endl;
            }
            while (front.size() > 0 && D_triang.actualPoints(front[0])) {
                front.erase(front.begin());
                delition_flag = true;
            }
            if (front.size() > 0 && D_triang.actualPoints(front[front.size() - 1])) {
                delition_flag = true;
            }
            if (front.size() < 2) { break; }
            if (delition_flag) {
                front.erase(front.end() - 1);
                auto t = front[0];
                front.erase(front.begin());
                front.push_back(t);
                delition_flag = false;
            }
            index = front.size();
            if (front.size() < 2) { break; }
            //std::cout<<index<<std::endl;
            insertion_point = InsertionPoint(surface, front, front[index - 1], front[index % front.size()], metr, eps);
            //if (front.size() == 3) {
                //D_triang.lastTriangle(front[2], front[0], front[1]);
               // break;
            //}
            for (auto item : front) {
                if (insertion_point[0] != -1.0 && insertion_point[1] != -1.0 && (item[0] - eps) < insertion_point[0] && (item[0] + eps) > insertion_point[0] && (item[1] - eps) < insertion_point[1] && (item[1] + eps) > insertion_point[1]) {
                     insertion_point[0] = -1.0;
                     insertion_point[1] = -1.0;
                     break;
                    }
            }
            //index++;
            if (insertion_point[0] == -1.0 && insertion_point[1] == -1.0) {
                insertion_point = D_triang.selectFromDT(front[index - 1], front[index % front.size()]);
                if (insertion_point[0] == 0.086475275566472526 && insertion_point[1] == 0.94230769230769229) {
                    std::cout << std::endl;
                }
            }
            else {
                D_triang.addPoint(insertion_point, false);
                D_triang.deleteBaseTriangle(front[index - 1], front[index % front.size()], insertion_point);
            }
            if (oneTriangleCheck(front, front[index - 1], front[index % front.size()], insertion_point)) {
                if (front.size() > 3 && front[front.size() - 1] == front[2]) {
                    front.erase(front.begin());
                    front.erase(front.begin());
                    front.erase(front.begin());
                }
                else {
                    front.erase(front.end() - 1);
                    front.erase(front.begin());
                    front.erase(front.begin());
                }
            }
            else if (front[index - 2] == insertion_point) {
                front.erase(front.begin() + (index - 1));
            }
            else if(front[(index % front.size()) + 1] == insertion_point) {
                front.erase(front.begin());
            }
            else {
                front.insert(front.end(), insertion_point);
            }
            insertion_point[0] = -1.0;
            insertion_point[1] = -1.0;
        }
        //D_triang.boundaryRecovery(front, false);
        auto coords = D_triang.getCoords();
        auto mesh = D_triang.getMesh();
        std::cout << std::endl << std::endl;
        for (auto i : coords) {
            std::cout << "(" << i[0] << ",  " << i[1] << ")," << std::endl;
        }
        for (auto i : mesh) {
            std::cout << "(" << i.second[0] << ", " << i.second[1] << ", " << i.second[2] << ")," << std::endl;
        }
        std::cout << std::endl << std::endl;
        auto quality = qualityEstimate(surface, coords, mesh, metr, stat_alpha, stat_sigma, stat_angle);
        abs_alpha += quality[0];
        abs_sigma += quality[2];
        abs_angle += quality[4];
        if (min_alpha > quality[1]) { min_alpha = quality[1]; }
        if (min_sigma > quality[3]) { min_sigma = quality[3]; }
        if (min_angle > quality[5]) { min_angle = quality[5]; }
        coords.erase(coords.begin());
        coords.erase(coords.begin());
        coords.erase(coords.begin());
        coords.erase(coords.begin());
        amount_coords += coords.size();
        amount_mesh += mesh.size();
        iteration_step += 1;
        std::cout << std::endl << std::endl << amount_coords << std::endl;
        for (auto i : coords) {
            double point[3];
            double* norm = new double[3];
            double* l1 = new double[3];
            double* l2 = new double[3];
            surface.getPoint(i[0], i[1], point);
            surface.getDer1Point(i[0], i[1], l1);
            surface.getDer2Point(i[0], i[1], l2);
            norm = VecMulV3(l1, l2);
            NormV3(norm);
            NormV3(l1);
            std::array<double, 9> temp_1 = { point[0],point[1],point[2] , norm[0],norm[1],norm[2], l1[0],l1[1],l1[2] };
            res_coords.push_back(temp_1);
            std::cout << " " << point[0] << "  " << point[1] << "  " << point[2]  << "  " << norm[0] << "  " << norm[1] << "  " << norm[2] << "  " << l1[0] << "  " << l1[1] << "  " << l1[2] << std::endl;
        }
        std::cout << 0 << std::endl;
        std::cout << amount_mesh << std::endl;
        for (auto i : mesh) {
            std::array<int, 4> temp_2 = { iteration_step, i.second[0] - 3 + porog ,i.second[1] - 3 + porog ,i.second[2] - 3 + porog };
            res_mesh.push_back(temp_2);
            std::cout << "   " << iteration_step << "         " << i.second[0] - 3 + porog << "     " << i.second[1] - 3 + porog << "     " << i.second[2] - 3 + porog << std::endl;
        }
        std::cout << std::endl;
        porog += coords.size();
    }
    abs_alpha = abs_alpha / iteration_step;
    abs_sigma = abs_sigma / iteration_step;
    abs_angle = abs_angle / iteration_step;
    std::cout << std::endl << std::endl << std::endl;
    for (auto i : stat_alpha) {
        std::cout << i << "," << std::endl;
    }
    std::cout << std::endl << std::endl << std::endl;
    for (auto i : stat_sigma) {
        std::cout << i << "," << std::endl;
    }
    std::cout << std::endl << std::endl << std::endl;
    for (auto i : stat_angle) {
        std::cout << i << "," << std::endl;
    }
    std::cout << std::endl << std::endl << std::endl;
    std::ofstream meshlog("meshlog.mesh", std::ios::out);
    //число вершин
    meshlog << res_coords.size() << std::endl;
    //цикл записи вершин
    for (auto i = res_coords.begin(); i < res_coords.end(); ++i) {
        meshlog << std::get<0>(*i) << "    " << std::get<1>(*i) << "    " << std::get<2>(*i) << std::endl; // << "    " << std::get<3>(*i)<< "    " << std::get<4>(*i) << "    " << std::get<5>(*i) << "    " << std::get<6>(*i) << "    " << std::get<7>(*i) << "    " << std::get<8>(*i) << std::endl;
    }
    //число внутренних  Ё
    meshlog << 0 << std::endl;
    //число  Ё
    meshlog << res_mesh.size() << std::endl;
    //цикл записи индексов  Ё
    for (auto i = res_mesh.begin(); i < res_mesh.end(); ++i) {
        meshlog << std::get<0>(*i) << "    " << std::get<1>(*i) << "    " << std::get<2>(*i) << "    " << std::get<3>(*i) << std::endl;
    }
    /*meshlog << 126 << std::endl;
    for (int i = 1; i < 17; i++) { meshlog << 1 << "    " << i << "    " << i + 1 << std::endl; }
    for (int i = 17; i < 32; i++) { meshlog << 2 << "    " << i << "    " << i + 1 << std::endl; }
    for (int i = 32; i < 48; i++) { meshlog << 3 << "    " << i << "    " << i + 1 << std::endl; }
    for (int i = 48; i < 63; i++) { meshlog << 4 << "    " << i << "    " << i + 1 << std::endl; }
    meshlog << 4 << "    " << 63 << "    " << 1 << std::endl;
    for (int i = 323; i < 339; i++) { meshlog << 5 << "    " << i << "    " << i + 1 << std::endl; }
    for (int i = 339; i < 354; i++) { meshlog << 6 << "    " << i << "    " << i + 1 << std::endl; }
    for (int i = 354; i < 370; i++) { meshlog << 7 << "    " << i << "    " << i + 1 << std::endl; }
    for (int i = 370; i < 385; i++) { meshlog << 8 << "    " << i << "    " << i + 1 << std::endl; }
    meshlog << 8 << "    " << 385 << "    " << 323 << std::endl;*/
    //закрытие потока
    meshlog.close();
}

int main()
{
    std::cout << "Hello World!\n";
    //std::string file_name = "Stift.IGS";
    //std::string file_name = "SphereGit.IGS";
    std::string file_name = "hard3.IGS";
    IGESLoader Loader;
    Loader.setGeoFileName(file_name);
    IGESData data;
    Loader.loadGeo();
    data = Loader.GetIGESModel(); 
    SubDivideBooudary(&data);
    auto c = data.getNURBSCurves();
    double res[3];
    const auto& s = data.getNURBSSurfaces();
    auto sbd = data.getSurfacesBoundaryData();
    c[4].getPoint(0.89, res);
    std::cout << res[0] << ' ' << res[1] << ' ' << res[2] << std::endl;
    s.at(46).getPoint(0.89, 0, res);
    std::cout << res[0] << ' ' << res[1] << ' ' << res[2] << std::endl;
    c[36].getPoint(0.89, res);
    std::cout << res[0] << ' ' << res[1] << ' ' << res[2] << std::endl;
    s.at(46).getPoint(1, 0.89, res);
    std::cout << res[0] << ' ' << res[1] << ' ' << res[2] << std::endl;
    c[26].getPoint(0.89, res);
    std::cout << res[0] << ' ' << res[1] << ' ' << res[2] << std::endl;
    s.at(46).getPoint(0.11, 1, res);
    std::cout << res[0] << ' ' << res[1] << ' ' << res[2] << std::endl;
    c[16].getPoint(0.89, res);
    std::cout << res[0] << ' ' << res[1] << ' ' << res[2] << std::endl;
    s.at(46).getPoint(0, 0.11, res);
    std::cout << res[0] << ' ' << res[1] << ' ' << res[2] << std::endl;
    //double res[18];
    /*“есты скал€рного и векторного произведени€*/
    //std::vector <double> l1(3, 1);
    //std::vector <double> l2(3, 3);
    /*double l1[3] = { 0, 1, 0 };
    double l2[3] = { 0, 0, -1 };
    double l3[3] = { 0, 0, 0 };*/

    PrincipalCurvatures(&data, 0, 0);
    PrincipalCurvatures(&data, 1.56, 3.14);
    PrincipalCurvatures(&data, 1, 0);
    PrincipalCurvatures(&data, -1, 0);
    PrincipalCurvatures(&data, -1, 1);
    //double l3[3] = { -1,0,0 };
    //std::vector <double> l3(3);
    ///l3 = VecMulV3(l1, l2);
    //std::cout << ScalarMulV3(l1, l2);
    //l3 = VecMulV3(l1, l2);
    //sur.second.
    //auto bd = data.getSurfacesBoundaryData();
    std::cout << "Hello World!\n";
}  

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
