#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

// Line class
class Line {
public:
    float gradient;
    float c_value;

    Line(float g, float c);
    Line();
};

// Point class
class Point {
public:
    static Point ref;
    float x, y;

    Point(float x_, float y_);
    Point();

    float segment() const;
    Point refers_to() const;
    float dot(const Point &p2) const;

    friend ostream &operator<<(ostream &os, const Point &obj);
};

// Utility functions
bool fileExists(const string &filename);
void handleFileError(const string &filename);
void loadPoints(const string &filename, vector<Point> &points);
int countLines(const string &filename);

bool point_comparator(Point a, Point b);
vector<Point> point_sorter(vector<Point> points, Point avg);

Point midpoint_func(Point A, Point B);
Point closest_point_to_center(const vector<Point> &points);
double gradient_func(const Point &p1, const Point &p2);
float perpendicular_gradient_func(float m);
void bucketSort(vector<Point>& points, int bucket_count);
Line eol_g_p(Point A, float gradient);
Line eol_pp(Point A, Point B);

bool similarity_check(Point A, Point B);
double calculate_error_margin(int m, double c_max);
double calculate_c_max(const vector<Point> &points);

bool affine_dependent_check(vector<Point> &A);
Point new_coordinate_func(Line A, Line B);
Line perpendicular_bisector_func(Point A, Point B);

#ifdef DEBUG
void run_tests();
#endif
