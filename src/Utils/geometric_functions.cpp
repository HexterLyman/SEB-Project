#include "../Header/geometric_functions.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <bits/stdc++.h>

using namespace std;

// Line class implementation
Line::Line(float g, float c) : gradient(g), c_value(c) {}
Line::Line() : gradient(0.0), c_value(0.0) {}


// Point class implementation
Point Point::ref = Point(0, 0);

Point::Point(float x_, float y_) : x(x_), y(y_) {}
Point::Point() : x(ref.x), y(ref.y) {}


float Point::segment() const {
    return sqrt((x - ref.x) * (x - ref.x) + (y - ref.y) * (y - ref.y));
}

Point Point::refers_to() const {
    return Point(x - ref.x, y - ref.y);
}

float Point::dot(const Point &p2) const {
    return x * p2.x + y * p2.y;
}

ostream &operator<<(ostream &os, const Point &obj) {
    os << "(" << obj.x << "," << obj.y << ")";
    return os;
}

ostream& operator<<(ostream& os, const vector<Point>& points) {
    os << "[";
    for (size_t i = 0; i < points.size(); ++i) {
        os << points[i];  // Use the existing operator<< for Point
        if (i != points.size() - 1) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}

// Utility functions
bool fileExists(const string &filename) {
    ifstream file(filename);
    return file.good();
}

void handleFileError(const string &filename) {
    cerr << "Error: Unable to open file '" << filename << "'. Please check if the file exists and is accessible." << endl;
}

void loadPoints(const string &filename, vector<Point> &points) {
    if (!fileExists(filename)) {
        handleFileError(filename);
        return;
    }

    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        line.erase(remove(line.begin(), line.end(), '{'), line.end());
        line.erase(remove(line.begin(), line.end(), '}'), line.end());
        line.erase(remove(line.begin(), line.end(), ','), line.end());

        istringstream iss(line);
        float x, y;
        if (iss >> x >> y) {
            points.emplace_back(x, y);
        } else {
            cerr << "Warning: Skipping malformed line: " << line << endl;
        }
    }
}

int countLines(const string &filename) {
    if (!fileExists(filename)) {
        handleFileError(filename);
        return -1;
    }

    ifstream file(filename);
    return count(istreambuf_iterator<char>(file), istreambuf_iterator<char>(), '\n');
}

// Additional function implementations
bool point_comparator(Point a, Point b) {
    return a.segment() > b.segment();
}



vector<Point> point_sorter(vector<Point> points, Point avg) {
    Point::ref = avg;

    float max_distance = 0;
    for (const Point &point : points) {
        max_distance = max(max_distance, point.segment());
    }

    const int bucket_count = 1000;
    vector<vector<Point>> buckets(bucket_count);

    for (const Point &point : points) {
        int bucket_index = static_cast<int>((point.segment() / max_distance) * (bucket_count - 1));
        buckets[bucket_index].push_back(point);
    }

    vector<Point> sorted_points;
    for (int i = bucket_count - 1; i >= 0; --i) {
        for (const Point &point : buckets[i]) {
            sorted_points.push_back(point);
        }
    }

    return sorted_points;
}

Point midpoint_func(Point A, Point B) {
    return Point((A.x + B.x) / 2, (A.y + B.y) / 2);
}

Point closest_point_to_center(const vector<Point> &points) {
    Point avg_point;
    size_t num_points = points.size();

    for (const auto &point : points) {
        avg_point.x += point.x;
        avg_point.y += point.y;
    }

    avg_point.x /= num_points;
    avg_point.y /= num_points;

    return avg_point;
}

double gradient_func(const Point &p1, const Point &p2) {
    if (p1.x == p2.x) {
        return numeric_limits<double>::infinity();
    }
    return (p2.y - p1.y) / (p2.x - p1.x);
}

float perpendicular_gradient_func(float m) {
    if (abs(m) < 1e-9) {
        return 0.0;
    }
    return -1.0 / m;
}

Line eol_g_p(Point A, float gradient) {
    Line output;
    float new_gradient = perpendicular_gradient_func(gradient);
    float c_gp = A.y - (new_gradient * A.x);

    output.gradient = new_gradient;
    output.c_value = c_gp;

    return output;
}

Line eol_pp(Point A, Point B) {
    Line output;
    float m_pp = gradient_func(A, B);
    float c_pp = (-m_pp * A.x) + A.y;

    output.gradient = m_pp;
    output.c_value = c_pp;

    return output;
}

bool similarity_check(Point A, Point B) {
    return abs(A.x - B.x) < 0.0001 && abs(A.y - B.y) < 0.0001;
}

double calculate_error_margin(int m, double c_max) {
    double term1 = 1.0 / (2 * pow(2, 2 + 2));
    double term2 = 1.0 / (25 * pow(m, 2) * pow(2 * c_max, 8 * m + 12));
    return max(term1, term2);
}

double calculate_c_max(const vector<Point> &points) {
    double c_max = 0;
    for (const auto &point : points) {
        c_max = std::max(c_max, static_cast<double>(abs(point.x)));
        c_max = std::max(c_max, static_cast<double>(abs(point.y)));
    }
    return c_max;
}


bool affine_dependent_check(vector<Point> &A) {
    if (A.size() < 3) return false;

    Line affine_checker = eol_pp(A[0], A[1]);
    int m = 2;
    double c_max = calculate_c_max(A);
    double error_margin = calculate_error_margin(m, c_max);

    double y_pred = affine_checker.gradient * A[2].x + affine_checker.c_value;
    double y_actual = A[2].y;
    return abs(y_pred - y_actual) < error_margin;
}

Point new_coordinate_func(Line A, Line B) {
    float x = (B.c_value - A.c_value) / (A.gradient - B.gradient);
    float y = B.gradient * x + B.c_value;

    return Point(x, y);
}

Line perpendicular_bisector_func(Point A, Point B) {
    Point mp = midpoint_func(A, B);
    float gradient = gradient_func(A, B);
    return eol_g_p(mp, gradient);
}

#ifdef DEBUG
#include <cassert>
void run_tests() {
    // Test perpendicular_gradient_func
    assert(abs(perpendicular_gradient_func(1) + 1.0) < 1e-9);
    assert(abs(perpendicular_gradient_func(-1) - 1.0) < 1e-9);
    assert(abs(perpendicular_gradient_func(0.5) + 2.0) < 1e-9);
    assert(abs(perpendicular_gradient_func(-0.5) - 2.0) < 1e-9);
    assert(abs(perpendicular_gradient_func(0) - 0.0) < 1e-9);

    // Test gradient_func
    assert(abs(gradient_func(Point(0, 0), Point(1, 1)) - 1.0) < 1e-9);
    assert(abs(gradient_func(Point(0, 0), Point(-1, -1)) - 1.0) < 1e-9);
    assert(abs(gradient_func(Point(0, 0), Point(1, 0)) - 0.0) < 1e-9); // Horizontal line
   assert(std::isinf(gradient_func(Point(0, 0), Point(0, 1))) && "Gradient should be infinity for vertical lines."); // Near-vertical line approximation

    // Test midpoint_func
    Point mid = midpoint_func(Point(0, 0), Point(2, 2));
    assert(abs(mid.x - 1.0) < 1e-9 && abs(mid.y - 1.0) < 1e-9);

    // Test new_coordinate_func (intersection of two lines)
    Line line1(1, 0);  // y = x
    Line line2(-1, 2); // y = -x + 2
    Point intersection = new_coordinate_func(line1, line2);
    assert(abs(intersection.x - 1.0) < 1e-9 && abs(intersection.y - 1.0) < 1e-9);

    // Test similarity_check
    assert(similarity_check(Point(0, 0), Point(0, 0)));
    assert(!similarity_check(Point(0, 0), Point(1, 1)));

    // Test point_sorter
    vector<Point> points = {Point(0, 0), Point(1, 1), Point(2, 2)};
    Point ref_point(0, 0);
    vector<Point> sorted = point_sorter(points, ref_point);
    assert(sorted[0].x == 2 && sorted[0].y == 2); // Farthest point first

    // Test affine_dependent_check
    vector<Point> collinear_points = {Point(0, 0), Point(1, 1), Point(2, 2)};
    assert(affine_dependent_check(collinear_points));
    vector<Point> non_collinear_points = {Point(0, 0), Point(1, 1), Point(1, 2)};
    assert(!affine_dependent_check(non_collinear_points));

    // Test perpendicular_bisector_func
    Line perp_bisector = perpendicular_bisector_func(Point(0, 0), Point(2, 0));
    assert(abs(perp_bisector.gradient - 0.0) < 1e-9); // Should be vertical

    std::cout << "All tests passed!" << std::endl;
}
#endif
