#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <sstream>
#include <iomanip>

// Function to decode a value from a given base to decimal
long long decode_value(int base, const std::string &value) {
    return std::stoll(value, nullptr, base);
}

// Lagrange interpolation to find the polynomial value at x = 0
double lagrange_interpolation(const std::vector<std::pair<int, double>> &points, int x) {
    double total = 0.0;
    int n = points.size();
    for (int i = 0; i < n; ++i) {
        double xi = points[i].first;
        double yi = points[i].second;
        double term = yi;
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                double xj = points[j].first;
                term *= (x - xj) / (xi - xj);
            }
        }
        total += term;
    }
    return total;
}

// Vandermonde method to find the polynomial coefficients
double vandermonde_method(const std::vector<std::pair<int, double>> &points) {
    int n = points.size();
    std::vector<std::vector<double>> V(n, std::vector<double>(n));
    std::vector<double> Y(n);

    // Fill the Vandermonde matrix
    for (int i = 0; i < n; ++i) {
        Y[i] = points[i].second;
        for (int j = 0; j < n; ++j) {
            V[i][j] = pow(points[i].first, j);
        }
    }

    // Solving V * C = Y for coefficients C using Gaussian elimination
    for (int i = 0; i < n; ++i) {
        // Partial pivoting
        for (int j = i + 1; j < n; ++j) {
            if (fabs(V[j][i]) > fabs(V[i][i])) {
                std::swap(V[i], V[j]);
                std::swap(Y[i], Y[j]);
            }
        }
        for (int j = i + 1; j < n; ++j) {
            double ratio = V[j][i] / V[i][i];
            for (int k = i; k < n; ++k) {
                V[j][k] -= ratio * V[i][k];
            }
            Y[j] -= ratio * Y[i];
        }
    }

    // Back substitution
    std::vector<double> coefficients(n);
    for (int i = n - 1; i >= 0; --i) {
        coefficients[i] = Y[i];
        for (int j = i + 1; j < n; ++j) {
            coefficients[i] -= V[i][j] * coefficients[j];
        }
        coefficients[i] /= V[i][i];
    }

    return coefficients[0]; // c is the constant term
}

int main() {
    // Example input
    std::map<int, std::pair<int, std::string>> data = {
        {1, {10, "28735619723837"}},
        {2, {16, "1A228867F0CA"}},
        {3, {12, "32811A4AA0B7B"}},
        {4, {11, "917978721331A"}},
        {5, {16, "1A22886782E1"}},
        {6, {10, "28735619654702"}},
        {7, {14, "71AB5070CC4B"}},
        {8, {9, "122662581541670"}},
        {9, {8, "642121030037605"}}
    };

    std::vector<std::pair<int, double>> roots;
    for (const auto &entry : data) {#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <sstream>
#include <iomanip>

// Function to decode a value from a given base to decimal
long long decode_value(int base, const std::string &value) {
    return std::stoll(value, nullptr, base);
}

// Lagrange interpolation to find the polynomial value at x = 0
double lagrange_interpolation(const std::vector<std::pair<int, double>> &points, int x) {
    double total = 0.0;
    int n = points.size();
    for (int i = 0; i < n; ++i) {
        double xi = points[i].first;
        double yi = points[i].second;
        double term = yi;
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                double xj = points[j].first;
                term *= (x - xj) / (xi - xj);
            }
        }
        total += term;
    }
    return total;
}

// Vandermonde method to find the polynomial coefficients
double vandermonde_method(const std::vector<std::pair<int, double>> &points) {
    int n = points.size();
    std::vector<std::vector<double>> V(n, std::vector<double>(n));
    std::vector<double> Y(n);

    // Fill the Vandermonde matrix
    for (int i = 0; i < n; ++i) {
        Y[i] = points[i].second;
        for (int j = 0; j < n; ++j) {
            V[i][j] = pow(points[i].first, j);
        }
    }

    // Solving V * C = Y for coefficients C using Gaussian elimination
    for (int i = 0; i < n; ++i) {
        // Partial pivoting
        for (int j = i + 1; j < n; ++j) {
            if (fabs(V[j][i]) > fabs(V[i][i])) {
                std::swap(V[i], V[j]);
                std::swap(Y[i], Y[j]);
            }
        }
        for (int j = i + 1; j < n; ++j) {
            double ratio = V[j][i] / V[i][i];
            for (int k = i; k < n; ++k) {
                V[j][k] -= ratio * V[i][k];
            }
            Y[j] -= ratio * Y[i];
        }
    }

    // Back substitution
    std::vector<double> coefficients(n);
    for (int i = n - 1; i >= 0; --i) {
        coefficients[i] = Y[i];
        for (int j = i + 1; j < n; ++j) {
            coefficients[i] -= V[i][j] * coefficients[j];
        }
        coefficients[i] /= V[i][i];
    }

    return coefficients[0]; // c is the constant term
}

int main() {
    // Example input
    std::map<int, std::pair<int, std::string>> data = {
        {1, {10, "28735619723837"}},
        {2, {16, "1A228867F0CA"}},
        {3, {12, "32811A4AA0B7B"}},
        {4, {11, "917978721331A"}},
        {5, {16, "1A22886782E1"}},
        {6, {10, "28735619654702"}},
        {7, {14, "71AB5070CC4B"}},
        {8, {9, "122662581541670"}},
        {9, {8, "642121030037605"}}
    };

    std::vector<std::pair<int, double>> roots;
    for (const auto &entry : data) {
        int base = entry.second.first;
        const std::string &value = entry.second.second;
        long long decoded_value = decode_value(base, value);
        roots.emplace_back(entry.first, static_cast<double>(decoded_value));
    }

    // Calculate c using both methods
    double c_lagrange = lagrange_interpolation(roots, 0);
    double c_vandermonde = vandermonde_method(roots);

    // Print results from all methods
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Lagrange Interpolation Result: c = " << c_lagrange << std::endl;
    std::cout << "Vandermonde Method Result: c = " << c_vandermonde << std::endl;

    // Check if results match
    if (fabs(c_lagrange - c_vandermonde) < 1e-9) {
        std::cout << "Both methods yield the same constant term: c = " << c_lagrange << std::endl;
    } else {
        std::cout << "Results do not match!" << std::endl;
    }

    return 0;
}

        int base = entry.second.first;
        const std::string &value = entry.second.second;
        long long decoded_value = decode_value(base, value);
        roots.emplace_back(entry.first, static_cast<double>(decoded_value));
    }

    // Calculate c using both methods
    double c_lagrange = lagrange_interpolation(roots, 0);
    double c_vandermonde = vandermonde_method(roots);

    // Print results from all methods
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Lagrange Interpolation Result: c = " << c_lagrange << std::endl;
    std::cout << "Vandermonde Method Result: c = " << c_vandermonde << std::endl;

    // Check if results match
    if (fabs(c_lagrange - c_vandermonde) < 1e-9) {
        std::cout << "Both methods yield the same constant term: c = " << c_lagrange << std::endl;
    } else {
        std::cout << "!" << std::endl;
    }

    return 0;
}

