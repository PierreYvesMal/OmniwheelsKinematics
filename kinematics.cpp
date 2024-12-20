#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>

const double PI = 3.141592;
const double DEG_TO_RAD = PI/180;
const double r = 30;
const double R = 116;

const double θ = 60*DEG_TO_RAD;
const double a[3] = {0, 120*DEG_TO_RAD, 240*DEG_TO_RAD};

//https://www.internationaljournalssrg.org/IJEEE/2019/Volume6-Issue12/IJEEE-V6I12P101.pdf
std::pair<std::vector<double>, std::vector<double>> kinematics_matrix_time(double wanted_distance, double wanted_angle, double wanted_rotation, double wanted_speed);

int main()
{
    double wanted_distance = 100;
    double wanted_angle = PI/2;
    double wanted_rotation = 0;
    double wanted_speed = 100;

    const auto [W, ticks] = kinematics_matrix_time(wanted_distance, wanted_angle, wanted_rotation, wanted_speed);
    std::vector<double> V;
    for(int i=0; i<3; ++i)
    {
        V.push_back(W[i]*r);
        std::cout << "W" << i << ": " << W[i] << "\t=>V" << i << ": " << V[i] << std::endl;
    }

    std::cout << "Verif:" << std::endl;

    double vx = V[0]*(-sin(θ)) + V[2]*(-sin(θ+a[2]));
    double vy = V[0]*cos(θ) - V[1] + V[2]*cos(θ+a[2]);
    double v = sqrt(vx*vx+vy*vy);
    std::cout << "vx: " << vx << std::endl;
    std::cout << "vy: " << vy << std::endl;
    std::cout << "v: " << v << std::endl;

    double rotation = R*(V[0]+V[1]+V[2]);
    std::cout << "Rotation: " << rotation << std::endl;
}

std::vector<double> dotProduct(const std::vector<std::vector<double>>& K, const std::vector<double>& U) {
    // Ensure the dimensions are correct
    if (K.size() != 3 || U.size() != 3) {
        throw std::invalid_argument("Matrix dimensions must be 3x3 for K and 1x3 for U");
    }
    for (const auto& row : K) {
        if (row.size() != 3) {
            throw std::invalid_argument("Matrix K must be 3x3");
        }
    }

    // Result vector for the 1x3 output
    std::vector<double> result(3, 0.0);

    // Perform the dot product calculation
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            result[i] += K[i][j] * U[j];
        }
    }

    return result;
}

std::pair<std::vector<double>, std::vector<double>>  kinematics_matrix_time(double wanted_distance, double wanted_angle, double wanted_rotation, double wanted_speed)
{
    // The global robot mouvements
    std::vector<double> U = {wanted_speed * cos(wanted_angle), wanted_speed * sin(wanted_angle), 0};//wanted_speed/R};

    // Computed using wolfram alpha
    std::vector<std::vector<double>> invK =
    {{-1/sqrt(3), (double)1/3, R/3}, {0, -(double)2/3, R/3}, {1/sqrt(3), (double)1/3, R/3}};

    // The wheel linear speeds V = invK * U
    std::vector<double> V = dotProduct(invK,U);
    std::cout << "V (mm/s):  " << std::endl;
    for(int i=0; i<3; ++i)
    {
        std::cout << "\t" << i << ": " << V[i] << std::endl;
    }

    // Wheel rotational speed W = V/r
    std::vector<double> W = V;
    for(int i=0; i<3; i++)
    {
        W[i]/=r;
    }

    // The number of ticks to rotate for each wheel
    std::vector<double> ticks;

    std::cout << "Ticks:  " << std::endl;
    for(int i=0; i<3; ++i)
    {
        // The wheel only participate in its plane, and simply drifts on the perpendicular direction
        // The distance that the wheel will travel is the combination of its rotation and drift

        // Portion of D in the rotation direction
        ticks.push_back( std::abs(wanted_distance * sin (θ + a[i] - wanted_angle) * (double) 200 / (2 * PI * r)));
        std::cout << "\t" << i << ": " << ticks[i] << std::endl;
    }

    return {W,ticks};   

}