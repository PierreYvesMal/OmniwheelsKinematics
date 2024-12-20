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
std::pair<std::vector<double>, std::vector<double>> kinematics(double wanted_distance, double wanted_angle, double wanted_rotation, double wanted_speed, double wanted_rotational_speed);

int main()
{
    /// INSTRUCTIONS
    double wanted_distance = 100;
    double wanted_angle = PI/3;
    double wanted_rotation = 2*PI;
    double wanted_speed = 100; //m/s
    double wanted_rotational_speed = 10; // rad/s 
    std::cout << "Instructions:" << std::endl;
    std::cout << "wanted_distance: " << wanted_distance << " (mm)" << std::endl;
    std::cout << "wanted_angle: " << wanted_angle << " (rad)" << std::endl;
    std::cout << "wanted_rotation: " << wanted_rotation << " (rad)" << std::endl;
    std::cout << "wanted_speed: " << wanted_speed << " (mm/s)" << std::endl;
    std::cout << "wanted_rotational_speed: " << wanted_rotational_speed << " (rad/s)" << std::endl;
    std::cout << std::endl;

    /// COMPUTATIONS
    const auto [W, ticks] = kinematics(wanted_distance, wanted_angle, wanted_rotation, wanted_speed, wanted_rotational_speed);
    std::cout << "Computed values:" << std::endl;
    std::vector<double> V;
    for(int i=0; i<3; ++i)
    {
        V.push_back(W[i]*r);
        std::cout << "W" << i << ": " << (int)W[i] << " (rad/s) \t =>V" << i << ": " << (int)V[i] << " (mm/s) \t for " << (int)ticks[i] << " ticks." << std::endl;
    }
    std::cout << std::endl;

    /// VERIFICATIONS
    std::cout << "Verification:" << std::endl;

    double vx = V[0]*(-sin(θ)) + V[2]*(-sin(θ+a[2]));
    double vy = V[0]*cos(θ) - V[1] + V[2]*cos(θ+a[2]);
    double v = sqrt(vx*vx+vy*vy);
    // std::cout << "vx: " << vx << std::endl;
    // std::cout << "vy: " << vy << std::endl;

    std::cout << "Robot translation speed: " << v << " (mm/s)" << std::endl;

    double rotation = (V[0]+V[1]+V[2])/R;
    std::cout << "Robot rotational speed: " << rotation  << " (rad/s)" << std::endl;
    std::cout << "--------------\n" << std::endl;
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

std::pair<std::vector<double>, std::vector<double>>  kinematics(double wanted_distance, double wanted_angle, double wanted_rotation, double wanted_speed, double wanted_rotational_speed)
{
    // The global robot mouvements
    std::vector<double> U = {wanted_speed * cos(wanted_angle), wanted_speed * sin(wanted_angle), wanted_rotational_speed};

    // Computed using wolfram alpha
    std::vector<std::vector<double>> invK =
    {{-1/sqrt(3), (double)1/3, R/3}, {0, -(double)2/3, R/3}, {1/sqrt(3), (double)1/3, R/3}};

    // The wheel linear speeds V = invK * U
    std::vector<double> V = dotProduct(invK,U);

    // Wheel rotational speed W = V/r
    std::vector<double> W = V;
    for(int i=0; i<3; i++)
    {
        W[i]/=r;
    }

    // The number of ticks to rotate for each wheel
    std::vector<double> ticks(3,0.0);

    for(int i=0; i<3; ++i)
    {
        // The wheel only participate in its plane, and simply drifts on the perpendicular direction
        // The distance that the wheel will travel is the combination of its rotation and drift

        // Portion of D in the rotation direction
        // Portion due to translation
        ticks[i] = std::abs(wanted_distance * sin (θ + a[i] - wanted_angle) * (double) 200 / (2 * PI * r));
        // Portion due to rotation
        ticks[i] += wanted_rotation * R / (r) * 200 / (2*PI);
    }

    return {W,ticks};   

}