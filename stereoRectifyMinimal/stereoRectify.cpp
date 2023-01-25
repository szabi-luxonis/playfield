#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigen>

namespace {
    Eigen::Vector3d rotationMatrixToEulerAngles(Eigen::Matrix3d R) {
    
        double sy = std::sqrt(R(0,0) * R(0,0) +  R(1,0) * R(1,0));
        bool singular = sy < 1e-6;
        double x,y,z;
        if (!singular) {
            x = std::atan2(R(2,1), R(2,2));
            y = std::atan2(-R(2,0), sy);
            z = std::atan2(R(1,0), R(0,0));
        } else {
            x = std::atan2(-R(1,2), R(1,1));
            y = std::atan2(-R(2,0), sy);
            z = 0;
        }
    
        return Eigen::Vector3d(x, y, z);
    }

    Eigen::Matrix3d eulerAnglesToRotationMatrix(Eigen::Vector3d rotationVector) {
    
        Eigen::Matrix3d R_x, R_y, R_z;
        R_x << 1, 0, 0, 
                0, std::cos(rotationVector(0)), -std::sin(rotationVector(0)),
                0, std::sin(rotationVector(0)), std::cos(rotationVector(0));
        
        R_y << std::cos(rotationVector(1)), 0, std::sin(rotationVector(1)),
                0, 1, 0,
                -std::sin(rotationVector(1)), 0, std::cos(rotationVector(1));
    
        R_z << std::cos(rotationVector(2)), -std::sin(rotationVector(2)), 0,
                std::sin(rotationVector(2)), std::cos(rotationVector(2)), 0,
                0, 0, 1;
    
        Eigen::Matrix3d R = R_z * (R_y * R_x);
    
        return R;
    }

    void stereoRectify(Eigen::Matrix3d R, Eigen::Vector3d T, Eigen::Matrix3d& R1, Eigen::Matrix3d& R2) {

        auto om = rotationMatrixToEulerAngles(R);
        om = om * -0.5;
        auto r_r = eulerAnglesToRotationMatrix(om);
        auto t = r_r * T;

        int idx = std::abs(t(0)) > std::abs(t(1)) ? 0 : 1;

        auto c = t(idx);
        auto nt = t.norm();
        Eigen::Vector3d uu = Eigen::Vector3d::Zero();
        uu(idx) = c > 0 ? 1 : -1;
        
        auto ww = t.cross(uu);
        auto nw = ww.norm();

        if (nw > 0) {
            auto scale = std::acos(std::abs(c)/nt)/nw;
            ww = ww * scale;
        }
            
        auto wR = eulerAnglesToRotationMatrix(ww);
        R1 = wR * r_r.transpose();
        R2 = wR * r_r;
    }
}

int main() {

    // reference values
    Eigen::Matrix3d refR1;
    refR1 <<  0.99994761, -0.00390605, -0.0094629,
            0.00391631,  0.99999177,  0.00106613,
            0.00945866, -0.00110314,  0.99995464;

    Eigen::Matrix3d refR2;
    refR2 << 0.99984103,  0.0033367 ,  0.01751487,
            -0.0033177 ,  0.99999386, -0.00111375,
            -0.01751848,  0.00105546,  0.99984598;

    // extrinsic data
    Eigen::Matrix3d R;
    R << 0.99960995, -0.00720377, -0.02698262,
        0.00726279,  0.99997145,  0.00208996, 
        0.02696679, -0.00228512,  0.99963373;

    Eigen::Vector3d T;
    T << -7.51322365, -0.02507336, -0.13161406;

    // minimal rectification code
    Eigen::Matrix3d R1, R2;
    stereoRectify(R, T, R1, R2);
    
    std::cout << std::fixed;

    std::cout << "refR1-R1:" << std::endl;
    std::cout << refR1-R1 << std::endl;
    std::cout << "refR2-R2:" << std::endl;
    std::cout << refR2-R2 << std::endl;

    return 0;

}