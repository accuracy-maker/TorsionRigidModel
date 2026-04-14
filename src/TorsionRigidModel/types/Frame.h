#pragma once

#include<array>
#include<cmath>
#include<ostream>
#include<iomanip>

/* Frame represent the world frame W(s) in the paper*/

namespace ctr {

    // 3D vector
    struct Vec3 {
        double x = 0.0, y = 0.0, z = 0.0;
    
        Vec3() = default;
        Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
    
        Vec3 operator+(const Vec3& v) const { return {x+v.x, y+v.y, z+v.z}; }
        Vec3 operator-(const Vec3& v) const { return {x-v.x, y-v.y, z-v.z}; }
        Vec3 operator*(double s) const { return {x*s, y*s, z*s}; }
    
        double norm() const { return std::sqrt(x*x + y*y + z*z); }
        double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
    };

    /**
    * 3x3 rotation matrix, stored row-major.
    *
    *   [m[0][0]  m[0][1]  m[0][2]]
    *   [m[1][0]  m[1][1]  m[1][2]]
    *   [m[2][0]  m[2][1]  m[2][2]]
    */
    
    struct Mat3 {
        double m[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};

        static Mat3 identity() {
            Mat3 I;
            return I;
        }

        // Matrix-vector product: R * v
        Vec3 operator*(const Vec3& v) const {
            return {
                m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z,
                m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z,
                m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z
            };
        }

        // matrix-matrix product: A * B
        Mat3 operator*(const Mat3& B) const {
            Mat3 C;
            for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++) {
                        C.m[i][j] = 0.0;
                        for (int k = 0; k < 3; k++)
                            C.m[i][j] += m[i][k] * B.m[k][j];
                    }
                return C;
        }

        // Transpose
            Mat3 transpose() const {
            Mat3 T;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    T.m[i][j] = m[j][i];
            return T;
        }

        // Rotate about z-axis angle theta
        static Mat3 rotZ(double theta) {
            double c = std::cos(theta), s = std::sin(theta);
            Mat3 R;
            R.m[0][0] =  c; R.m[0][1] = -s; R.m[0][2] = 0;
            R.m[1][0] =  s; R.m[1][1] =  c; R.m[1][2] = 0;
            R.m[2][0] =  0; R.m[2][1] =  0; R.m[2][2] = 1;
            return R;
        }
    };

    /**
    * Rigid body frame in SE(3): rotation R + position p.
    *
    * Represents the 4x4 homogeneous transform:
    *   [R  p]
    *   [0  1]
    *
    * In the paper:
    *   W(s) is the needle frame at arc length s
    *   R = orientation of the cross-section
    *   p = position of the cross-section center
    *   z-axis of R is tangent to the backbone
    */
    struct Frame {
        Mat3 R; // orientation
        Vec3 p; // position

        Frame() : R(Mat3::identity()), p(0, 0, 0) {}

        Frame(const Mat3& R, const Vec3& p) : R(R), p(p) {}

        static Frame identity() {
            return Frame();
        }

        // W(s) = T1 * T2 * T3 * Tm, where m is the number of sections
        // T1 * T2 = [R1  p1] * [R2  p2] = [R1R2  R1p2 + p1]
        //           [0    1]   [0    1]   [0          1   ]
        Frame operator*(const Frame& other) const {
            return Frame(
                R * other.R,
                p + R * other.p
            );
        }

        // Get the z-axis (tangent direction) of this frame
        Vec3 tangent() const {
            return {R.m[0][2], R.m[1][2], R.m[2][2]};
        }
    
        // Get the x-axis of this frame
        Vec3 xAxis() const {
            return {R.m[0][0], R.m[1][0], R.m[2][0]};
        }
    
        // Get the y-axis of this frame
        Vec3 yAxis() const {
            return {R.m[0][1], R.m[1][1], R.m[2][1]};
        }

        // Print frame
        friend std::ostream& operator<<(std::ostream& os, const Frame& f) {
            os << std::fixed << std::setprecision(6);
            os << "p=(" << f.p.x << ", " << f.p.y << ", " << f.p.z << ")";
            return os;
        }
    };


} // namespace ctr
