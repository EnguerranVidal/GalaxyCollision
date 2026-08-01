#pragma once

#include <cmath>


struct Vector3
{
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;

    constexpr Vector3() = default;

    constexpr Vector3(float x, float y, float z): x(x), y(y), z(z) {}

    constexpr Vector3 operator+(const Vector3& other) const {return {x + other.x, y + other.y, z + other.z};}

    constexpr Vector3 operator-(const Vector3& other) const {return {x - other.x, y - other.y, z - other.z};}

    constexpr Vector3 operator-() const {return {-x, -y, -z};}

    constexpr Vector3 operator*(float scalar) const {return {x * scalar, y * scalar, z * scalar};}

    constexpr Vector3 operator/(float scalar) const {return {x / scalar, y / scalar, z / scalar};}

    Vector3& operator+=(const Vector3& other)
    {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    Vector3& operator-=(const Vector3& other)
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    Vector3& operator*=(float scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    [[nodiscard]] float squaredNorm() const{return x*x + y*y + z*z;}

    [[nodiscard]] float norm() const {return std::sqrt(squaredNorm());}

    void normalize()
    {
        float length = norm();
        if (length > 0.0f)
        {
            x /= length;
            y /= length;
            z /= length;
        }
    }

    [[nodiscard]] Vector3 normalized() const
    {
        Vector3 result = *this;
        result.normalize();
        return result;
    }
};

constexpr Vector3 operator*(float scalar, const Vector3& vector) {return vector * scalar;}


inline float dot(const Vector3& a, const Vector3& b) {return a.x * b.x + a.y * b.y + a.z * b.z;}

inline Vector3 cross(const Vector3& a, const Vector3& b) {return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};}