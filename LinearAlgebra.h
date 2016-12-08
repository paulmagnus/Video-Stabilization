#ifndef LINEAR_ALGEBRA
#define LINEAR_ALGEBRA
#include "Output.h"
#include <math.h>

typedef struct vec2 {
	int x;
	int y;

	vec2(int x = 0, int y = 0) : x(x), y(y) {}

	bool operator==(const vec2& v) const {
		return x == v.x && y == v.y;
	}

	vec2 operator+(const vec2& v) const {
		return vec2(x + v.x, y + v.y);
	}

	vec2 operator+=(const vec2& v) {
		x += v.x;
		y += v.y;
		return *this;
	}

	vec2 operator-(const vec2& v) const {
		return vec2(x - v.x, y - v.y);
	}

	vec2 operator-=(const vec2& v) {
		x -= v.x;
		y -= v.y;
		return *this;
	}

	/* Scalar multiplication of the vector */
	vec2 operator*(const int s) const {
		return vec2(x * s, y * s);
	}

	/* Scalar multiplication of the vector */
	vec2 operator*=(const int s) {
		x *= s;
		y *= s;
		return *this;
	}

	/* Elementwise multiplication of the vector */
	vec2 operator*(const vec2& v) const {
		return vec2(x * v.x, y * v.y);
	}

	/* Elementwise multiplication of the vector */
	vec2 operator*=(const vec2& v) {
		x *= v.x;
		y *= v.y;
		return *this;
	}

	/* Prints the vector for debugging purposes */
	/*
	void print() {
		::print("[%d, %d]", x, y);
	}
	*/

	double length() const {
		return sqrt(pow((float)x, 2) + pow((float)y, 2));
	}

	/* The dot product of two vectors*/
	float dot(const vec2& v) const {
		return x * v.x + y * v.y;
	}

	vec2 normalize() {
		return *this * (1 / length());
	}
};

typedef struct vec3 {
	double x, y, z;

	vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

	/* Addition of two matricies */
	vec3 operator+(const vec3& v) const {
		return vec3(x + v.x, y + v.y, z + v.z);
	}

	/* Addition of two matricies */
	vec3 operator+=(const vec3& v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	/* Scalar multiplication of the vector */
	vec3 operator*(const double s) const {
		return vec3(x * s, y * s, z * s);
	}

	/* Scalar multiplication of the vector */
	vec3 operator*=(const double s) {
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}

	/* Elementwise multiplication of the vector */
	vec3 operator*(const vec3& v) const {
		return vec3(x * v.x, y * v.y, z * v.z);
	}

	/* Elementwise multiplication of the vector */
	vec3 operator*=(const vec3& v) {
		x *= v.x;
		y *= v.y;
		z *= v.z;
		return *this;
	}

	/* Prints the vector for debugging purposes */
	/*
	void print() {
		::print("[%f, %f, %f]\n", x, y, z);
	}
	*/
};

struct mat3 {
	double m[3][3];

	mat3(double m00 = 0.0, double m01 = 0.0, double m02 = 0.0,
		double m10 = 0.0, double m11 = 0.0, double m12 = 0.0,
		double m20 = 0.0, double m21 = 0.0, double m22 = 0.0) {
		m[0][0] = m00; m[0][1] = m01; m[0][2] = m02;
		m[1][0] = m10; m[1][1] = m11; m[1][2] = m12;
		m[2][0] = m20; m[2][1] = m21; m[2][2] = m22;
	}

	vec3 operator*(const vec3& v) const {
		return vec3(m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
			m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
			m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z);
	}

	mat3 translation(double dx, double dy) {
		return mat3(
			1, 0, dx,
			0, 1, dy,
			0, 0, 1);
	}


	void print() const {
		printf("begin matrix\n");
		printf("%f %f %f\n", m[0][0], m[0][1], m[0][2]);
		printf("%f %f %f\n", m[1][0], m[1][1], m[1][2]);
		printf("%f %f %f\n", m[2][0], m[2][1], m[2][2]);
		printf("end matrix\n\n");
	}
};
#endif
