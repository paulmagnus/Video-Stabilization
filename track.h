#ifndef TRACK
#define TRACK
#include <vector>
#include "LinearAlgebra.h"

struct track {
	std::vector<vec2> positions;

	track() : positions(0) {}

	void addPosition(const vec2 & pos) {
		positions.push_back(pos);
	}

	vec2 getPosition(int frame) {
		return positions[frame];
	}

	int size() {
		return positions.size();
	}

	// This returns true if the given feature is found in the given frame or false if it is not
	bool isValid(int frame) {
		return positions[frame].x > 0 && positions[frame].y > 0;
	}
	/*
	void printOut() {
		print("[");
		for (int i = 0; i < positions.size(); i++) {
			print("(%f, %f)", positions[i].x, positions[i].y);
			if (i != positions.size() - 2)
				print(", ");
		}
		print("]\n");
	}
	*/
};
#endif