#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

typedef struct _Point {
	_Point() {} 
	_Point(int x, int y) : x(x), y(y) {}
	int x;
	int y;
	bool operator<(const _Point &other) const {
/*		if (x != other.x) {*/
/*			return x < other.x;*/
/*		}*/
/*		return y < other.y;*/
		return std::tie(x, y) < std::tie(other.x, other.y);
	}
	bool operator==(const _Point &other) const {
		return x == other.x && y == other.y;
	}
	bool operator!=(const _Point &other) const {
		return !(*this == other);
	}
	friend std::ostream &operator<<(std::ostream &out, const _Point &p) {
		out << p.x << "," << p.y;
		return out;
	}
} Point;

#endif
