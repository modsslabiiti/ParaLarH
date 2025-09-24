/*
 * placement.h
 *
 *  Created on: 24 Mar, 2014
 *      Author: chinhau5
 */

#ifndef PLACEMENT_H_
#define PLACEMENT_H_

#include <string>
#include <map>

using namespace std;

class Placement {
	public:
		Placement(const string &filename) { load(filename); }

		pair<int, int> getBlockPosition(const string &block_name) const;
		pair<int, int> getGridSize() const { return pair<int, int>(nx, ny); }

	private:
		void load(const string &filename);

	private:
		typedef map<string, pair<int, int>> PositionContainer;
		PositionContainer block_positions;
		int nx, ny;

};



#endif /* PLACEMENT_H_ */
