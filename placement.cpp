/*
 * placement.cpp
 *
 *  Created on: 24 Mar, 2014
 *      Author: chinhau5
 */

#include "pch.hpp"
#include "placement.hpp"

void Placement::load(const string &filename)
{

	FILE *file;
	char buffer[256];
	char block_name[256];
	int block_number;
	int x, y, z;

	file = fopen(filename.c_str(), "r");

	fgets(buffer, sizeof(buffer), file);

	fscanf(file, "Array size : %d x %d logic blocks", &nx, &ny);
	nx += 2;
	ny += 2;

	fgets(buffer, sizeof(buffer), file);
	fgets(buffer, sizeof(buffer), file);
	fgets(buffer, sizeof(buffer), file);
	fgets(buffer, sizeof(buffer), file);

	while (fscanf(file, "%s %d %d %d #%d ", block_name, &x, &y, &z, &block_number) == 5) {
		assert(block_positions.find(block_name) == block_positions.end());
/*		cout << block_name << " " << x << " " << y << endl;*/
		block_positions.insert(pair<string, pair<int, int>>(block_name, pair<int, int>(x, y)));
	}

	fclose(file);
}

pair<int, int> Placement::getBlockPosition(const string &block_name) const
{
	pair<int, int> position;

	PositionContainer::const_iterator iter = block_positions.find(block_name);

	if (iter != block_positions.end()) {
		position = iter->second;
	} else {
		position.first = -1;
		position.second = -1;
	}

	return position;
}


