#include "pch.hpp"
#include "netlist.hpp"
using namespace std;
using namespace boost;

/*#define VERBOSE*/

/*typedef Graph<adjacency_list<vec*/
 

vector<Net> read_netlist(const string &filename)
{
	fstream file(filename, ios_base::in);
	string line;
	vector<Net> nets;
	int num_nets = 0;
	
	int line_num = 0;
	while (!file.eof()) {
		char_separator<char> sep(" ");
		getline(file, line);
		tokenizer<char_separator<char>> tokens(line, sep);

		if (line.empty()) {
			continue;
		}
#ifdef VERBOSE
		printf("line: %s\n", line.c_str());
#endif
		if (line_num == 0) {
			num_nets = stoi(line);
			nets.reserve(num_nets);
#ifdef VERBOSE
			printf("num_nets: %d\n", num_nets);
#endif
		} else {
			int i = 0;
			Net net;
			net.index = nets.size();
			for (const auto &token : tokens) {
#ifdef VERBOSE
				printf("token: %s\n", token.c_str());
#endif
				switch (i) {
					case 0: /* net name */
						net.name = token;
						break;
					case 1: /* num points */
						{
							int num_points = stoi(token);
#ifdef VERBOSE
							printf("num_points: %d\n", num_points);
#endif
							net.points.reserve(num_points);
							break;
						}
					default: /* points */
						{
							tokenizer<char_separator<char>> point_tokens(token, char_separator<char>(","));
							tokenizer<char_separator<char>>::const_iterator iter = point_tokens.begin();
							int x = stoi(*iter++);
							int y = stoi(*iter++);
							assert(iter == point_tokens.end());
#ifdef VERBOSE
							printf("x=%d y=%d\n", x, y);
#endif
							net.points.emplace_back(Point(x, y));
							break;
						}
						break;
				}
				++i;
			}
			nets.emplace_back(net);
		}
		++line_num;
	}
#ifdef VERBOSE
	printf("nets.size %d num_nets %d\n", nets.size(), num_nets);
#endif
	assert(nets.size() == num_nets);
	return nets;
}
