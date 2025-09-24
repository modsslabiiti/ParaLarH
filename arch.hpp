#ifndef ARCH_H
#define ARCH_H

using namespace std;
using namespace pugi;

class Arch {
	public:
	void load(const string &filename);
	void parse_complex_block_list(const xml_node &cb_node);
	void parse_pb_type(const xml_node &pb_node, struct Pb *pb, int level);
	void parse_delays(const xml_node &interconnect_node, struct Interconnect *interconnect, int level);
	void parse_interconnect(const xml_node &interconnect_node, struct Mode *mode, int level);
	void process_lut_class(struct Pb *lut_pb);

	map<string, Pb *> pbs;
};

enum { DELAY_CONSTANT, DELAY_MATRIX };

typedef struct DelayPort {
	string instance;
	int instance_low;
	int instance_high;
	string port;
	int port_low;
	int port_high;
} DelayPort;

typedef struct Delay {
	int type;
	float **matrix;
	bool matches(const string &in_instance, int in_instance_index, const string &in_port, int in_pin, const string &out_instance, int out_instance_index, const string &output_port, int out_pin) const;
/*	bool matches_input(const string &instance_name, const string &port_name) const;*/
	DelayPort in_port;
	DelayPort out_port;
} Delay;

typedef struct Interconnect {
	string name;
	vector<Delay *> delays;
} Interconnect;

typedef struct Mode {
	string name;
	map<string, struct Pb *> pbs;
	vector<Interconnect *> interconnects;
} Mode;

typedef struct Pb {
	string name;
	string class_name;
	string blif_model; 
	map<string, Mode *> modes;
} Pb; 

void parse_delay_input_output(const string &in_port, const string &out_port, Delay *delay);


#endif

