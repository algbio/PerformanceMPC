#include <fstream>
#include <algorithm>
#include <string>
#include "cassert"
#include <mpc/utils.h>

long long mem_peak() {
	std::ifstream fs("/proc/self/status", std::ios::in);
	if(!fs.good()) {
		assert(false);
		return 0;
	}
	if(fs.is_open()) {
		for(std::string s; std::getline(fs, s);) {
			if(prefix(s, "VmPeak:"))
				return std::stoll(s.substr(7, s.size()-2-7));
		}
	}
	assert(false);
	return 0;
}

bool prefix(std::string &s, std::string prefix) {
	for(int i=0; i<prefix.size(); i++)
		if(s[i] != prefix[i])
			return false;
	return true;
}

void log_time(stopwatch::time_used t, nlohmann::json &j) {
	j["sys"] = t.sys;
	j["usr"] = t.usr;
	j["real"] = t.real;
}
