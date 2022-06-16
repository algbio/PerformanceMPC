#include <string>
#include <nlohmann/json.hpp>
#include <chrono>
#include <sys/resource.h>
#include <iostream>
#include <nlohmann/json.hpp>

bool prefix(std::string &s, std::string prefix);
long long mem_peak();

enum Type {sw_self, sw_child};
struct stopwatch {
	std::chrono::time_point<std::chrono::high_resolution_clock> t_real, t_lap_real;
	long long t_usr, t_lap_usr, t_sys, t_lap_sys;
	int who;

	stopwatch(Type t=sw_self) {
		if(t == sw_self) {
			who = RUSAGE_SELF;
		} else if(t == sw_child) {
			who = RUSAGE_CHILDREN;
		} else {
			assert(false);
		}
		reset();
	}

	long long tvals(timeval tv) {
		return tv.tv_sec*1000000+tv.tv_usec;
	}

	void reset() {
		t_real = std::chrono::high_resolution_clock::now();
		t_lap_real = t_real;
		rusage u;
		int ret = getrusage(who, &u);
		assert(ret == 0);
		t_usr = tvals(u.ru_utime);
		t_lap_usr = tvals(u.ru_utime);
		t_sys = tvals(u.ru_stime);
		t_lap_sys = tvals(u.ru_stime);
	}

	struct time_used {
		long long real, usr, sys; // In microseconds
	};

	time_used lap() {
		auto real_t = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t_lap_real).count();
		rusage u;
		int ret = getrusage(who, &u);
		assert(ret == 0);

		int usr_t = tvals(u.ru_utime)-t_lap_usr;
		int sys_t = tvals(u.ru_stime)-t_lap_sys;

		t_lap_real = std::chrono::high_resolution_clock::now();
		t_lap_usr = tvals(u.ru_utime);
		t_lap_sys = tvals(u.ru_stime);

		return {real_t, usr_t, sys_t};
	}

	time_used total() {
		auto real_t = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t_real).count();
		rusage u;
		int ret = getrusage(who, &u);
		assert(ret == 0);

		int usr_t = tvals(u.ru_utime)-t_usr;
		int sys_t = tvals(u.ru_stime)-t_sys;

		return {real_t, usr_t, sys_t};
	}
};

void log_time(stopwatch::time_used t, nlohmann::json &j);
