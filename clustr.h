#ifdef CLUST_EXPORTS
#define CLUST_EXPORTS_API __declspec(dllexport)
#else
#define CLUST_EXPORTS_API __declspec(dllimport)
#endif



#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <random>

using namespace std;


class DTT_parameters {
public:
	double sim_thrsh;
	double belta = 0.01;
	double Fy = 2.23;
	double pf;

	int dim;

	int init_val;
	double upper_thrsh;

	void compute_parameters() {
		init_val = (int)(log(belta) / log(1 - pf) + 1);
		upper_thrsh = (init_val*pf + Fy * sqrt(init_val*pf*(1 - pf)))*init_val;
	}

};




class Dic_member {
public:
	double * value;
	int num_member;

};

class Table_member {
public:
	double * value;
	int num_member;
	int active_value;

};


class Dictionary {
public:
	vector<Dic_member> centers;

	int dim;
	double sim_thrsh;

	int match(double * value);  //returns the index of matched entry (if success), returns a nagative integer if there is no match. The index starts from 0.
	void update(int indx, double * p);
	void Dic_write(string filename);
	void Assign(vector<int> &A, vector<double *> &D);
	void SetDictionary(DTT_parameters &Parameters);

};


class Temperate_Table {
public:
	vector<Table_member> centers;
	double sim_thrsh;
	int int_val;
	double upper_thresh;
	int dim;
	int match(double * value);
	void update(int indx, double * p, Dictionary &Dic); // indx is the return value of match(), it indicates which entry has been matched, -1 if no such mathc
	void decrease();
	void SetTemperate_Table(DTT_parameters Parameters);

};



class WangClust {
public:
	DTT_parameters parameters;
	Dictionary kernels;
	Temperate_Table temp_kernels;
	int max_numc = INT_MAX;

	WangClust(double r, double th, int dim, int numc) {
		parameters.sim_thrsh = th;
		parameters.pf = r;
		parameters.dim = dim;
		parameters.compute_parameters();
		kernels.SetDictionary(parameters);
		temp_kernels.SetTemperate_Table(parameters);
		max_numc = numc;
	}

	void feed_data(vector <double *> &D);
	void copy_centers(vector<double *> &c);
	void process_one_point(double *p);
};

double cosine_similarity(double *a, double *b, int dim);
double similarity_measure(double *x, double *y, string opt, int dim);

extern "C" CLUST_EXPORTS_API  void* clust_new(double r, double th, int dim, int numc);
extern "C" CLUST_EXPORTS_API  void  clust_feed(void * indata, int r, int c, void * handle);
extern "C" CLUST_EXPORTS_API  void  clust_get(void *out_data, void *handle);
extern "C" CLUST_EXPORTS_API  int   num_cluster(void *handle);