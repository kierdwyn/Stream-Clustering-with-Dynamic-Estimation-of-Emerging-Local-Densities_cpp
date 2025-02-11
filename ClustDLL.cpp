// ClustDLL.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"

#include <iostream>
#include <fstream>
//#include <climits.h>

#include "clustr.h"
using namespace std;

//void DTT_clustering (Dictionary &Dic, DTT_parameters &Parameters, vector <double *> &D ){
//    Temperate_Table Table(Parameters);
//    unsigned seed = (unsigned) chrono::system_clock::now().time_since_epoch().count();
//    shuffle (D.begin(), D.end(), default_random_engine(seed));
//    for (int i=0; i<D.size(); i++) {
//        process_one_point(Dic, Table, D[i]);
//    }
//}

void WangClust::process_one_point(double *p) {
	temp_kernels.decrease();
	int k = kernels.match(p);
	if (k >= 0) {
		kernels.update(k, p);
	}
	else {
		if (kernels.centers.size() < max_numc) {
			k = temp_kernels.match(p);
			temp_kernels.update(k, p, kernels);
		}
	}
}

void WangClust::feed_data(vector <double *> &D) {
	for (int i = 0; i<D.size(); i++) {
		process_one_point(D[i]);
	}
}

void WangClust::copy_centers(vector<double *> &c) {
	for (int i = 0; i < kernels.centers.size(); i++) {
		double * p = new double[parameters.dim];
			for (int j = 0; j < parameters.dim; j++)
				p[j] = kernels.centers[i].value[j];
		c.push_back(p);
	}
}


double similarity_measure(double *x, double *y, string opt, int dim) {
	if (opt.compare("euclidean") == 0) {
		double s = 0;
		for (int i = 0; i<dim; i++) {
			s += (x[i] - y[i])*(x[i] - y[i]);
		}
		s = sqrt(s);
		return -s;
	}
	else if (opt.compare("cosine") == 0) {
		double s = cosine_similarity(x, y, dim);
		return s;
	}
	else {
		cout << "no such option" << endl;
		return -1;
	}
}

double cosine_similarity(double *a, double *b, int dim) {

	double sum1 = 0;
	double sum_a = 0;
	double sum_b = 0;

	for (int i = 0; i<dim; i++) {
		sum1 += a[i] * b[i];
		sum_a += a[i] * a[i];
		sum_b += b[i] * b[i];
	}
	return (double)sum1 / (sqrt((double)sum_a)*sqrt((double)sum_b));
}






//-------------------------------------------Dictionary methods----------------------------------------------------
void Dictionary::SetDictionary(DTT_parameters & Parameters) {
	dim = Parameters.dim;
	sim_thrsh = Parameters.sim_thrsh;
}


int Dictionary::match(double * p) {


	if (centers.size() == 0) {
		return INT_MIN;
	}

	double max_sim_score = -DBL_MAX;
	int max_sim_indx = INT_MIN;

	// find the best match in the Dictionary
	for (int k = 0; k<centers.size(); k++) {
		double similarity_score = similarity_measure(p, centers[k].value, "cosine", dim);
		if (similarity_score>max_sim_score) {
			max_sim_score = similarity_score;
			max_sim_indx = k;
		}
	}
	//cout << max_sim_score << endl;
	if (max_sim_score>sim_thrsh) {
		return max_sim_indx;
	}
	else if (max_sim_indx<0)
		return max_sim_indx;
	else
		return -max_sim_indx - 1;

}
void Dictionary::update(int k, double * p) {
	for (int l = 0; l<dim; l++)
		centers[k].value[l] = (centers[k].value[l] * centers[k].num_member + p[l]) / (centers[k].num_member + 1);
	centers[k].num_member++;
}

void Dictionary::Dic_write(string filename) {
	ofstream c;
	//    centers.open("/Users/WangZiyin/Documents/DTT_centers.txt");
	c.open(filename);
	for (int i = 0; i<centers.size(); i++) {
		for (int j = 0; j<dim; j++) {
			c << centers[i].value[j] << ' ';
		}
		c << endl;
	}
	c.close();
}


//-------------------------------------Temperate Table methods----------------------------------------------------
void Temperate_Table::SetTemperate_Table(DTT_parameters Parameters) {
	dim = Parameters.dim;
	sim_thrsh = Parameters.sim_thrsh;
	int_val = Parameters.init_val;
	upper_thresh = Parameters.upper_thrsh;
}


int Temperate_Table::match(double * p) {
	double max_sim_score = -DBL_MAX;
	int max_sim_indx = INT_MIN;

	for (int k = 0; k<centers.size(); k++) {
		double similarity_score = similarity_measure(p, centers[k].value, "cosine", dim);
		//cout << similarity_score << endl;
		if (similarity_score>max_sim_score) {
			max_sim_score = similarity_score;
			max_sim_indx = k;
		}
	}
	//cout << max_sim_score << endl;
	if (max_sim_score>sim_thrsh)
		return max_sim_indx;
	else return -1;
}

void Temperate_Table::update(int k, double * p, Dictionary &Dic) {
	if (k >= 0) {
		for (int l = 0; l<dim; l++) {
			centers[k].value[l] = (centers[k].value[l] * centers[k].num_member + p[l]) / (centers[k].num_member + 1);
		}
		centers[k].num_member++;
		centers[k].active_value += int_val;// inrease active value by the initial value N_0

		if (centers[k].active_value>upper_thresh) {// if the active value increased higher than threshold

			Dic_member new_dic;
			new_dic.num_member = centers[k].num_member;
			new_dic.value = centers[k].value; // assign
			Dic.centers.push_back(new_dic);
			centers.erase(centers.begin() + k);
		}

	}
	else {
		Table_member a;
		a.num_member = 1;
		a.value = new double[dim];
		for (int i = 0; i < dim; i++)
			a.value[i] = p[i];
		a.active_value = int_val;
		centers.push_back(a);
	}
}

void Temperate_Table::decrease() {
	if (centers.size()>0) {

		for (int i = 0; i<centers.size(); i++) {
			centers[i].active_value -= 1;

			if (centers[i].active_value<0) {
				delete[] centers[i].value;
				centers.erase(centers.begin() + i);
				i--;
			}
		}
	}
}


void Dictionary::Assign(vector<int> &A, vector<double *> &D) {
	for (int i = 0; i<D.size(); i++) {
		int k = match(D[i]);
		if (k<0) {
			k = abs(k) - 1;
		}
		A.push_back(k);
	}
}


double * cpvector2array2D(vector<double *> x, int c) {
	int r = x.size();
	double * D = new double[r*c];
	for (int i = 0; i<r; i++) {
		double *p = new double[c];
		for (int j = 0; j < c; j++)
			D[i*c+j] = x[i][j];
	}
	return D;
}

void cparray2vector2D(vector<double*> &x, double * D, int r, int c) {
	for (int i = 0; i < r; i++) {	 
		double * p = new double[c];
		for (int j = 0; j < c; j++)
			p[j] = D[i*c+j];
		x.push_back(p);
	}
}


void* clust_new(double r, double th, int dim, int numc) {
	WangClust* cc = new WangClust(r, th, dim, numc);
	cout << dim << endl;
	cout << cc->parameters.dim << endl;
	cout << cc->parameters.pf << endl;
	cout << cc->parameters.sim_thrsh << endl;
	cout << cc->parameters.init_val << endl;
	cout << cc->parameters.upper_thrsh << endl;
	cout << cc->max_numc << endl;
	return (void *)cc;
}

void clust_feed(void * indata, int r, int c, void * handle) {
	double * D = (double *)indata;
	WangClust * clust = (WangClust *)handle;
	vector<double *> x;
	
	cparray2vector2D(x, D, r, c);
	cout << "data fed, num of instances "<< x.size() << endl;

	cout << "data converted" << endl;
	cout << clust->parameters.dim << endl;
	clust->feed_data(x);
	//cout << clust->kernels.centers.size() << endl;
	for (int i = 0; i < x.size(); i++)
		delete[] x[i];
}

void clust_get(void *out_data, void *handle) {
	WangClust * clust = (WangClust *)handle;
	double * C = (double *)out_data;
	vector<double *> centers;
	clust->copy_centers(centers);
	int dim = clust->parameters.dim;
	for (int i = 0; i < centers.size(); i++)
		for (int j = 0; j < dim; j++)
			C[i*dim + j] = centers[i][j];

	for (int i = 0; i < centers.size(); i++)
		delete[] centers[i];
}

int num_cluster(void *handle) {
	WangClust * clust = (WangClust *)handle;
	//cout << (int)clust->kernels.centers.size() << endl;
	return (int)clust->kernels.centers.size();
}

















