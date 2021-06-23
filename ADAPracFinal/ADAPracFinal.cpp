#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<numeric>
#include <queue>

using namespace std;

bool leerArgumentos(int argc, char* argv[], string& fichero) {

	bool f = false;
	string s = "";
	for (int i = 1; i < argc; i++) {
		s = argv[i];
		if (s == "-f") {
			if (i < argc) {
				f = true;
				i++;
				fichero = argv[i];
			}
			else {
				cout << "Debe introducir nombre de fichero despues de -f" << endl;
				return false;
			}
		}
		else {
			cout << "Argumento desconocido: " << argv[i] << endl;
			cout << "Los argumentos que se aceptan son: potter -f ficheroDeEntrada" << endl;
			return false;
		}
	}
	if (!f) {
		cout << "-f es argumento obligatorio" << endl;
		return false;
	}
	return true;
}

bool leerFichero(string nombre, int& n, double& T, vector<double>& t, vector<double>& v, vector<unsigned>& m) {

	bool abierto = false;
	double aux = 0;
	ifstream f(nombre);
	if (f.is_open()) {
		abierto = true;
		f >> n;
		f >> T;
		//cout << n << " " << T << endl;
		for (int i = 0; i < n; i++) {
			f >> aux;
			t.push_back(aux);
			//cout<<t[i]<<endl;
		}
		for (int i = 0; i < n; i++) {
			f >> aux;
			v.push_back(aux);
			//cout<<v[i]<<endl;
		}
		for (int i = 0; i < n; i++) {
			f >> aux;
			m.push_back(aux);
			//cout<<m[i]<<endl;
		}
	}
	else {
		cout << "El fichero no existe" << endl;
		abierto = false;
		f.close();
	}
	return abierto;
}

ostream& operator<<(ostream& os, const vector<unsigned>& vect) {
	for (unsigned i = 0; i < vect.size(); i++) {
		os << vect[i] << " ";
	}
	return os;
}

double knapsack_c(
	const vector<double>& v, // values
	const vector<double>& w, // weights
	const vector<unsigned>& m,
	unsigned k,
	double W // knapsack weight limit
) {
	/*vector<size_t> idx(w.size()); // objects sorted by value density
	for (size_t i = 0; i < idx.size(); i++) idx[i] = i;

	sort(begin(idx), end(idx),
		[&v, &w](size_t x, size_t y) { // function "bigger than"
			return v[x] / w[x] > v[y] / w[y]; // sorts form bigger to lower
		}
	);*/

	/*for (auto i : idx) {
		for (int j = m[i]; j > 0; j--) {
			if (W > 0) {//
				if (w[i] * j >= W) { //>=
					acc_v += (W / (w[i])) * v[i];
					W = 0;
					j = 0;
					break;
				}
				else {
					acc_v += v[i] * j;
					W -= w[i] * j;
					j = 0;
				}
			}
		}
	}
	return acc_v;
	*/
	int n = m.size();
	double acc_v = 0.0;
	for (int j = k; j < n; ++j) {
		if (W > 0) {//
			if (w[j] * m[j] >= W) { //>=
				acc_v += (W / (w[j])) * v[j];
				break;
			}

			acc_v += v[j] * m[j];
			W -= w[j] * m[j];
		}
	}
	return acc_v;
}


double knapsack_d(
	const vector<double>& v,
	const vector<double>& w,
	const vector<unsigned>& m,
	int k,
	double W
) {
	vector<size_t> idx(w.size());
	for (size_t i = 0; i < idx.size(); i++) idx[i] = i;

	sort(idx.begin(), idx.end(), [&v, &w](size_t x, size_t y) {
		return v[x] / w[x] > v[y] / w[y]; });

	double acc_v = 0.0;
	/*for (auto i : idx) {
		for (int j = m[i]; j > 0; j--) {
			if (w[i] * j <= W) {
				acc_v += v[i] * j;
				W -= w[i] * j;
				j = 0;

			}
		}
	}*/
	

	for (int j = k; j < m.size(); ++j) {
		if (w[j] * m[j] <= W) {
			acc_v += v[j] * m[j];
			W -= w[j] * m[j];
		}
	}


	return acc_v;
}

double weight(const vector<double>& w, size_t k, const vector<unsigned>& x) {
	double acc_w = 0.0;
	for (size_t i = 0; i < k; i++)
		acc_w += x[i] * w[i];
	return acc_w;
}

double value(const vector<double>& v, const vector<unsigned>& x) {
	double r = 0.0;
	for (size_t i = 0; i < v.size(); i++) {
		r += v[i] * x[i];
	}
	return r;
}

float knapsack(const vector<double>& v, const vector<double>& w, const vector<unsigned>& m, double W) {
	using Sol = vector<short>;
	using Node = tuple<double, double, Sol, int>; // value, weight, vector, k
	priority_queue< Node > pq; // A priority_queue is a max-heap

	double best_val;// = knapsack_d(v, w, 0, W);// updating best current solution
	pq.emplace(0.0, 0.0, Sol(v.size()), 0); // insert initial node

	while (!pq.empty()) {

		auto [value, weight, x, k] = pq.top(); // structured auto (c++17)
		pq.pop();
		if (k == v.size()) { // base case
			best_val = max(value, best_val);
			continue;

		}

		for (int j = m[k]; j >= 0; j--) { // expanding
			x[k] = j;

			double new_weight = weight + x[k] * w[k]; // updating weight
			double new_value = value + x[k] * v[k]; // updating value

			if (new_weight <= W) { // is feasible
				 // pessimistic bound
				double pes_bound = new_value + knapsack_d(v, w, m, k + 1, W - new_weight);
				best_val = max(best_val, pes_bound);

				double opt_bound = new_value + knapsack_c(v, w, m, k + 1, W - new_weight);
				if (opt_bound > best_val) // is promising
					pq.emplace(new_value, new_weight, x, k + 1);

			}
		}
		return best_val;

	}
}

void vueltaatras(const vector<double>& v, const vector<double>& w, double W, const vector<unsigned>& m,
	size_t k, vector<unsigned>& x, double& best_v, vector<unsigned> sol, double acc_w, double acc_v
) {
	if (k == x.size()) { // It is a leaf
		best_v = max(best_v, acc_v);
		//cout << "k==x.size" << endl;
		return;
	}
	for (int j = m[k]; j >= 0; j--) {
		x[k] = j;
		double present_w = acc_w + x[k] * w[k];
		double present_v = acc_v + x[k] * v[k];
		if (present_w <= W && present_v + knapsack_c(v, w, m, k + 1, W - present_w) > best_v) { // if it is feasible
			vueltaatras(v, w, W, m, k + 1, x, best_v, sol, present_w, present_v); // expand
		}

	}

}


double vueltaatras1(const vector<double>& v, const vector<double>& w, const vector<unsigned>& m, double W) {

	vector<size_t> idx(v.size()); // index vector
	iota(begin(idx), end(idx), 0);

	sort(begin(idx), end(idx),
		[&v, &w](size_t i, size_t j) {
			return v[i] / w[i] > v[j] / w[j];
		}
	);

	vector<double> s_v(v.size()), s_w(w.size());
	vector<unsigned> s_m(w.size());
	for (size_t i = 0; i < v.size(); i++) {
		s_v[i] = v[idx[i]]; // sorted values
		s_w[i] = w[idx[i]]; // sorted weights
		s_m[i] = m[idx[i]];
	}

	vector<unsigned> x(v.size()), sol(v.size());
	double best_v = knapsack_d(s_v, s_w, s_m, 0, W);
	vueltaatras(s_v, s_w, W, s_m, 0, x, best_v, sol, 0, 0);
	return best_v;

}

void imprimir(double resul, double tiempo, const vector<int> copias) {
	cout << resul << endl;
	for (int i = 0; i < copias.size(); i++)
		cout << copias[i] << " ";
	cout << endl;
	cout << tiempo;
}

int main(int argc, char* argv[]) {
	string nombreFichero = "";

	//leemos argumentos de entrada
	if (!leerArgumentos(argc, argv, nombreFichero)) {
		return 0;
	}
	else {
		int n = -1;
		double T = -1; //num de objetos, maximo de tiempo
		vector<double> v, t;//tiempos
		vector<unsigned> m; // maximo, valor


		nombreFichero = "potter_n15.def";
		//nombreFichero = "mipotter.def.txt";

		if (!leerFichero(nombreFichero, n, T, t, v, m))
			return 0;
		vector<int> copias(n, 0);//contador de copias de cada objeto
		double tiempoTotal = 0;
		double resul = -1;
		vector<unsigned> sol;
		clock_t start = clock();
		//resul= vueltaatras1(v, t, T, sol);
		//knapsack(v, m, t, T);
		vector<size_t> idx(v.size()); // index vector
		iota(begin(idx), end(idx), 0);

		sort(begin(idx), end(idx),
			[&v, &t](size_t i, size_t j) {
				return v[i] / t[i] > v[j] / t[j];
			}
		);

		vector<double> s_v(v.size()), s_t(t.size());
		vector<unsigned> s_m(t.size());
		for (size_t i = 0; i < v.size(); i++) {
			s_v[i] = v[idx[i]]; // sorted values
			s_t[i] = t[idx[i]]; // sorted weights
			s_m[i] = m[idx[i]];
		}
		resul = knapsack_d(s_v, s_t, s_m, 0, T);
		cout << "_d:" << resul << endl;
		resul = knapsack_c(s_v, s_t, s_m, 0, T);
		cout << "_c:" << resul << endl;
		resul = vueltaatras1(v, t, m, T);
		cout << "vueltaatras: " << resul << endl;
		clock_t end = clock();
		cout << (double(end - start) / ((clock_t)1000)) << "s" << endl;

	}
	return 0;
}
