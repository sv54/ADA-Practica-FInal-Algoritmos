#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<numeric>
#include <queue>
#include <tuple>
#include <ctime>


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
		for (int i = 0; i < n; i++) {
			f >> aux;
			t.push_back(aux);
		}
		for (int i = 0; i < n; i++) {
			f >> aux;
			v.push_back(aux);
		}
		for (int i = 0; i < n; i++) {
			f >> aux;
			m.push_back(aux);
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

ostream& operator<<(ostream& os, const vector<double>& vect) {
	for (unsigned i = 0; i < vect.size(); i++) {
		os << vect[i] << " ";
	}
	return os;
}

double knapsack_c(
	const vector<double>& v, // values
	const vector<double>& t, // time
	const vector<unsigned>& m,
	unsigned k,
	double T // time limit
) {
	int n = m.size();
	double acc_v = 0.0;
	for (int j = k; j < n; ++j) {
		if (T > 0) {
			if (t[j] * m[j] >= T) {
				acc_v += (T / (t[j])) * v[j];
				break;
			}

			acc_v += v[j] * m[j];
			T -= t[j] * m[j];
		}
	}
	return acc_v;
}


double knapsack_d(
	const vector<double>& v,
	const vector<double>& t,
	const vector<unsigned>& m,
	int k,
	double T
) {
	int n = m.size();
	double acc_v = 0.0;

	for (int j = k; j < n; ++j) {
		for (int i = m[j]; i > 0; --i) {
			if (t[j] * i <= T) {
				acc_v += v[j] * i;
				T -= t[j] * i;
				i = 0;
			}
		}
	}

	return acc_v;
}


double time(const vector<double>& t, size_t k, const vector<unsigned>& x) {
	double acc_t = 0.0;
	for (size_t i = 0; i < k; i++)
		acc_t += x[i] * t[i];
	return acc_t;
}

double value(const vector<double>& v, const vector<unsigned>& x) {
	double r = 0.0;
	for (size_t i = 0; i < v.size(); i++) {
		r += v[i] * x[i];
	}
	return r;
}
void vueltaatras(const vector<double>& v, const vector<double>& t, double T, const vector<unsigned>& m,
	size_t k, vector<unsigned>& x, double& best_v, vector<unsigned> sol, double acc_t, double acc_v, vector<unsigned>& estadisticas
) {
	if (k == x.size()) { // It is a leaf
		estadisticas[2] = estadisticas[2] + 1;
		best_v = max(best_v, acc_v);
		//cout << "k==x.size" << endl;
		return;
	}
	for (int j = m[k]; j >= 0; j--) {
		estadisticas[0] = estadisticas[0] + 1;
		x[k] = j;
		double present_t = acc_t + x[k] * t[k];
		double present_v = acc_v + x[k] * v[k];
		if (present_t <= T) { //factible?
			if (present_v + knapsack_c(v, t, m, k + 1, T - present_t) > best_v) { // prometedor?
				estadisticas[1] = estadisticas[1] + 1;
				vueltaatras(v, t, T, m, k + 1, x, best_v, sol, present_t, present_v, estadisticas); // expand
			}
			else {
				estadisticas[4] = estadisticas[4] + 1;
			}
		}
		else {
			estadisticas[3] = estadisticas[3] + 1;
		}
	}

}

double add_rest(const vector<double>& v, const vector<unsigned>& m, size_t k) {
	double r = 0.0;
	for (size_t i = k; i < v.size(); i++) r += v[i] * m[i];
	return r;
}


double vueltaatras1(const vector<double>& v, const vector<double>& t, const vector<unsigned>& m, double T, vector<unsigned>& estadisticas) {

	vector<unsigned> x(v.size()), sol(v.size());
	double best_v = knapsack_d(v, t, m, 0, T);
	vueltaatras(v, t, T, m, 0, x, best_v, sol, 0, 0, estadisticas);
	return best_v;

}


float knapsack1(const vector<double>& v, const vector<double>& t, const vector<unsigned>& m, double T, vector<unsigned>& estadisticas) {
	using Sol = vector<short>;
	using Node = tuple<double, double, Sol, int>; // value, time, vector, k
	priority_queue< Node > pq; // A priority_queue is a max-heap

	double best_val = knapsack_d(v, t, m, 0, T);
	//double opt_bound = knapsack_c(v, t, m, 0, T);
	pq.emplace(0.0, 0.0, Sol(v.size()), 0); // insert initial node

	while (!pq.empty()) {
		estadisticas[1] = estadisticas[1] + 1; //nodos visitados
		auto [value, time, x, k] = pq.top(); // structured auto (c++17)
		pq.pop();

		if (value < best_val)
			estadisticas[5] = estadisticas[5] + 1; //eran prometedores pero fueron descartados

		if (k == v.size()) { // base case

			estadisticas[2] = estadisticas[2] + 1;// nodo hoja visitados
			if (value > best_val) {
				best_val = value;
				estadisticas[6] = estadisticas[6] + 1; //actualizado a partir nodo completado
			}
			continue;

		}

		for (int j = m[k]; j >= 0; j--) { // expanding
			x[k] = j;
			estadisticas[0] = estadisticas[0] + 1; //nodos explorados

			double new_time = time + x[k] * t[k]; // updating time
			double new_value = value + x[k] * v[k]; // updating value

			if (new_time <= T) { // es factible?

				// cota pesimista, no se utiliza
				//double pes_bound = new_value + knapsack_d(v, t, m, k + 1, T - new_time);
				//if (pes_bound > best_val) {
				//	best_val = pes_bound;
				//	estadisticas[7] = estadisticas[7] + 1; //actualizado a partir de la cota pesimista
				//}

				double opt_bound = new_value + knapsack_c(v, t, m, k + 1, T - new_time);
				if (opt_bound > best_val) { // is promising
					pq.emplace(new_value, new_time, x, k + 1);
				}
				else {
					estadisticas[4] = estadisticas[4] + 1; //nodos descartados por no ser prometedores
				}
			}
			else {
				estadisticas[3] = estadisticas[3] + 1; //nodos descartados por no ser factibles
			}
		}
		


	}
	return best_val;
}

float knapsack(const vector<double>& v, const vector<double>& t, const vector<unsigned>& m, double T, vector<unsigned>& estadisticas) {
	using Sol = vector<short>;
	using Node = tuple<double, double, double, Sol, int>; // cota optimista, value, time, vector, k
	priority_queue< Node > pq; // A priority_queue is a max-heap

	double best_val = knapsack_d(v, t, m, 0, T);
	double opt_bound = knapsack_c(v, t, m, 0, T);
	pq.emplace(opt_bound, 0.0, 0.0, Sol(v.size()), 0); // insert initial node

	while (!pq.empty()) {
		estadisticas[1] = estadisticas[1] + 1; //nodos visitados
		auto [opt_bound1, value, time, x, k] = pq.top(); // structured auto (c++17)
		pq.pop();

		if (opt_bound < best_val)
			estadisticas[5] = estadisticas[5] + 1; //eran prometedores pero fueron descartados
		else{
		if (k == v.size()) { // base case

			estadisticas[2] = estadisticas[2] + 1;// nodo hoja visitados
			if (value > best_val) {
				best_val = value;
				estadisticas[6] = estadisticas[6] + 1; //actualizado a partir nodo completado
			}
			continue;

		}

		for (int j = m[k]; j >= 0; j--) { // expanding
			x[k] = j;
			estadisticas[0] = estadisticas[0] + 1; //nodos explorados

			double new_time = time + x[k] * t[k]; // updating time
			double new_value = value + x[k] * v[k]; // updating value

			if (new_time <= T) { // es factible?

				// cota pesimista, no se utiliza
				//double pes_bound = new_value + knapsack_d(v, t, m, k + 1, T - new_time);
				//if (pes_bound > best_val) {
				//	best_val = pes_bound;
				//	estadisticas[7] = estadisticas[7] + 1; //actualizado a partir de la cota pesimista
				//}

				double opt_bound = new_value + knapsack_c(v, t, m, k + 1, T - new_time);
				if (opt_bound > best_val) // is promising
					pq.emplace(opt_bound, new_value, new_time, x, k + 1);
				else {
					estadisticas[4] = estadisticas[4] + 1; //nodos descartados por no ser prometedores
				}
			}
			else {
				estadisticas[3] = estadisticas[3] + 1; //nodos descartados por no ser factibles
			}
		}
		}


	}
	return best_val;
}




int main(int argc, char* argv[]) {
	string nombreFichero = "";

	//leemos argumentos de entrada
	if (!leerArgumentos(argc, argv, nombreFichero)) {
		return 0;
	}
	else {
		int n = -1; //num objetos
		double T = -1; //maximo de tiempo
		vector<double> v, t;//valor, tiempos
		vector<unsigned> m; // maximo de copias
		nombreFichero = "potter_n30.def";
		if (!leerFichero(nombreFichero, n, T, t, v, m))
			return 0;

		double resul = -1;
		//clock_t start = clock();

		vector<size_t> idx(v.size()); // index vector
		iota(begin(idx), end(idx), 0);

		sort(begin(idx), end(idx),
			[&v, &t](size_t i, size_t j) {
				return v[i] / t[i] > v[j] / t[j];
			}
		);

		vector<double> s_v(v.size()), s_t(t.size());
		vector<unsigned> s_m(t.size()), estadisticasRyP(8), estadisticasVA(5);
		for (size_t i = 0; i < v.size(); i++) {
			s_v[i] = v[idx[i]]; // sorted values
			s_t[i] = t[idx[i]]; // sorted times
			s_m[i] = m[idx[i]];	// cantidad de copias
		}

		double start = clock();
		resul = knapsack1(s_v, s_t, s_m, T, estadisticasRyP); //ramificacion y poda
		double end = clock();


		cout << resul << endl;
		for (int i = 0; i < estadisticasRyP.size()-1; i++) {
			cout << estadisticasRyP[i] << " ";
		}
		//Al no calcular cotas pesimistas no se incluye
		//el ultimo valor en las estadisticas
		cout << "-" << endl;

		cout << ((end - start) / (double)CLOCKS_PER_SEC) * 1000.0 << endl;

	}
	return 0;
}
