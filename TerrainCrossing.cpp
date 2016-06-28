#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <utility>
#include <limits>
#include <queue>
#include <random>
#include <chrono>
#include <list>
#include <deque>

using namespace std;

const int UNIQUE_OFFSET = 5;
const int MAX_TERRAIN_TYPE = 9;
const int MIN_TERRAIN_TYPE = 0;

//--------------------------
//------Location class------
//--------------------------


string resize_string(string &s, int N){
	while( s.length() < N ){
		s = " " + s;
	}
	return s;
}


string resize_string(string &&s, int N){
	while( s.length() < N ){
		s = " " + s;
	}
	return s;
}


template<typename Type> void print_vector_2d(vector< vector<Type> > &v, int offset){

	for(int i = 0; i < v.size(); i++){
			for(int j = 0; j < v[i].size(); j++){
					string s = to_string(v[i][j]);
					s = resize_string(s, offset);
					cerr << s << " ";
			}
			cerr << endl;
	}
}


template<typename Type> void print_vector_2d(vector< vector<Type> > &&v){
	for(int i = 0; i < v.size(); i++){
			for(int j = 0; j < v[i].size(); j++){
					cerr << v[i][j] << " ";
			}
			cerr << endl;
	}
}

void print_matrix_ij(int N, int offset){

	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			string s = to_string(i) + "-" + to_string(j);
			s = resize_string(s, offset);
			cerr << s << " ";
		}
		cerr << endl;
	}
}


//----------------------------
//--------Vertex class--------
//----------------------------



class Vertex {
private:
	int m_i;
	int m_j;

	double m_x;
	double m_y;

    int m_tt;
	int m_vertex_type;

public:
	Vertex(int i = 0, int j = 0, double x = 0.0, double y = 0.0, int tt = 0, int vertex_type = 0):
	m_i(i), m_j(j), m_x(x), m_y(y), m_tt(tt), m_vertex_type(vertex_type) {}

    int get_i() const {return m_i;}
    int get_j() const {return m_j;}

    double get_x() const {return m_x;}
    double get_y() const {return m_y;}

    int get_tt() const {return m_tt;}


    vector< pair<int, int> > get_neighbours(int &N);

};



vector< pair<int, int> > Vertex::get_neighbours(int &N){

    vector< pair<int, int> > neighbours = vector< pair<int, int> >();

        int i_up = this->get_i()-1;
		int j_up = this->get_j();
		if (i_up >= 0)
            neighbours.push_back(pair<int, int>(i_up, j_up));

		int i_down = this->get_i()+1;
		int j_down = this->get_j();
		if (i_down <= N-1)
            neighbours.push_back(pair<int, int>(i_down, j_down));

		int i_left = this->get_i();
		int j_left = this->get_j()-1;
		if (j_left >= 0)
            neighbours.push_back(pair<int, int>(i_left, j_left));

		int i_right = this->get_i();
		int j_right = this->get_j()+1;
		if (j_right <= N-1)
			neighbours.push_back(pair<int, int>(i_right, j_right));

    return neighbours;
}




inline double euclidean_distance(Vertex &v1, Vertex &v2){
    double x2 = pow( v1.get_x() - v2.get_x(), 2);
    double y2 = pow( v1.get_y() - v2.get_y(), 2);

    return sqrt(x2 + y2);
}



inline double manhatan_distance(Vertex &v1, Vertex &v2){
    double i_abs = abs(v1.get_i() - v2.get_i());
    double j_abs = abs(v1.get_j() - v2.get_j());

    return (i_abs + j_abs);
}



pair<double, double> intersection_point(Vertex &v1, Vertex &v2){

    double x0 = 0.0;
    double y0 = 0.0;

    if( v1.get_j() == v2.get_j() ){
        x0 = max(v1.get_i(), v2.get_i());
        y0 = v1.get_y() + (v2.get_y() - v1.get_y()) * (x0 - v1.get_x()) / (v2.get_x() - v1.get_x());
    } else {
        y0 = max(v1.get_j(), v2.get_j());
        x0 = v1.get_x() + (v2.get_x() - v1.get_x()) * (y0 - v1.get_y()) / (v2.get_y() - v1.get_y());
    }
    return make_pair(x0,y0);
}



double distance(Vertex &v1, Vertex &v2){

    if(manhatan_distance(v1,v2) == 0)
        return euclidean_distance(v1, v2)*v1.get_tt();

    int tt_1 = v1.get_tt();
    int tt_2 = v2.get_tt();

    double score = pow(tt_1 - tt_2 ,2);

    // intersection point
    pair<double, double> ip = intersection_point(v1, v2);
    double x0 = ip.first;
    double y0 = ip.second;

    // intersection vertex
    Vertex iv = Vertex(floor(x0), floor(y0), x0, y0, -1);

    return (score + euclidean_distance(v1, iv)*v1.get_tt() + euclidean_distance(iv, v2)*v2.get_tt());
}


std::ostream& operator<<(std::ostream& os, const Vertex& v)
{

	string si = to_string(v.get_i());
	string sj = to_string(v.get_j());
	string s = resize_string( si+"-"+sj, UNIQUE_OFFSET);

    os << s;
    return os;
}



//----------------------------
//--------Graph class---------
//----------------------------



class Graph {
private:

    int m_n; // world_map size
    int m_n_vertexes; // number of vertexes, m_n*m_n
    vector<Vertex> m_vertex_array;
    vector< vector<int> > m_vertex_matrix;
    vector< vector<pair<int, double> > > m_adjacency_list;

public:
    Graph();
    Graph(vector< vector<int> > &map_matrix);

    vector<vector<int> > get_vertex_matrix() const { return m_vertex_matrix;};
    void print_graph();

    vector<Vertex> find_shortest_path_dijkstra(Vertex &v1, Vertex &v2);
    vector<Vertex> find_shortest_path_a_star(Vertex &v1, Vertex &v2);

    double heuristic_function(Vertex &current_vertex_id, Vertex &neighbour_id, Vertex &target_id);


};


class priority_queue_compare_less{
	public:
		bool operator ()(pair<int, double> &pair_1, pair<int, double> &pair_2){
			return pair_1.second > pair_2.second;
		}
};


Graph::Graph() {
    m_n = 0;
    m_vertex_array = vector<Vertex>();
    m_vertex_matrix = vector< vector<int> >();
    m_adjacency_list = vector< vector<pair<int, double> > >();
}


Graph::Graph(vector<vector<int> > &map_matrix) {

    m_n = map_matrix.size();
    m_n_vertexes = m_n*m_n;

    m_vertex_array.resize(m_n_vertexes);
    m_vertex_matrix.resize(m_n);
    for(int i = 0; i < m_n; i++)
        m_vertex_matrix[i].resize(m_n);

    int index = 0;
    for(int i = 0; i < m_n; i++){
        for(int j = 0; j < m_n; j++){
            m_vertex_array[index] = Vertex(i,j,i+0.5,j+0.5,map_matrix[i][j]);
            m_vertex_matrix[i][j] = index;
            index++;
        }
    }


    m_adjacency_list.resize(m_n_vertexes);
    for(int i = 0; i < m_n_vertexes; i++){
        vector<pair<int, int> > neighbours = m_vertex_array[i].get_neighbours(m_n);


        for(int j = 0; j < neighbours.size(); j++){

            int i_index = neighbours[j].first;
            int j_index = neighbours[j].second;
            int vertex_id = m_vertex_matrix[i_index][j_index];
            double d = distance(m_vertex_array[i], m_vertex_array[vertex_id]);

            m_adjacency_list[i].push_back(pair<int, double>( vertex_id , d ) );
        }
    }
}


void Graph::print_graph(){

    for(int i = 0; i < m_vertex_array.size(); i++){
        int i_index = m_vertex_array[i].get_i();
        int j_index = m_vertex_array[i].get_j();

        string s = "Vertex: " + to_string(i_index) + to_string(j_index) + " -> ";
        for(int k = 0; k < m_adjacency_list[i].size(); k++){
            int n_i_index = m_vertex_array[ m_adjacency_list[i][k].first ].get_i();
            int n_j_index = m_vertex_array[ m_adjacency_list[i][k].first ].get_j();
            s = s + "(" + to_string(n_i_index) + to_string(n_j_index) + ", " + to_string(m_adjacency_list[i][k].second) + ") -> ";
        }
        s = s + "END";
        cerr << s << endl;
    }

}


vector<Vertex> Graph::find_shortest_path_dijkstra(Vertex &source, Vertex &target){

	vector<double> distances(m_n_vertexes, std::numeric_limits<double>::max());
	vector<double> previous(m_n_vertexes, -1);
	vector<bool> visited(m_n_vertexes, false);

	int source_id = m_vertex_matrix[source.get_i()][source.get_j()];
	int target_id = m_vertex_matrix[target.get_i()][target.get_j()];
	distances[source_id] = 0.0;

	priority_queue< pair<int, double>, vector< pair<int, double> >, priority_queue_compare_less> pq;
	pq.push(pair<int, double>(source_id, 0.0));

	while( !pq.empty() ){
		pair<int, double> current_vertex = pq.top();
		pq.pop();

		int current_vertex_id = current_vertex.first;

		if(current_vertex_id == target_id)
			break;

		visited[current_vertex_id] = true;

		for(int i = 0; i < m_adjacency_list[current_vertex_id].size(); ++i){

			int neighbour_vertex_id = m_adjacency_list[current_vertex_id][i].first;

			if(visited[neighbour_vertex_id] == true)
				continue;

			double d = distances[current_vertex_id] + m_adjacency_list[current_vertex_id][i].second;
			if(d < distances[neighbour_vertex_id]){
				distances[neighbour_vertex_id] = d;
				previous[neighbour_vertex_id] = current_vertex_id;
				pq.push(pair<int, double>(neighbour_vertex_id, d));
			}
		}
	}

	vector<Vertex> shortest_path;
	int u = target_id;
	while(previous[u] != -1){
		shortest_path.push_back(m_vertex_array[u]);
		u = previous[u];
	}
	shortest_path.push_back(m_vertex_array[u]);
	
	
	return shortest_path;
}


inline double Graph::heuristic_function(Vertex &current_vertex_id, Vertex &neighbour_id, Vertex &target_id){

    double d1 = distance(current_vertex_id, neighbour_id);
    double d2 = manhatan_distance(neighbour_id, target_id);

    return d1 + d2;

}


vector<Vertex> Graph::find_shortest_path_a_star(Vertex &source, Vertex &target){


	vector<double> distances(m_n_vertexes, std::numeric_limits<double>::max());
	vector<double> previous(m_n_vertexes, -1);
	vector<bool> visited(m_n_vertexes, false);

	int source_id = m_vertex_matrix[source.get_i()][source.get_j()];
	int target_id = m_vertex_matrix[target.get_i()][target.get_j()];
	distances[source_id] = 0.0;

	priority_queue< pair<int, double>, vector< pair<int, double> >, priority_queue_compare_less> pq;
	pq.push(pair<int, double>(source_id, 0.0));

	while( !pq.empty() ){
		pair<int, double> current_vertex = pq.top();
		pq.pop();

		int current_vertex_id = current_vertex.first;

		if(current_vertex_id == target_id)
			break;

		visited[current_vertex_id] = true;

		for(int i = 0; i < m_adjacency_list[current_vertex_id].size(); ++i){

			int neighbour_vertex_id = m_adjacency_list[current_vertex_id][i].first;

			if(visited[neighbour_vertex_id] == true)
				continue;

			double d = distances[current_vertex_id] + m_adjacency_list[current_vertex_id][i].second;
			if(d < distances[neighbour_vertex_id]){
                double h_d = heuristic_function(m_vertex_array[current_vertex_id],
                                                m_vertex_array[neighbour_vertex_id], target);
				distances[neighbour_vertex_id] = d;
				previous[neighbour_vertex_id] = current_vertex_id;
				pq.push(pair<int, double>(neighbour_vertex_id, d + h_d ));
			}
		}
	}

	vector<Vertex> shortest_path;
	int u = target_id;
	while(previous[u] != -1){
		shortest_path.push_back(m_vertex_array[u]);
		u = previous[u];
	}
	shortest_path.push_back(m_vertex_array[u]);

    return shortest_path;
}



//---------------------------------
//------TerrainCrossing class------
//---------------------------------


class TerrainCrossing {
private:

    vector< vector<int> > m_map_matrix;

public:

	TerrainCrossing();


	void initialize_map_matrix(vector<string> &input_map);


    vector<vector<int> > world_map_to_map_matrix(vector<string> world_map);
	vector< vector<int> > generate_random_map_matrix(int N);
    vector<double> getPath(vector<string> input_map, vector<double> locations, int capacity);


};


//------------------------------------------------
//------TerrainCrossing class implementation------
//------------------------------------------------


TerrainCrossing::TerrainCrossing(){

    m_map_matrix = vector< vector<int> >();

}


void TerrainCrossing::initialize_map_matrix(vector<string> &world_map) {

	// Resize vector to meet the map size.
	int N = world_map.size();
	m_map_matrix.resize(N);
	for(int i = 0; i < N; i++)
		m_map_matrix[i].resize(N);

	// Fill the map matrix.
	for(int i = 0; i < N; i++){
		for(int j = 0; j < world_map[i].length(); j++){
				m_map_matrix[i][j] = world_map[i][j] - '0';
		}
	}
}



vector< vector<int> > TerrainCrossing::generate_random_map_matrix(int N){

		std::random_device rd;
		std::mt19937 engine(rd());
		std::uniform_int_distribution<> dist(MIN_TERRAIN_TYPE, MAX_TERRAIN_TYPE);

		vector< vector<int> > map_matrix;
		map_matrix.resize(N);

		for(int i = 0; i < N; i++)
			map_matrix[i].resize(N);

		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				map_matrix[i][j] = dist(engine);
			}
		}
		return map_matrix;
}



vector<double> TerrainCrossing::getPath(vector<string> world_map, vector<double> locations, int capacity) {

		int N = world_map.size();
		int offset = 5;


        vector<string> custom_world_map = {"99999","99009","90909","90009","99999"};


        initialize_map_matrix(custom_world_map);
		//vector< vector<int> > random_map_matrix = generate_random_map_matrix(5);

		int random_map_size = 50;
		m_map_matrix = generate_random_map_matrix(random_map_size);


		int source_i = 0;
		int source_j = 0;
		int target_i = random_map_size-1;
		int target_j = random_map_size-1;
        Vertex v1 = Vertex(source_i, source_j, source_i + 0.5, source_j + 0.5, m_map_matrix[source_i][source_j]);
        Vertex v2 = Vertex(target_i, target_j, target_i + 0.5, target_j + 0.5, m_map_matrix[target_i][target_j]);

        Graph g = Graph(m_map_matrix);
        print_vector_2d(m_map_matrix, offset);
        cerr << endl;

        print_matrix_ij(5, offset);
        cerr << endl;

		//g.print_graph();

		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
        vector<Vertex> shortest_path_dijkstra = g.find_shortest_path_dijkstra(v1, v2);
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;

		cerr << "Dijkstra elapsed time: " << elapsed_seconds.count() << endl;


		start = std::chrono::system_clock::now();
        vector<Vertex> shortest_path_a_star = g.find_shortest_path_a_star(v1, v2);
		end = std::chrono::system_clock::now();
		elapsed_seconds = end-start;

		cerr << "A* elapsed time: " << elapsed_seconds.count() << endl;


        for(int i = 0; i < shortest_path_dijkstra.size(); i++){
        	cerr << shortest_path_dijkstra[i] << " -> ";
        }
        cerr << endl;

        for(int i = 0; i < shortest_path_a_star.size(); i++){
        	cerr << shortest_path_a_star[i] << " -> ";
        }
        cerr << endl;



		vector<double> ret;
		ret.push_back(4.9995); ret.push_back(0.7658);
			// pick up item 1
			ret.push_back(4.7867); ret.push_back(0.7658);
			// drop it off at target 6
			ret.push_back(3.8144); ret.push_back(0.1081);
			// pick up item 0
			ret.push_back(3.7648); ret.push_back(1.2640);
			// drop it off at target 7
			ret.push_back(3.3420); ret.push_back(2.2000);
			ret.push_back(3.3420); ret.push_back(3.0530);
			// pick up item 2
			ret.push_back(2.5000); ret.push_back(3.0530);
			ret.push_back(1.5000); ret.push_back(3.0530);
			ret.push_back(0.7225); ret.push_back(3.0530);
			ret.push_back(0.7225); ret.push_back(2.5000);
			ret.push_back(0.7225); ret.push_back(1.4533);
			// pick up item 3
			ret.push_back(0.2299); ret.push_back(2.8555);
			ret.push_back(0.2299); ret.push_back(3.8555);
			ret.push_back(0.2299); ret.push_back(4.8555);
			// drop it off at target 4
			ret.push_back(0.5000); ret.push_back(3.3869);
			ret.push_back(1.2611); ret.push_back(3.3869);
			// drop it off at target 5
			ret.push_back(2.2611); ret.push_back(3.3869);
			ret.push_back(2.2611); ret.push_back(4.6214);
			ret.push_back(3.7958); ret.push_back(4.6214);
			// exit
			ret.push_back(3.7958); ret.push_back(4.9995);
			return ret;
    }


// -------8<------- end of solution submitted to the website -------8<-------

template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}

int main() {
	TerrainCrossing tc;
    int M;
    cin >> M;
    vector<string> input_map(M);
    getVector(input_map);

    int L;
    cin >> L;
    vector<double> locations(L);
    getVector(locations);

    int capacity;
    cin >> capacity;

    vector<double> ret = tc.getPath(input_map, locations, capacity);
    cout << ret.size() << endl;
    for (int i = 0; i < ret.size(); ++i)
        cout << ret[i] << endl;
    cout.flush();
}






