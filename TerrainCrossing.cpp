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
#include <unordered_map>
#include <tuple>
#include <array>

using namespace std;

const int UNIQUE_OFFSET = 5;
const int MAX_TERRAIN_TYPE = 9;
const int MIN_TERRAIN_TYPE = 0;
const double EPSILON = 5E-4;
const double TWO_TIMES_EPSILON = 0.001;

//-----------------------------
//------Support functions------
//-----------------------------


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
				cerr << v[j][i] << " ";
		}
		cerr << endl;
	}
}

void print_matrix_ij(int N, int offset){

	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			string s = to_string(j) + "-" + to_string(i);
			s = resize_string(s, offset);
			cerr << s << " ";
		}
		cerr << endl;
	}
}


bool double_compare(double &d1, double &d2){
	return abs(d1 - d2) < 1E-6;
}


bool double_compare(double &&d1, double &&d2){
	return abs(d1 - d2) < 1E-6;
}


vector<double> linspace(double a, double b, int N){

	vector<double> grid(N);
	double step = (b-a)/(double(N-1));

	for(int i = 0; i < N; i++){
		grid[i] = a + i*step;
	}
	return grid;
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
	int m_vertex_type; // 0 - transfer vertex, 1 - item vertex, 2 - target vertex

	int m_id;

	bool m_visited;


public:
	Vertex(int i = 0, int j = 0, double x = 0.0, double y = 0.0, int tt = 0, int vertex_type = 0, int id = -1):
	m_i(i), m_j(j), m_x(x), m_y(y), m_tt(tt), m_vertex_type(vertex_type), m_id(id) {
		m_visited = false;
	}

	int get_id() const {return m_id;}

    int get_i() const {return m_i;}
    int get_j() const {return m_j;}

    double get_x() const {return m_x;}
    double get_y() const {return m_y;}

    void set_x(double x) { m_x = x;}
    void set_y(double y) { m_y = y;}

    int get_tt() const {return m_tt;}
	int get_vertex_type() const {return m_vertex_type;}

	bool get_visited() const {return m_visited;}

	void set_as_visited() {m_visited = true;}

    vector< pair<int, int> > get_neighbours(int &N);


};



vector< pair<int, int> > Vertex::get_neighbours(int &N){

    vector< pair<int, int> > neighbours = vector< pair<int, int> >();

	if (this->get_vertex_type() != 0)
		neighbours.push_back(pair<int, int>(this->get_i(), this->get_j()));

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



double triple_path_distance(Vertex &v1, Vertex &v2, Vertex &v3){

	return distance(v1,v2) + distance(v2, v3);
}



std::ostream& operator<<(std::ostream& os, const Vertex& v) {

	string sindex = to_string(v.get_id());
	string si = to_string(v.get_i());
	string sj = to_string(v.get_j());
	string sx = to_string(v.get_x());
	string sy = to_string(v.get_y());
	string stt = to_string(v.get_tt());
	string svt = to_string(v.get_vertex_type());

	string s = resize_string( si+"-"+sj+" x: "+sx+" y: "+sy+" tt: "+stt+" vertex type: "+svt+" index: "+sindex, UNIQUE_OFFSET);

    os << s << endl;
    return os;
}


class VertexCompare {
public:
	bool operator()(const Vertex& v1, const Vertex& v2){

		bool i_true = (v1.get_i() < v2.get_i());
		bool j_true = (v1.get_j() < v2.get_j());

		bool x_true = (v1.get_x() < v2.get_x());
		bool y_true = (v1.get_y() < v2.get_y());

		bool tt_true = (v1.get_tt() < v2.get_tt());
		bool vertex_type_true = (v1.get_vertex_type() < v2.get_vertex_type());

		return (i_true || j_true || x_true || y_true || tt_true || vertex_type_true);

	}

};




//-----------------------------------------------------
//--------Vertex class ofr stochastic approach---------
//-----------------------------------------------------




class StochasticVertex {
private:
    int m_vertex_id;
    int m_capacity;
    StochasticVertex* m_next;
    StochasticVertex* m_prev;

    double m_distance_to_prev;
    int m_vertex_type;
    
public:
    StochasticVertex(int id = -1, int c = 0, StochasticVertex* vn = nullptr, StochasticVertex* vp = nullptr, double d = 0.0, int t = 0) :
                     m_vertex_id(id), m_next(vn), m_prev(vp), m_distance_to_prev(d), m_vertex_type(t) {}

    void set_id(int id) {m_vertex_id = id;}
    void set_capacity(int c) {m_capacity = c;}
    void set_next(StochasticVertex *vn) {m_next = vn;}
    void set_prev(StochasticVertex *vp) {m_prev = vp;}
    void set_distance_to_prev(double d) {m_distance_to_prev = d;}
    void set_vertex_type(int t) {m_vertex_type = t;}

    int get_id() const {return m_vertex_id; }
    int get_capacity() const {return m_capacity;}
    StochasticVertex* get_next() const {return m_next;}
    StochasticVertex* get_prev() const {return m_prev;}
    double get_distance_to_prev() const {return m_distance_to_prev;}
    int get_vertex_type() const {return m_vertex_type;}
};





//----------------------------
//--------Graph class---------
//----------------------------



class Graph {
private:

    int m_n; // world_map size
    int m_n_vertexes; // number of vertexes, m_n*m_n
    int m_items;
    vector<Vertex> m_vertex_array;
    vector< vector<int> > m_vertex_matrix;
    vector< vector<pair<int, double> > > m_adjacency_list;

	vector<double> m_distances;
	vector<double> m_previous;
	vector<bool> m_visited;


public:
    Graph();
    Graph(vector< vector<int> > &map_matrix, vector<double> &locations);

	Vertex& get_vertex_in_graph(int &vertex_id) {return m_vertex_array[vertex_id];}
	Vertex& get_vertex_in_graph(int &&vertex_id) {return m_vertex_array[vertex_id];}
	
    vector<vector<int> > get_vertex_matrix() const { return m_vertex_matrix;};
    int return_vertex_id(int &i, int &j) {return m_vertex_matrix[i][j];}
    double return_vertex_x(int &vertex_id) {return m_vertex_array[vertex_id].get_x();}
    double return_vertex_y(int &vertex_id) {return m_vertex_array[vertex_id].get_y();}

    void print_graph();
	vector<Vertex> find_shortest_path_a_star(int &source_id, int &target_id);
	vector<Vertex> find_shortest_path_a_star(int &&source_id, int &&target_id);

    double path_cost(vector<Vertex> &path);
    double path_cost(vector<Vertex> &&path);
    double many_path_cost(vector< vector<Vertex> > &path_vectors);

    double heuristic_function(Vertex &current_vertex_id, Vertex &neighbour_id, Vertex &target_id);

};


//-------------------------------------------
//--------Graph class implementation---------
//-------------------------------------------


class PriorityQueueCompareLess{
	public:
		inline bool operator ()(const pair<int, double> &pair_1, const pair<int, double> &pair_2) const {
			return pair_1.second > pair_2.second;
		}
};


Graph::Graph() {
    m_n = 0;
    m_vertex_array = vector<Vertex>();
    m_vertex_matrix = vector< vector<int> >();
    m_adjacency_list = vector< vector<pair<int, double> > >();
}


Graph::Graph(vector<vector<int> > &map_matrix, vector<double> &locations) {

	m_items = locations.size()/4;

    m_n = map_matrix.size();
    m_n_vertexes = m_n*m_n + 2*m_items;

    m_vertex_array.resize(m_n_vertexes);
    m_vertex_matrix.resize(m_n);
    for(int i = 0; i < m_n; i++)
        m_vertex_matrix[i].resize(m_n);

	// add item vertexes to graph
    int index = 0;
	for(int i = 0; i < m_items; i = i + 1){
		double x_item = locations[2*i];
		double y_item = locations[2*i + 1];

		int i_item = floor(x_item);
		int j_item = floor(y_item);

		m_vertex_array[index] = Vertex(i_item, j_item, x_item, y_item, map_matrix[i_item][j_item], 1, index);
		index++;
	}

	// add target vertexes to graph
	for(int i = 0; i < m_items; i = i + 1){
		double x_target = locations[2*i + 2*m_items];
		double y_target = locations[2*i + 1 + 2*m_items];

		int i_target = floor(x_target);
		int j_target = floor(y_target);

		m_vertex_array[index] = Vertex(i_target, j_target, x_target, y_target, map_matrix[i_target][j_target], 2, index);
		index++;
	}


	for(int i = 0; i < m_n; i++){
        for(int j = 0; j < m_n; j++){
            m_vertex_array[index] = Vertex(i,j,i+0.5,j+0.5,map_matrix[i][j], 0, index);
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

			// item and target vertexes have five neighbours including the one that is located in their square
			if(m_vertex_array[i].get_vertex_type() != 0)
				m_adjacency_list[vertex_id].push_back(pair<int, double>( i , d ) );
		}
    }

    m_distances.resize(m_n_vertexes, std::numeric_limits<double>::max());
    m_previous.resize(m_n_vertexes, -1);
    m_visited.resize(m_n_vertexes, false);

}



void Graph::print_graph(){

    for(int i = 0; i < m_vertex_array.size(); i++){
        int i_index = m_vertex_array[i].get_i();
        int j_index = m_vertex_array[i].get_j();

		int vertex_id = m_vertex_array[i].get_id();
        string s = "Vertex: " + to_string(i_index) + "-" + to_string(j_index) + "," + to_string(vertex_id) + " -> ";
        for(int k = 0; k < m_adjacency_list[i].size(); k++){
            int n_i_index = m_vertex_array[ m_adjacency_list[i][k].first ].get_i();
            int n_j_index = m_vertex_array[ m_adjacency_list[i][k].first ].get_j();
            s = s + "(" + to_string(n_i_index) + "-" + to_string(n_j_index) + ", " + to_string(m_adjacency_list[i][k].second) + ", " + to_string(m_adjacency_list[i][k].first) + ") -> ";
        }
        s = s + "END";
        cerr << s << endl;
    }

}



double Graph::many_path_cost(vector< vector<Vertex> > &path_vectors){

	double total_cost = 0.0;
	for(int i = 0; i < path_vectors.size(); i++){
		total_cost = total_cost + path_cost(path_vectors[i]);
	}
	return total_cost;
}



inline double Graph::heuristic_function(Vertex &current_vertex_id, Vertex &neighbour_id, Vertex &target_id){
    return manhatan_distance(neighbour_id, target_id);

}



vector<Vertex> Graph::find_shortest_path_a_star(int &source_id, int &target_id){

    for(int i = 0; i < m_n_vertexes; i++){
        m_distances[i] = std::numeric_limits<double>::max();
        m_previous[i] = -1;
        m_visited[i] = false;
    }
    
    m_distances[source_id] = 0.0;

    priority_queue< pair<int, double>, vector< pair<int, double> >, PriorityQueueCompareLess> pq;
    pq.push(pair<int, double>(source_id, 0.0));

    while( !pq.empty() ){
        pair<int, double> current_vertex = pq.top();
        pq.pop();

        int current_vertex_id = current_vertex.first;

        if(current_vertex_id == target_id)
			break;

        m_visited[current_vertex_id] = true;

        for(int i = 0; i < m_adjacency_list[current_vertex_id].size(); ++i){

            int neighbour_vertex_id = m_adjacency_list[current_vertex_id][i].first;

            if(m_visited[neighbour_vertex_id] == true)
                continue;

            double d = m_distances[current_vertex_id] + m_adjacency_list[current_vertex_id][i].second;
            if(d < m_distances[neighbour_vertex_id]){
                m_distances[neighbour_vertex_id] = d;
                m_previous[neighbour_vertex_id] = current_vertex_id;

                // there is a problem with TerrainCrossingVis.java if you visit the same item/target more than once
                // this if statement avoids going so
                if((neighbour_vertex_id != target_id) && (m_vertex_array[neighbour_vertex_id].get_vertex_type() != 0 ))
                	continue;
                else
                	pq.push(pair<int, double>(neighbour_vertex_id, d ));
                	// effectivelly we are using the Djikstra algorithm since we are not using the heuristic function
                	//pq.push(pair<int, double>(neighbour_vertex_id, d + heuristic_function(m_vertex_array[neighbour_vertex_id], m_vertex_array[target_id]) ));
            }
        }
    }

    vector<Vertex> shortest_path;
    int u = target_id;
    while(m_previous[u] != -1){
        shortest_path.push_back(m_vertex_array[u]);
        u = m_previous[u];
    }
    shortest_path.push_back(m_vertex_array[u]);

    return shortest_path;
}




vector<Vertex> Graph::find_shortest_path_a_star(int &&source_id, int &&target_id){

    for(int i = 0; i < m_n_vertexes; i++){
        m_distances[i] = std::numeric_limits<double>::max();
        m_previous[i] = -1;
        m_visited[i] = false;
    }
    m_distances[source_id] = 0.0;

    priority_queue< pair<int, double>, vector< pair<int, double> >, PriorityQueueCompareLess> pq;
    pq.push(pair<int, double>(source_id, 0.0));

    while( !pq.empty() ){
        pair<int, double> current_vertex = pq.top();
        pq.pop();

        int current_vertex_id = current_vertex.first;

        if(current_vertex_id == target_id)
			break;

        m_visited[current_vertex_id] = true;

        for(int i = 0; i < m_adjacency_list[current_vertex_id].size(); ++i){

            int neighbour_vertex_id = m_adjacency_list[current_vertex_id][i].first;

            if(m_visited[neighbour_vertex_id] == true)
                continue;

            double d = m_distances[current_vertex_id] + m_adjacency_list[current_vertex_id][i].second;
            if(d < m_distances[neighbour_vertex_id]){
                m_distances[neighbour_vertex_id] = d;
                m_previous[neighbour_vertex_id] = current_vertex_id;

                // there is a problem with TerrainCrossingVis.java if you visit the same item/target more than once
                // this if statement avoids going so
                if((neighbour_vertex_id != target_id) && (m_vertex_array[neighbour_vertex_id].get_vertex_type() != 0 ))
                	continue;
                else
                	pq.push(pair<int, double>(neighbour_vertex_id, d ));
                	// effectivelly we are using the Djikstra algorithm since we are not using the heuristic function
                	// finding an effective heurestic function was harder than I enticipated,
                	//non of the standard ones worked (were visably faster than standard Dijkstra)
                	//pq.push(pair<int, double>(neighbour_vertex_id, d + heuristic_function(m_vertex_array[neighbour_vertex_id], m_vertex_array[target_id]) ));
            }
        }
    }

    vector<Vertex> shortest_path;
    int u = target_id;
    while(m_previous[u] != -1){
        shortest_path.push_back(m_vertex_array[u]);
        u = m_previous[u];
    }
    shortest_path.push_back(m_vertex_array[u]);

    return shortest_path;
}



double Graph::path_cost(vector<Vertex> &path){

    double cost = 0.0;
    for(int i = 1; i < path.size(); i++){

        int current_id = path[i].get_id();
        int previous_id = path[i-1].get_id();
        for(int j = 0; j < m_adjacency_list[current_id].size(); j++){
            if (m_adjacency_list[current_id][j].first == previous_id)
                cost = cost + m_adjacency_list[current_id][j].second;
        }
    }

    return cost;
}


double Graph::path_cost(vector<Vertex> &&path){

    double cost = 0.0;
    for(int i = 1; i < path.size(); i++){

        int current_id = path[i].get_id();
        int previous_id = path[i-1].get_id();
        for(int j = 0; j < m_adjacency_list[current_id].size(); j++){
            if (m_adjacency_list[current_id][j].first == previous_id)
                cost = cost + m_adjacency_list[current_id][j].second;
        }
    }

    return cost;
}



//---------------------------------
//------TerrainCrossing class------
//---------------------------------


class TerrainCrossing {
private:

    int m_number_of_items;
    int m_number_of_border_vertexes;
    int m_map_size;
    int m_number_of_items_and_targets;
	int m_max_capacity;

    vector< vector<int> > m_map_matrix;

    vector< tuple<int, double, double, bool> > m_items_vector;
    vector< tuple<int, double, double, bool> > m_targets_vector;

    vector<int> m_border_vertex_ids;

public:

    TerrainCrossing();

    vector< vector<int> > generate_random_map_matrix(int N);

    void initialize_map_matrix(vector<string> &input_map);
    void initialize_item_and_target_vectors(vector<double> &locations);
    void initialize_border_vertexes(Graph &g);

    // Greedy approach methods

    int find_outermost_item(int metric_type);

    int get_closest_vertex_euclidean(int &current_item_id);
    pair<int, double> get_closest_item_a_star(Graph &g, int &source_id, int metric_type);
    pair<int, double> get_closest_target_a_star(Graph &g, int &source_id, int metric_type);
    int get_closest_border_vertex(Graph &g, int &source_id, int metric_type);
    vector<double> get_final_path(vector< vector<Vertex> > &path_vectors);
    vector<Vertex> get_final_path_as_vertex_vector(vector< vector<Vertex> >  &path_vectors);

    vector<vector<Vertex> > search_for_path_the_random_greedy_way(Graph &g, int &capacity, int &smetric_type);


    void print_path(vector<Vertex> &path);
    void print_final_path(Graph &g, vector< vector<Vertex> > &path_vectors);

    void optimize_final_vertex_path_random(vector<Vertex> &path_vectors, int n_iterations);
    void optimize_final_vertex_path_exact(vector<Vertex> &path_vectors, int grid_size);
    vector<double> convert_final_vertex_path_to_normal(vector<Vertex> &path_vectors);

    // Stochastic approach methods

    vector<StochasticVertex> construct_initial_state(Graph &g);
    int make_swap(Graph &g, vector<StochasticVertex> &state, int &greater, int &smaller);
    pair<int, int> get_smaller_and_greater();
    double get_cost_of_stochastic_vertex_vector(vector<StochasticVertex> &state);
	double recalculate_previous_distance_in_stochastic_vertex_vector(Graph &g, vector<StochasticVertex> &state);
    vector<StochasticVertex> metropolis_alg(Graph &g, vector<StochasticVertex> &state, int n_iterations);

    vector<double> getPath(vector<string> input_map, vector<double> locations, int capacity);
};



//------------------------------------------------
//------TerrainCrossing class implementation------
//------------------------------------------------



TerrainCrossing::TerrainCrossing(){

    m_number_of_items = 0;
    m_number_of_border_vertexes = 0;
    m_map_size = 0;
	
	m_number_of_items_and_targets = 0;
	m_max_capacity = 0;

    m_map_matrix = vector< vector<int> >();
    m_items_vector = vector< tuple<int, double, double, bool> >();
    m_targets_vector = vector< tuple<int, double, double, bool> >();

    m_border_vertex_ids = vector<int>();

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
				m_map_matrix[i][j] = world_map[j][i] - '0';
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



void TerrainCrossing::initialize_item_and_target_vectors(vector<double> &locations){

	m_items_vector.resize(m_number_of_items);
	m_targets_vector.resize(m_number_of_items);

	int index = 0;
	for(int i = 0; i < m_number_of_items; i++){
		get<0>(m_items_vector[i]) = index;
		get<1>(m_items_vector[i]) = locations[2*i];;
		get<2>(m_items_vector[i]) = locations[2*i + 1];
		get<3>(m_items_vector[i]) = false;
		index++;
	}

	for(int i = 0; i < m_number_of_items; i++){
		get<0>(m_targets_vector[i]) = index;
		get<1>(m_targets_vector[i]) = locations[2*i + 2*m_number_of_items];;
		get<2>(m_targets_vector[i]) = locations[2*i + 1 + 2*m_number_of_items];
		get<3>(m_targets_vector[i]) = false;
		index++;
	}
}



void TerrainCrossing::initialize_border_vertexes(Graph &g){

    m_number_of_border_vertexes = 4*m_map_size - 4;

    vector< vector<int> > border_vertexes;
    border_vertexes.resize(m_number_of_border_vertexes);
    for(int i = 0; i < border_vertexes.size(); i++)
        border_vertexes[i].resize(2);

    int index = 0;
    int last_row = m_map_size-1;
    for(int i = 0; i < m_map_size; i++){
		//first row 00 01 02 03 04 ...
        border_vertexes[index][0] = 0;
        border_vertexes[index][1] = i;
        index++;

		// last row 40 41 42 43 ... if map 5x5
        border_vertexes[index][0] = last_row;
        border_vertexes[index][1] = i;
        index++;
    }

    for(int i = 1; i < m_map_size-1; i++){
		// first column 10 20 30 40 ...
        border_vertexes[index][0] = i;
        border_vertexes[index][1] = 0;
        index++;

		// last column 14 24 34 44 ...
        border_vertexes[index][0] = i;
        border_vertexes[index][1] = last_row;
        index++;
    }

    m_border_vertex_ids.resize(m_number_of_border_vertexes);
    for(int k = 0; k < m_number_of_border_vertexes; k++){
        m_border_vertex_ids[k] = g.return_vertex_id(border_vertexes[k][0], border_vertexes[k][1]);
    }

}



int TerrainCrossing::find_outermost_item(int metric_type = 0){

	double x_y_center = (double)m_map_size/2.0;

	double outermost_item = -1;
	double d = std::numeric_limits<double>::max();
	for(int i = 0; i < m_number_of_items; i++){
		double xt = get<1>(m_items_vector[i]);
		double yt = get<2>(m_items_vector[i]);

		double current_d = sqrt( pow(x_y_center - xt, 2) + pow(x_y_center - xt, 2) );
		if(current_d < d){
			d = current_d;
			outermost_item = i;
		}
	}
	return outermost_item;
}



int TerrainCrossing::get_closest_vertex_euclidean(int &current_item_id){

	double xi = get<1>(m_items_vector[current_item_id]);
	double yi = get<2>(m_items_vector[current_item_id]);

	double d = std::numeric_limits<double>::max();
	int closest_target_id = -1;


	for(int i = 0; i < m_number_of_items; i++){

		if(get<3>(m_targets_vector[i]) == true)
			continue;

		double xt = get<1>(m_targets_vector[i]);
		double yt = get<2>(m_targets_vector[i]);

		double current_d = sqrt( pow(xi - xt, 2) + pow(yi - xt, 2) );
		if (current_d < d){
			d = current_d;
			closest_target_id = i;
		}
	}

	return closest_target_id;
}



pair<int, double> TerrainCrossing::get_closest_item_a_star(Graph &g, int &source_id, int metric_type = 0){

    double d = std::numeric_limits<double>::max();
	int closest_item_id = -1;


	for(int i = 0; i < m_number_of_items; i++){
		if(get<3>(m_items_vector[i]) == true)
			continue;

		double current_d = std::numeric_limits<double>::max();
		if(metric_type == 0){
			vector<Vertex> path = g.find_shortest_path_a_star(source_id, get<0>(m_items_vector[i]));
			current_d = g.path_cost(path);
		} else {

			double xi = g.return_vertex_x(source_id);
			double yi = g.return_vertex_y(source_id);

			double xt = get<1>(m_items_vector[i]);
			double yt = get<2>(m_items_vector[i]);

			if (metric_type == 1)
				current_d = sqrt( pow(xi - xt, 2) + pow(yi - yt, 2) );
			if (metric_type == 2)
				current_d = abs(xi - xt) + abs(yi - yt);
		}

		if (current_d < d){
			d = current_d;
			closest_item_id = get<0>(m_items_vector[i]);
		}
	}

	return pair<int, double>(closest_item_id, d);

}



pair<int, double> TerrainCrossing::get_closest_target_a_star(Graph &g, int &source_id, int metric_type = 0){

    double d = std::numeric_limits<double>::max();
	int closest_target_id = -1;


	for(int i = 0; i < m_number_of_items; i++){
		if(get<3>(m_targets_vector[i]) == true)
			continue;

		double current_d = std::numeric_limits<double>::max();
		if(metric_type == 0){
			vector<Vertex> path = g.find_shortest_path_a_star(source_id, get<0>(m_targets_vector[i]));
			current_d = g.path_cost(path);
		} else {

			double xi = g.return_vertex_x(source_id);
			double yi = g.return_vertex_y(source_id);

			double xt = get<1>(m_targets_vector[i]);
			double yt = get<2>(m_targets_vector[i]);

			if (metric_type == 1)
				current_d = sqrt( pow(xi - xt, 2) + pow(yi - yt, 2) );
			if (metric_type == 2)
				current_d = abs(xi - xt) + abs(yi - yt);
		}


		if (current_d < d){
			d = current_d;
			closest_target_id = get<0>(m_targets_vector[i]);
		}
	}

	return pair<int, double>(closest_target_id, d);

}



int TerrainCrossing::get_closest_border_vertex(Graph &g, int &source_id, int metric_type = 0){

	double d = std::numeric_limits<double>::max();
	int closest_border_vertex_id = -1;

	for(int i = 0; i < m_border_vertex_ids.size(); i++){

		vector<Vertex> path = g.find_shortest_path_a_star(source_id, m_border_vertex_ids[i]);
		double current_d = g.path_cost(path);
		if (current_d < d){
			d = current_d;
			closest_border_vertex_id = m_border_vertex_ids[i];
		}

	}
	return closest_border_vertex_id;
}



vector<double> TerrainCrossing::get_final_path(vector< vector<Vertex> > &path_vectors){


	vector<double> final_path;
	int n_paths = path_vectors.size();

	// add the entry point
	Vertex start_vertex = path_vectors[0][path_vectors[0].size()-1];
	int si = start_vertex.get_i();
	int sj = start_vertex.get_j();

	double sx = start_vertex.get_x();
	double sy = start_vertex.get_y();

	if(si == 0)
                sx = EPSILON;

	if((sj == 0) && (si > 0) && (si < m_map_size-1))
                sy = EPSILON;

	if((sj == m_map_size-1) && (si > 0) && (si < m_map_size-1))
                sy = sy + 0.5 - EPSILON;

	if(si == m_map_size-1)
                sx = m_map_size - EPSILON;

	final_path.push_back(sx);
        final_path.push_back(sy);


	for(int i = 0; i < path_vectors.size(); i++){
		for(int j = path_vectors[i].size()-2; j >= 0; j--){
            // the small offset is required to meet the 1E-3 condition of
            // visiting an item/target
			double so = 0.0005; // small offset

			double px = path_vectors[i][j].get_x();
			double py = path_vectors[i][j].get_y();

			if(path_vectors[i][j].get_vertex_type() == 0){
				final_path.push_back(px);
				final_path.push_back(py);
				continue;
			}


			int px_i = floor(px);
			int py_j = floor(py);

			double cx = px_i + 0.5;
			double cy = py_j + 0.5;

			if( px > cx && py > cy ){
				px = px - so;
				py = py - so;
				final_path.push_back(px);
				final_path.push_back(py);
				continue;
			}

			if( px < cx && py < cy ){
				px = px + so;
				py = py + so;
				final_path.push_back(px);
				final_path.push_back(py);
				continue;
			}

			if( px < cx && py > cy ){
				px = px + so;
				py = py - so;
				final_path.push_back(px);
				final_path.push_back(py);
				continue;
			}

			if( px > cx && py < cy ){
				px = px - so;
				py = py + so;
				final_path.push_back(px);
				final_path.push_back(py);
				continue;
			}


		}
	}

	// add the departure point from the map
	Vertex final_vertex = path_vectors[path_vectors.size()-1][0];
	int fi = final_vertex.get_i();
	int fj = final_vertex.get_j();

	double fx = final_vertex.get_x();
	double fy = final_vertex.get_y();

    // these ifs are mutually exclussive
	if(fi == 0)
		fx = EPSILON;

	if((fj == 0) && (fi > 0) && (fi < m_map_size-1))
		fy = EPSILON;

	if((fj == m_map_size-1) && (fi > 0) && (fi < m_map_size-1))
		fy = fy + 0.5 - EPSILON;

	if(fi == m_map_size-1)
		fx = m_map_size - EPSILON;

	final_path.push_back(fx);
	final_path.push_back(fy);

	return final_path;
}




vector<Vertex> TerrainCrossing::get_final_path_as_vertex_vector(vector< vector<Vertex> >  &path_vectors){

	// this is similar as get_final_path but returns a vector of Vertexes

	vector<Vertex> final_path;
	int n_paths = path_vectors.size();

	// add the entry point
	Vertex start_vertex = path_vectors[0][path_vectors[0].size()-1];
	int si = start_vertex.get_i();
	int sj = start_vertex.get_j();

	double sx = start_vertex.get_x();
	double sy = start_vertex.get_y();

	int stt = start_vertex.get_tt();
	int s_vertex_type = start_vertex.get_vertex_type();
	int s_id = start_vertex.get_id();

	if(si == 0)
                sx = EPSILON;

	if((sj == 0) && (si > 0) && (si < m_map_size-1))
                sy = EPSILON;

	if((sj == m_map_size-1) && (si > 0) && (si < m_map_size-1))
                sy = sy + 0.5 - EPSILON;

	if(si == m_map_size-1)
                sx = m_map_size - EPSILON;


	final_path.push_back(Vertex(si, sj, sx, sy, stt, s_vertex_type, s_id));


	for(int i = 0; i < path_vectors.size(); i++){
		for(int j = path_vectors[i].size()-2; j >= 0; j--){
            // the small offset is required to meet the 1E-3 condition of
            // visiting an item/target
			double so = 0.0005; // small offset

			int px_i = path_vectors[i][j].get_i();
			int py_j = path_vectors[i][j].get_j();

			double px = path_vectors[i][j].get_x();
			double py = path_vectors[i][j].get_y();

			int ptt = path_vectors[i][j].get_tt();
			int p_vertex_type = path_vectors[i][j].get_vertex_type();
			int p_id = path_vectors[i][j].get_id();


			if(path_vectors[i][j].get_vertex_type() == 0){
				final_path.push_back(Vertex(px_i, py_j, px, py, ptt, p_vertex_type, p_id));
				continue;
			}


			double cx = px_i + 0.5;
			double cy = py_j + 0.5;

			if( px > cx && py > cy ){
				px = px - so;
				py = py - so;
				final_path.push_back(Vertex(px_i, py_j, px, py, ptt, p_vertex_type, p_id));
				continue;
			}

			if( px < cx && py < cy ){
				px = px + so;
				py = py + so;
				final_path.push_back(Vertex(px_i, py_j, px, py, ptt, p_vertex_type, p_id));
				continue;
			}

			if( px < cx && py > cy ){
				px = px + so;
				py = py - so;
				final_path.push_back(Vertex(px_i, py_j, px, py, ptt, p_vertex_type, p_id));
				continue;
			}

			if( px > cx && py < cy ){
				px = px - so;
				py = py + so;
				final_path.push_back(Vertex(px_i, py_j, px, py, ptt, p_vertex_type, p_id));
				continue;
			}


		}
	}
	
	// add the departure point from the map
	Vertex final_vertex = path_vectors[path_vectors.size()-1][0];
	int fi = final_vertex.get_i();
	int fj = final_vertex.get_j();

	double fx = final_vertex.get_x();
	double fy = final_vertex.get_y();

	int ftt = final_vertex.get_tt();
	int f_vertex_type = final_vertex.get_vertex_type();
	int f_id = final_vertex.get_id();

        // these ifs are mutually exclussive
	if(fi == 0)
		fx = EPSILON;

	if((fj == 0) && (fi > 0) && (fi < m_map_size-1))
		fy = EPSILON;

	if((fj == m_map_size-1) && (fi > 0) && (fi < m_map_size-1))
		fy = fy + 0.5 - EPSILON;

	if(fi == m_map_size-1)
		fx = m_map_size - EPSILON;

	final_path.push_back(Vertex(fi, fj, fx, fy, ftt, f_vertex_type, f_id));


	return final_path;

}



vector< vector<Vertex> > TerrainCrossing::search_for_path_the_random_greedy_way(Graph &g, int &capacity, int &metric_type){

        std::random_device rd;
        std::mt19937 engine(rd());
        std::uniform_int_distribution<> border_vertex_dist(0, m_number_of_border_vertexes - 1);

        vector< vector<Vertex> > path_vectors;

        int number_of_unvisited_targets = m_number_of_items; // number_of_targets == number_of_items
        int number_of_unvisited_items = m_number_of_items;
        int number_of_items_in_bag = 0;
        int rbv = m_border_vertex_ids[border_vertex_dist(engine)]; // random border vertex, this is the starting vertex of the path

        pair<int, double> start_item = get_closest_item_a_star(g, rbv, metric_type);
        int starting_point_id = start_item.first;
        double start_item_dist = start_item.second;

        vector<Vertex> first_path = g.find_shortest_path_a_star(rbv, starting_point_id);
        path_vectors.push_back(first_path);

        get<3>(m_items_vector[starting_point_id]) = true;
        number_of_items_in_bag++;
        number_of_unvisited_items--;

        while(number_of_unvisited_targets != 0){

            pair<int, double> closest_item;
            int closest_item_id = -1;
            double closest_item_dist = std::numeric_limits<double>::max();

            if((number_of_items_in_bag == 0) && (number_of_unvisited_items != 0)){

                pair<int, double> closest_item = get_closest_item_a_star(g, starting_point_id, metric_type);
                closest_item_id = closest_item.first;
                closest_item_dist = closest_item.second;

                vector<Vertex> middle_path = g.find_shortest_path_a_star(starting_point_id, closest_item_id);
                path_vectors.push_back(middle_path);

                get<3>(m_items_vector[closest_item_id]) = true;
                starting_point_id = closest_item_id;

                number_of_items_in_bag++;
                number_of_unvisited_items--;
                continue;
            }

            if((number_of_items_in_bag < capacity) && (number_of_unvisited_items != 0)){
                pair<int, double> closest_item = get_closest_item_a_star(g, starting_point_id, metric_type);
                closest_item_id = closest_item.first;
                closest_item_dist = closest_item.second;
            }

            pair<int, double> closest_target;
            int closest_target_id = -1;
            double closest_target_dist = std::numeric_limits<double>::max();

            if(number_of_items_in_bag != 0){
                pair<int, double> closest_target = get_closest_target_a_star(g, starting_point_id, metric_type);
                closest_target_id = closest_target.first;
                closest_target_dist = closest_target.second;
            }

            if(closest_item_dist < closest_target_dist){
                vector<Vertex> middle_path = g.find_shortest_path_a_star(starting_point_id, closest_item_id);
                path_vectors.push_back(middle_path);

                get<3>(m_items_vector[closest_item_id]) = true;
                starting_point_id = closest_item_id;

                number_of_items_in_bag++;
                number_of_unvisited_items--;
                continue;
            } else {
                vector<Vertex> middle_path = g.find_shortest_path_a_star(starting_point_id, closest_target_id);
                path_vectors.push_back(middle_path);

                get<3>(m_targets_vector[closest_target_id - m_number_of_items]) = true;
                starting_point_id = closest_target_id;

                number_of_unvisited_targets--;
                number_of_items_in_bag--;
                continue;
            }
        }

        // path found, reset if items were visited
        for(int i = 0; i < m_number_of_items; i++){
        	get<3>(m_items_vector[i]) = false;
        	get<3>(m_targets_vector[i]) = false;
        }

        // find final vertex of the path
        int final_vertex_id = get_closest_border_vertex(g, starting_point_id);
        vector<Vertex> final_path = g.find_shortest_path_a_star(starting_point_id, final_vertex_id);
        path_vectors.push_back(final_path);

        return path_vectors;

}



void TerrainCrossing::print_path(vector<Vertex> &path){

	for(int i = 0; i < path.size(); i++)
		cerr << path[i] << endl;
	cerr << endl;
}



void TerrainCrossing::optimize_final_vertex_path_random(vector<Vertex> &path_vectors, int n_iterations){

	std::random_device rd;
	std::mt19937 engine(rd());
	//std::uniform_real_distribution<> real_dist(-0.495, 0.495);

    std::vector<double> intervals = {-0.495,  -0.28, 0.28, 0.495};
    std::vector<double> weights = {  1, 0,  1};
    std::piecewise_constant_distribution<> real_dist(intervals.begin(), intervals.end(), weights.begin());

	for(int i = 1; i < path_vectors.size()-1; i++){

		int vertex_type = path_vectors[i].get_vertex_type();
		if(vertex_type != 0)
			continue;

		double initial_d = triple_path_distance(path_vectors[i-1], path_vectors[i], path_vectors[i+1]);

		double px = path_vectors[i].get_x();
		double py = path_vectors[i].get_y();

		double optimal_shift_x = 0.0;
		double optimal_shift_y = 0.0;

		for(int j = 0; j < n_iterations; j++){

			double shift_x = real_dist(engine);
			double shift_y = real_dist(engine);

			path_vectors[i].set_x(px + shift_x);
			path_vectors[i].set_y(py + shift_y);

			double modified_d = triple_path_distance(path_vectors[i-1], path_vectors[i], path_vectors[i+1]);

			if(modified_d < initial_d){
				
				double next_x = path_vectors[i+1].get_x();
				double next_y = path_vectors[i+1].get_y();

				double prev_x = path_vectors[i-1].get_x();
				double prev_y = path_vectors[i-1].get_y();
				
				double dn =  sqrt( pow( px + shift_x - next_x ,2) + pow( py + shift_y - next_y ,2) );
				double dp =  sqrt( pow( px + shift_x - prev_x ,2) + pow( py + shift_y - prev_y ,2) );				
				
				
				optimal_shift_x = shift_x;
				optimal_shift_y = shift_y;
				
			}
		}

		path_vectors[i].set_x(px + optimal_shift_x);
		path_vectors[i].set_y(py + optimal_shift_y);
	}
}




void TerrainCrossing::optimize_final_vertex_path_exact(vector<Vertex> &path_vectors, int grid_size){

	// optimize_final_vertex_path_random provides better results on average

	vector<double> grid = linspace(0.05, 0.95, grid_size);
	for(int i = 1; i < path_vectors.size()-1; i++){

		int vertex_type = path_vectors[i].get_vertex_type();
		if(vertex_type != 0)
			continue;

		double initial_d = triple_path_distance(path_vectors[i-1], path_vectors[i], path_vectors[i+1]);

		double px_i = (double)path_vectors[i].get_i();
		double py_j = (double)path_vectors[i].get_j();

		double px = path_vectors[i].get_x();
		double py = path_vectors[i].get_y();

		double optimal_shift_x = 0.0;
		double optimal_shift_y = 0.0;

		bool shift_found = false;
		for(int in = 0; in < grid_size; in++){
			for(int jn = 0; jn < grid_size; jn++){

				path_vectors[i].set_x(px_i + grid[in]);
				path_vectors[i].set_y(py_j + grid[jn]);

				double modified_d = triple_path_distance(path_vectors[i-1], path_vectors[i], path_vectors[i+1]);

				if(modified_d < initial_d){
					optimal_shift_x = grid[in];
					optimal_shift_y = grid[jn];
					shift_found = true;
				}
			}
		}

		if(shift_found == true){
			path_vectors[i].set_x(px_i + optimal_shift_x);
			path_vectors[i].set_y(py_j + optimal_shift_y);
		} else {
			path_vectors[i].set_x(px);
			path_vectors[i].set_y(py);
		}

	}
}




vector<double> TerrainCrossing::convert_final_vertex_path_to_normal(vector<Vertex> &path_vectors){

	vector<double> ret(2*path_vectors.size());

	int index = 0;
	for(int i = 0; i < path_vectors.size(); i++){
		ret[index] = path_vectors[i].get_x();
		index++;
		ret[index] = path_vectors[i].get_y();
		index++;
	}
	return ret;
}



vector<StochasticVertex> TerrainCrossing::construct_initial_state(Graph &g){


    vector<StochasticVertex> initial_state_vector(2*m_number_of_items, StochasticVertex());

    int index = 0;
    for(int i = 0; i < m_number_of_items; i++){
        initial_state_vector[index] = StochasticVertex(get<0>(m_items_vector[i]));
        initial_state_vector[index].set_capacity(1);
        initial_state_vector[index].set_vertex_type(1);
        index++;
        initial_state_vector[index] = StochasticVertex(get<0>(m_targets_vector[i]));
        initial_state_vector[index].set_capacity(0);
        initial_state_vector[index].set_vertex_type(2);
        index++;
    }

    
    initial_state_vector[0].set_next(&initial_state_vector[1]);
    

    for(int i = 1; i < initial_state_vector.size()-1; i++){
    	initial_state_vector[i].set_next(&initial_state_vector[i+1]);
    	initial_state_vector[i].set_prev(&initial_state_vector[i-1]);
    
    	int source = initial_state_vector[i].get_id();
    	int target = initial_state_vector[i-1].get_id();
    	
    	vector<Vertex> path = g.find_shortest_path_a_star(source, target);
    	double d = g.path_cost(path);
		
		//double d = euclidean_distance( g.get_vertex_in_graph( source ), g.get_vertex_in_graph( target ) );
    	
    	initial_state_vector[i].set_distance_to_prev(d);
    }
    
    initial_state_vector[initial_state_vector.size()-1].set_prev(&initial_state_vector[initial_state_vector.size()-2]);
       
    
	int source = initial_state_vector[initial_state_vector.size()-1].get_id();
	int target = initial_state_vector[initial_state_vector.size()-2].get_id();
    
    initial_state_vector[initial_state_vector.size()-1].set_distance_to_prev( g.path_cost( g.find_shortest_path_a_star(source, target) ) );
    
    return initial_state_vector;
}




int TerrainCrossing::make_swap(Graph &g, vector<StochasticVertex> &state, int &greater, int &smaller){
	
	
	int last = m_number_of_items_and_targets-1;
	if((smaller > 0) && (greater < last)){
	

		if( state[greater].get_vertex_type() == 2 && (state[smaller].get_prev())->get_capacity() == 0 )
			return 0;

		if( state[greater].get_vertex_type() == 1 && (state[smaller].get_prev())->get_capacity() == m_max_capacity)
			return 0;
		
		
		StochasticVertex* gn = state[greater].get_next();
		StochasticVertex* sp = state[smaller].get_prev();

		double d1 = g.path_cost( g.find_shortest_path_a_star(state[smaller].get_id(), gn->get_id() ) );
		double d2 = g.path_cost( g.find_shortest_path_a_star(state[greater].get_id(), sp->get_id() ) );
		
		//double d1 = euclidean_distance( g.get_vertex_in_graph( state[smaller].get_id() ), g.get_vertex_in_graph( gn->get_id() ) );
		//double d2 = euclidean_distance( g.get_vertex_in_graph( state[greater].get_id() ), g.get_vertex_in_graph( sp->get_id() ) );
		
		gn->set_distance_to_prev(d1);
		state[smaller].set_distance_to_prev(d2);
		
		
		int gtype = state[greater].get_vertex_type();
		int stype = state[smaller].get_vertex_type();
		
		if(stype == 1 && gtype == 2){
			int x = state[smaller].get_capacity();
			
			state[smaller].set_capacity(x-2);
			state[greater].set_capacity(x-1);
			
			state[smaller].set_vertex_type(gtype);
			state[greater].set_vertex_type(stype);
			
		} else if (stype == 2 && gtype == 1) {
			int x = state[greater].get_capacity();
			
			state[smaller].set_capacity(x+1);
			state[greater].set_capacity(x);
			
			state[smaller].set_vertex_type(gtype);
			state[greater].set_vertex_type(stype);
		}

		
		int temp_id = state[greater].get_id();
		state[greater].set_id( state[smaller].get_id() );
		state[smaller].set_id( temp_id );
		
	}

	
	if(smaller == 0){
		
		StochasticVertex* nv = state[smaller].get_next();
		if(nv->get_vertex_type() == 2)
			return 0;
		
		
		int temp_id = state[smaller].get_id();
		state[smaller].set_id( nv->get_id() );
		nv->set_id( temp_id );
		
		double d = g.path_cost( g.find_shortest_path_a_star(nv->get_id(), (nv->get_next())->get_id() ) );
		//double d = euclidean_distance( g.get_vertex_in_graph( nv->get_id() ), g.get_vertex_in_graph( (nv->get_next())->get_id() ) );
		
		(nv->get_next())->set_distance_to_prev(d);
		
	}

	

	if(greater == last){
		
		StochasticVertex* pv = state[greater].get_prev();
		if(pv->get_vertex_type() == 1)
			return 0;
		
		int temp_id = state[greater].get_id();
		state[greater].set_id( pv->get_id() );
		pv->set_id( temp_id );
		
		
		double d = g.path_cost( g.find_shortest_path_a_star( pv->get_id(), (pv->get_prev())->get_id()) );
		//double d = euclidean_distance( g.get_vertex_in_graph( pv->get_id() ), g.get_vertex_in_graph( (pv->get_prev())->get_id() ) );

		pv->set_distance_to_prev(d);
		
	}
	
	
	return 1;
}



pair<int, int> TerrainCrossing::get_smaller_and_greater(){
	
	std::random_device rd;
	std::mt19937 engine(rd());
	std::uniform_int_distribution<> dist(0, m_number_of_items_and_targets-1);
	
	int swap_position = dist(engine);

	std::uniform_int_distribution<> left_right(0, 1);
	int lr = left_right(engine);
	
	if(lr == 0)
		lr = -1;
	
	int next_to_swap_position = swap_position + lr;
	
	int greater = 0;
	int smaller = 0;
	
	if( swap_position > next_to_swap_position){
		greater = swap_position;
		smaller = next_to_swap_position;
	} else {
		smaller = swap_position;
		greater = next_to_swap_position;
	}
	
	return pair<int, int>(greater, smaller);
}



double TerrainCrossing::recalculate_previous_distance_in_stochastic_vertex_vector(Graph &g, vector<StochasticVertex> &state){
	
	for(int i = 1; i < state.size(); i++){
		double d =  g.path_cost( g.find_shortest_path_a_star( state[i].get_id(), state[i-1].get_id()) );
		state[i].set_distance_to_prev(d);
	}
	
}



double TerrainCrossing::get_cost_of_stochastic_vertex_vector(vector<StochasticVertex> &state){

	double cost = 0.0;
	for(int i = 0; i < state.size(); i++)
		cost = cost + state[i].get_distance_to_prev();
	return cost;
	
}



vector<StochasticVertex> TerrainCrossing::metropolis_alg(Graph &g, vector<StochasticVertex> &state, int n_iterations){
	
	std::random_device rd;
	std::mt19937 engine(rd());
	uniform_real_distribution<double> uni(0.0, 1.0);
	
	double min_energy = get_cost_of_stochastic_vertex_vector(state);
	
	vector<StochasticVertex> min_state = state;
	
	bool recalc = false;
	double E1 = min_energy;
	for(int i = 0; i < n_iterations; i++){
		
		if(recalc == true)
			E1 = get_cost_of_stochastic_vertex_vector(state);
		
		
		pair<int, int> gs_pair = get_smaller_and_greater();
		
		int swap = make_swap(g, state, gs_pair.first, gs_pair.second);
		
		if (swap == 0)
			continue;
		
		double E2 = get_cost_of_stochastic_vertex_vector(state);
		
		double T = 1.0;
		double A = exp((E1-E2)/T);
				
		double p = uni(engine);
		
		if(p < A){
			if( E2 < min_energy ){
				min_state = state;
				recalc = true;
				continue;
			}
		} else {
			make_swap(g, state, gs_pair.first, gs_pair.second);
			recalc = false;
		}
	}
	return min_state;
	
}



vector<double> TerrainCrossing::getPath(vector<string> world_map, vector<double> locations, int capacity) {

        std::chrono::time_point<std::chrono::system_clock> start, end;

        m_number_of_items = locations.size()/4;
        m_number_of_items_and_targets = 2*m_number_of_items;
        
        m_map_size = world_map.size();
		m_max_capacity = capacity;
        int offset = 5;

        initialize_map_matrix(world_map);

        Graph g = Graph(m_map_matrix, locations);

        initialize_item_and_target_vectors(locations);
        initialize_border_vertexes(g);
		
        start = std::chrono::system_clock::now();

        int metric_type = 2;


        vector<double> ret;
        vector< vector<Vertex> > optimal_path_vectors;
        double min_cost = std::numeric_limits<double>::max();
		
		int n_iter = 0;
		
		if(m_map_size <= 10)
        	n_iter = 1000;
        else if (m_map_size <= 20)
        	n_iter = 500;
        else if (m_map_size <= 30)
        	n_iter = 250;
        else if (m_map_size <= 40)
        	n_iter = 100;
        else
        	n_iter = 50;
		
		
        for(int i = 0; i < n_iter; i++){

        	vector< vector<Vertex> > path_vectors = search_for_path_the_random_greedy_way(g, capacity, metric_type);
			double current_cost =  g.many_path_cost(path_vectors);
			if(current_cost < min_cost){
				min_cost = current_cost;
				optimal_path_vectors = path_vectors;

			}
        }


        vector<Vertex> final_path_vertex_vector = get_final_path_as_vertex_vector(optimal_path_vectors);

        if(m_map_size <= 10)
        	optimize_final_vertex_path_random(final_path_vertex_vector, 200000);
        else if (m_map_size <= 20)
        	optimize_final_vertex_path_random(final_path_vertex_vector, 10000);
        else if (m_map_size <= 30)
        	optimize_final_vertex_path_random(final_path_vertex_vector, 2500);
        else if (m_map_size <= 40)
        	optimize_final_vertex_path_random(final_path_vertex_vector, 2500);
        else
        	optimize_final_vertex_path_random(final_path_vertex_vector, 2500);


        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds_in_main;
        elapsed_seconds_in_main = end-start;
        //cerr << "get closest item elapsed time: " << elapsed_seconds_in_main.count() << endl;


        
        //check is any two points are less than 0.001 apart
        // if yes then shift them appart so they are at least 0.001 appart
		for(int i = 1; i < final_path_vertex_vector.size(); i++){
		
			double next_x = final_path_vertex_vector[i].get_x();
			double next_y = final_path_vertex_vector[i].get_y();

			double prev_x = final_path_vertex_vector[i-1].get_x();
			double prev_y = final_path_vertex_vector[i-1].get_y();
				
			double d =  sqrt( pow( next_x - prev_x ,2) + pow( next_y - prev_y ,2) );
			
			if(d < TWO_TIMES_EPSILON){
				
				while(d < TWO_TIMES_EPSILON){
					
					std::random_device rd;
					std::mt19937 engine(rd());
					uniform_real_distribution<double> dist(-0.1, 0.1);
					
					double xx = dist(engine);
					double yy = dist(engine);
					
					double xi = final_path_vertex_vector[i].get_i();
					double yj = final_path_vertex_vector[i].get_j();
					
					double bup = xi + TWO_TIMES_EPSILON;
					double bdown = xi + 1 - TWO_TIMES_EPSILON;
					
					double bleft = yj + TWO_TIMES_EPSILON;
					double bright = yj + 1 - TWO_TIMES_EPSILON;
					
					if(final_path_vertex_vector[i].get_vertex_type() == 0){
						
						double xs = next_x + xx;
						double ys = next_y + yy;
						if( xs < bup || xs > bdown || ys < bleft || ys > bright)
							continue;
						
						final_path_vertex_vector[i].set_x( xs );
						final_path_vertex_vector[i].set_y( ys );
					}
					
					if(final_path_vertex_vector[i-1].get_vertex_type() == 0){
						
						double xs = prev_x + xx;
						double ys = prev_y + yy;
						if( xs < bup || xs > bdown || ys < bleft || ys > bright)
							continue;
						
						final_path_vertex_vector[i-1].set_x( prev_x + xx );
						final_path_vertex_vector[i-1].set_y( prev_y + yy );
					}
					
					d =  sqrt( pow( final_path_vertex_vector[i].get_x() - final_path_vertex_vector[i-1].get_x() ,2) + pow( final_path_vertex_vector[i].get_y() - final_path_vertex_vector[i-1].get_y() ,2) );					
				}
			}
				
		}
		
		ret = convert_final_vertex_path_to_normal(final_path_vertex_vector);
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






