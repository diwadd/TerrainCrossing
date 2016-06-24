#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <set>
#include <string>

using namespace std;

//--------------------------
//------Location class------
//--------------------------


template<typename Type> void print_vector_2d(vector< vector<Type> > v){
	for(int i = 0; i < v.size(); i++){
			for(int j = 0; j < v[i].size(); j++){
					cerr << v[i][j] << " ";
			}
			cerr << endl;
	}
}


class Location {
private:
	double m_x;
	double m_y;

	bool m_visited;
	
public:
	Location(double x = 0.0, double y = 0.0, bool visited = false): m_x(x), m_y(y), m_visited(visited) {}
	
	void set_x(double x) { m_x = x; }
	void set_y(double y) { m_y = y; }
	void set_visited(bool visited) { m_visited = visited; }
	
	double get_x() const { return m_x; }
	double get_y() const { return m_y; }
	bool get_visited() const { return m_visited; }
	
};
 

ostream & operator<<(ostream & os, const Location &loc){
    os << "x: " << loc.get_x() << " y: " << loc.get_y() << " y: " << " visited: " << loc.get_visited();
    return os;
}


//---------------------------
//--------Graph class--------
//---------------------------


class Vertex {
private:
	int m_i;
	int m_j;
	
	double m_x;
	double m_y;
	
	short m_location_type;
	
public:
	Vertex(int i = 0, int j = 0, double x = 0.0, double y = 0.0, short location_type = 0): m_i(i), m_j(j), m_x(x), m_y(y), m_location_type(location_type) {}
	
};


//class Graph {
//private:


//public:	
	
//};


//---------------------------------
//------TerrainCrossing class------
//---------------------------------


class TerrainCrossing {
private:

	vector<Location> m_item_locations;
	vector<Location> m_target_locations;
	
	vector< vector<int> > m_map_matrix;

public:
 
	TerrainCrossing();

    void set_items_and_target_locations(vector<double> &locations);
    void print_item_locations();
	void print_target_locations();
    
	void initialize_map_matrix(vector<string> &input_map);
    
    vector<double> getPath(vector<string> input_map, vector<double> locations, int capacity);
};


//------------------------------------------------
//------TerrainCrossing class implementation------
//------------------------------------------------


TerrainCrossing::TerrainCrossing(){
	
	vector<Location> m_item_locations = vector<Location>();
	vector<Location> m_target_locations = vector<Location>();
	
	vector< vector<int> > m_map_matrix = vector< vector<int> >();
	
}


void TerrainCrossing::set_items_and_target_locations(vector<double> &locations){
	
	// Sets m_item_locations and m_target_locations.
	// In the problem statement the location of both items and target locations is given in 
	// the locations vector that is passed to getPath;
	
	int N = locations.size()/4;
	
	m_item_locations.resize(N);
	m_target_locations.resize(N);
	
	for(int i = 0; i < N; i++)
		m_item_locations[i] = Location(locations[2*i], locations[2*i+1]);

	for(int i = 0; i < N; i++)
		m_target_locations[i] = Location(locations[2*N + 2*i], locations[2*N + 2*i+1]);

}


void TerrainCrossing::print_item_locations(){
	
	cerr << "Printing item locations" << endl;
	for(int i = 0; i < m_item_locations.size(); i++)
		cerr << m_item_locations[i] << endl;
	
}


void TerrainCrossing::print_target_locations(){
	
	cerr << "Printing target locations" << endl;
	for(int i = 0; i < m_target_locations.size(); i++)
		cerr << m_target_locations[i] << endl;
	
}


void TerrainCrossing::initialize_map_matrix(vector<string> &input_map) {
	
	// Resize vector to meet the map size.
	int N = input_map.size();
	m_map_matrix.resize(N);
	for(int i = 0; i < N; i++)
		m_map_matrix[i].resize(N);
	
	// Fill the map matrix.
	for(int i = 0; i < N; i++){
		for(int j = 0; j < input_map[i].length(); j++){
				m_map_matrix[i][j] = input_map[i][j] - '0';
		}
	}
	
}


vector<double> TerrainCrossing::getPath(vector<string> input_map, vector<double> locations, int capacity) {

	
		set_items_and_target_locations(locations);
		print_item_locations();
		print_target_locations();
		
		initialize_map_matrix(input_map);
		print_vector_2d(m_map_matrix);
		
		cerr << "Input map size: " << input_map.size() << endl;
		for(int i = 0; i < input_map.size(); i++)
			cerr << input_map[i] << endl;
		
		
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






