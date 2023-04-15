#include <igl/opengl/glfw/Viewer.h>
#ifndef HALFEDGE_DS_HEADER
#define HALFEDGE_DS_HEADER
#include "HalfedgeDS.cpp"
#endif

#include "distance.h"


using namespace Eigen;
using namespace std;

bool vector_contains(const vector<int>& v, int key);
bool vector_contains(const vector<double>& v, double key);

vector<vector<int>> anchor_points(HalfedgeDS& he, const MatrixXi& Regions, const MatrixXd& V, const MatrixXd& Proxies, double treshold);
vector<int> find_vertex_proxies(HalfedgeDS& he, int vertex_index, const MatrixXi& Regions);