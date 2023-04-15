#include <igl/opengl/glfw/Viewer.h>
#ifndef HALFEDGE_DS_HEADER
#define HALFEDGE_DS_HEADER
#include "HalfedgeDS.cpp"
#endif

#include <queue>
#include "partitioning.h"


using namespace Eigen;
using namespace std;

MatrixXi color_region (const MatrixXi& R, int region, const vector<vector<int>>& anchors, const MatrixXd& V, HalfedgeDS& he);
pair<MatrixXi,MatrixXi> triangulation (MatrixXi& Regions, const vector<vector<int>>& Anchors, const MatrixXd& Vertices, const MatrixXi& Faces, HalfedgeDS& he);