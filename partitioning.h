#include <igl/opengl/glfw/Viewer.h>
#include <queue>
#include <map>

using namespace Eigen;
using namespace std;

MatrixXi face_adjacency(const MatrixXi& Faces, int n);
void tcolor(MatrixXi &R);
void fcolor(MatrixXd &Cf, MatrixXi Ad);
void distance_color(MatrixXd &Cf, MatrixXi F, MatrixXd V, int norme);
void initial_partition(int proxy_num, MatrixXi &Regions, const MatrixXd& Vertices, const MatrixXi& Faces, const MatrixXi& FaceAdjacencies, int norme);
void initial_partition2(int p, MatrixXi &R, MatrixXd V, MatrixXi F, MatrixXi Ad, int norme);
VectorXi find_best_triangles(const MatrixXi& Regions, const MatrixXd& Proxies, const MatrixXd& Vertices, const MatrixXi& Faces, int norme);
void proxy_fitting(MatrixXi &Regions, const MatrixXd& Proxies, const MatrixXd& Vertices, const MatrixXi& Faces, const MatrixXi& FaceAdjacencies, int norme);




