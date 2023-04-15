#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;
using namespace std;

//******************Distortion error******************
double orthogonal_distance(Vector3d X, Vector3d N, Vector3d M);
Vector3d get_center(int i);
Vector3d get_normal(int i);
double get_area(int i);


double distance(int face_index, const Vector3d& proxy_barycenter, const Vector3d& proxy_normal, const MatrixXd& Vertices, int norme);

double global_distortion_error(const MatrixXi& Regions, const MatrixXd& Proxies, const MatrixXd& Vertices, const MatrixXi& Faces, int norme);

double distance_projection(const MatrixXd& V, const MatrixXd& Proxies, const int anchor1, const int anchor2, const int v, const int r1, const int r2);

void initialize_normals_areas(const MatrixXi& F, const MatrixXd& V);