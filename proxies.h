#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;
using namespace std;

MatrixXd new_proxies_L_2(const MatrixXi& Regions, const MatrixXi& Faces, const MatrixXd& Vertices, int proxy_num);
MatrixXd new_proxies_L_2_1(const MatrixXi& Regions, const MatrixXi& Faces, const MatrixXd& Vertices, int proxy_num);
MatrixXd new_proxies(const MatrixXi& Regions, const MatrixXi& Faces, const MatrixXd& V, int k, int norme);