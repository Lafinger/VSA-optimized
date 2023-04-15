#include <igl/opengl/glfw/Viewer.h>
#include "anchors.h"

using namespace Eigen;
using namespace std;

map<int,int> renumber(MatrixXi& F);
MatrixXd new_V(HalfedgeDS& he, const MatrixXd& V, const MatrixXd& Proxies, const MatrixXi& R, map<int, int>& index);
Vector3d projection(const Vector3d& z, const MatrixXd& Proxies, int k);