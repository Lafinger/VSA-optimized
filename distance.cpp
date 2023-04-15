#include "distance.h"

//******************Distortion error******************

VectorXd Face_area;
MatrixXd Face_normal;
MatrixXd Face_center;
MatrixXi Face;


Vector3d get_center(int i) {
  return Face_center.row(i);
}

Vector3d get_normal(int i) {
  return Face_normal.row(i);
}

double get_area(int i) {
  return Face_area(i);
}
double triangle_area(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3){

  double l1 = (v2-v1).norm();
  double l2 = (v3-v2).norm();
  double l3 = (v1-v3).norm();
  double p = (l1+l2+l3)/2.;
  double S = pow(p*(p-l1)*(p-l2)*(p-l3),0.5);

  return S;
}

double orthogonal_distance(Vector3d X, Vector3d N, Vector3d M){
  
  double n = N.norm();
  return fabs((N).dot(M-X))/n;

}

Vector3d triangle_normal(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3){

  Vector3d N = (v2-v1).cross(v3-v1);
  double n = N.norm();

  return N/n;

}

// optimized
Vector3d triangle_center(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3){
  return (v1 + v2 + v3) / 3.0;
}

Vector3d triangle_normal(Vector3i T, MatrixXd V){
  return triangle_normal(V.row(T(0)),V.row(T(1)),V.row(T(2)));
}


Vector3d triangle_center(Vector3i T, MatrixXd V){
  return (V.row(T(0)) + V.row(T(1)) + V.row(T(2))) / 3.0;
}

double distance_L_2(int i, Vector3d X, Vector3d N, MatrixXd V){

  Vector3d v1 = V.row(Face(i,0));
  Vector3d v2 = V.row(Face(i,1));
  Vector3d v3 = V.row(Face(i,2));

  double area = triangle_area(v1,v2,v3);
  double d1 = orthogonal_distance(X,N,v1);
  double d2 = orthogonal_distance(X,N,v2);
  double d3 = orthogonal_distance(X,N,v3);

  return area*(d1*d1 + d2*d2 + d3*d3 + d1*d2 + d1*d3 + d2*d3) / 6.0;

}

double distance_L_2_1(int face_index, const Vector3d& proxy_normal) {
  Vector3d face_normal = Face_normal.row(face_index);
  // Vector3d nf = Vector3d(Face_normal(face_index,0),Face_normal(face_index,1),Face_normal(face_index,2));
  double norm = (face_normal - proxy_normal).norm();
  return Face_area(face_index) * norm * norm;
}

double distance(int face_index, const Vector3d& proxy_barycenter, const Vector3d& proxy_normal, const MatrixXd& Vertices, int norme){
  if (norme == 0){
    return distance_L_2(face_index, proxy_barycenter, proxy_normal, Vertices);
  }
  else {
    return distance_L_2_1(face_index, proxy_normal);
  }
}

double global_distortion_error(const MatrixXi& Regions, const MatrixXd& Proxies, const MatrixXd& Vertices, const MatrixXi& Faces, int norme){

  int proxy_num = Proxies.rows()/2; // number of proxies
  int face_num = Faces.rows(); // number of faces
  double E = 0.; // global error

  // sums the distortion errors of each triangle with its proxy
  int proxy_index; // index of the region/proxy
  Vector3d proxy_barycenter; // barycenter of the proxy
  Vector3d proxy_normal; // normal of the proxy
  Vector3i T; // triangle
  double e; // error of the triangle T with its proxy X,N
  
  for (int face_index = 0; face_index < face_num; ++face_index){

    proxy_index = Regions(face_index,0);
    proxy_barycenter = Proxies.row(proxy_index);
    proxy_normal = Proxies.row(proxy_index + proxy_num);
    e = distance(face_index, proxy_barycenter, proxy_normal, Vertices, norme);
    E += e;

  }

  return E;

}



double distance_projection(const MatrixXd& V, const MatrixXd& Proxies, const int anchor1, const int anchor2, const int v, const int r1, const int r2){
  int p = Proxies.rows()/2;
  // sin angle
  Vector3d p1 = Proxies.row(p+r1);
  Vector3d p2 = Proxies.row(p+r2);
  double ang = (p1.cross(p2)).norm(); // proxies are normalized

  // distance to segment
  Vector3d x = V.row(v)-V.row(anchor1);
  Vector3d s = V.row(anchor2)-V.row(anchor1); //segment
  if (s.norm()==0) return 0.;
  double t = x.dot(s)/s.norm();

  return (x-t*s/s.norm()).norm()/s.norm()*ang;
}


void initialize_normals_areas(const MatrixXi& F, const MatrixXd& V) {
  int face_num = F.rows();
  Face = F;
  Face_center.setZero(face_num,3);
  Face_area.setZero(face_num);
  Face_normal.setZero(face_num,3);
  for (int i=0 ; i<face_num ; i++){
    Vector3d v1 = V.row(F(i,0));
    Vector3d v2 = V.row(F(i,1));
    Vector3d v3 = V.row(F(i,2));
    Face_area(i) = triangle_area(v1,v2,v3);
    Face_normal.row(i) = triangle_normal(v1,v2,v3);
    // Face_center.row(i) = triangle_center(F.row(i),V);
    Face_center.row(i) = triangle_center(v1,v2,v3); // optimized
  }
}