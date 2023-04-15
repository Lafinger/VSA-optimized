#include "proxies.h"
#include "distance.h"

Vector3d g(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3){

  return (v1+v2+v3)/3.0;

};

MatrixXd M(Vector3d v1,Vector3d v2,Vector3d v3){

  MatrixXd M = MatrixXd::Zero(3,3);
  M.row(0) = v2-v1;
  M.row(1) = v3-v1;

  return M;

}; 

Vector3d new_Xi_L_2 (const MatrixXi& Regions, int proxy_index, const MatrixXi& Faces, const MatrixXd& Vertices) {

  Vector3d proxy_barycenter(0.0, 0.0, 0.0);
  double w = 0.0;

  Vector3i vertex_index_in_face;
  Vector3d v1;
  Vector3d v2;
  Vector3d v3;
  Vector3d face_barycenter;
  double face_area;

  for (int face_index = 0; face_index < Regions.rows(); ++face_index){
    //we only add the triangles that belong to the region i
    if (Regions(face_index, 0) == proxy_index){

      vertex_index_in_face = Faces.row(face_index);
      v1 = Vertices.row(vertex_index_in_face(0));
      v2 = Vertices.row(vertex_index_in_face(1));
      v3 = Vertices.row(vertex_index_in_face(2));

      face_barycenter = g(v1,v2,v3);
      face_area = get_area(face_index);

      proxy_barycenter += face_area*face_barycenter;
      w += face_area;
    }
  }
  
  return proxy_barycenter/w;

};

Vector3d new_Ni_L_2 (MatrixXi R, int i, MatrixXi F, MatrixXd V){

  //Compute the Covariance Matrix Ci
  MatrixXd Ci = MatrixXd::Zero(3,3); 
  double w = 0.;
  MatrixXd Xi = new_Xi_L_2(R,i,F,V);

  MatrixXd A(3,3);
  A(0,0) = 10; A(0,1) = 7; A(0,2) = 0;
  A(1,0) = 7; A(1,1) = 10; A(1,2) = 0;
  A(2,0) = 0; A(2,1) = 0.; A(2,2) = 0;

  Vector3i T;
  Vector3d v1;
  Vector3d v2;
  Vector3d v3;
  Vector3d gT;
  MatrixXd MT;
  double s;

  for (int f=0 ; f<R.rows() ; f++){
    //we only add the triangles that belong to the region i
    if (R(f,0) == i){

      T = F.row(f);
      v1 = V.row(T(0));
      v2 = V.row(T(1));
      v3 = V.row(T(2));

      gT = g(v1,v2,v3);
      MT = M(v1,v2,v3);
      s = get_area(f);
      
      Ci += (2./72.)*s*MT*A*MT.transpose() + s*gT*gT.transpose();
      w += s;
    } 
  }
  
  Ci = Ci - w*Xi*Xi.transpose();
  
  //Find the eigenvector of the min eigenvalue of Ci
  EigenSolver<MatrixXd> es(Ci);

  MatrixXd valp;
  MatrixXd vectp;
  valp = es.eigenvalues().real();
  vectp = es.eigenvectors().real();

  double min_valp;
  MatrixXd::Index minRow, minCol;
  Vector3d Ni;
  min_valp = valp.minCoeff(&minRow,&minCol);
  Ni = vectp.col(minRow);

  return Ni.normalized();

};

Vector3d new_Xi_L_2_1 (const MatrixXi& Regions, int proxy_index, const MatrixXi& Faces, const MatrixXd& Vertices) {

  return new_Xi_L_2(Regions, proxy_index, Faces, Vertices);

};


MatrixXd compute_N (const MatrixXi& Regions, int proxy_num) {
  MatrixXd Normals;
  Normals.setZero(proxy_num, 3);

  double face_area;
  Vector3d face_normal;

  for (int face_index = 0; face_index < Regions.rows(); ++face_index){
    //we only add the triangles that belong to the region i
    int proxy_index = Regions(face_index,0);

    face_area = get_area(face_index);
    face_normal = get_normal(face_index); 

    Normals.row(proxy_index) += face_area*face_normal;
  }

  return Normals;
}


Vector3d new_Ni_L_2_1 (MatrixXi R, int i, MatrixXi F, MatrixXd V){

  Vector3d Ni(0.,0.,0.);

  Vector3i T;
  Vector3d v1;
  Vector3d v2;
  Vector3d v3;
  double s;
  Vector3d nT;

  for (int f=0 ; f<R.rows() ; f++){
    //we only add the triangles that belong to the region i
    if (R(f,0) == i){
      s = get_area(f);
      nT = get_normal(f); 
      Ni += s*nT;
    }
  }

  return Ni.normalized();

};

//proxy_num is the number of regions/proxies of the partition R
MatrixXd new_proxies_L_2(const MatrixXi& Regions, const MatrixXi& Faces, const MatrixXd& Vertices, int proxy_num) {

  MatrixXd P(2*proxy_num,3);

  Vector3d Xi;
  Vector3d Ni;

  for (int proxy_index = 0; proxy_index < proxy_num; ++proxy_index){
    Xi = new_Xi_L_2(Regions, proxy_index, Faces, Vertices);
    Ni = new_Ni_L_2(Regions, proxy_index, Faces, Vertices);
    P.row(proxy_index) = Xi;
    P.row(proxy_num + proxy_index) = Ni.normalized();
  }

  return P;

};

MatrixXd new_proxies_L_2_1(const MatrixXi& Regions, const MatrixXi& Faces, const MatrixXd& Vertices, int proxy_num) {

  MatrixXd P(2*proxy_num,3);

  Vector3d barycenter;
  MatrixXd Normals = compute_N(Regions, proxy_num);

  for (int proxy_index = 0; proxy_index < proxy_num; ++proxy_index){
    barycenter = new_Xi_L_2_1(Regions, proxy_index, Faces, Vertices);
    // barycenter of each proxy
    P.row(proxy_index) = barycenter;
    // normal of each proxy
    P.row(proxy_num + proxy_index) = Normals.row(proxy_index).normalized();
  }

  return P;

};

MatrixXd new_proxies(const MatrixXi& Regions, const MatrixXi& Faces, const MatrixXd& Vertices, int proxy_num, int norme) {
  if (norme == 0){
    return new_proxies_L_2(Regions,Faces,Vertices,proxy_num);
  }
  else if (norme == 1){
    return new_proxies_L_2_1(Regions,Faces,Vertices,proxy_num);
  }
  else {
    cout<<"wrong norme parameter"<<endl;
  }
};
