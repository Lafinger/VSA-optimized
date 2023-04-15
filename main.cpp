#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/jet.h>
#include <iostream>
#include <ostream>
#include <chrono>

#include "HalfedgeBuilder.cpp"

#include "partitioning.h"
#include "distance.h"
#include "proxies.h"
#include "anchors.h"
#include "triangulation.h"
#include "renumbering.h"

// M_PI error
#define M_PI 3.14159265358979323846

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd Vertices; // matrix storing vertex coordinates of the input mesh (n rows, 3 columns)
MatrixXi Faces; // incidence relations between faces and edges (f columns)
MatrixXi Regions; // index of proxy each face 
MatrixXd newV; // matrix storing vertex coordinates of the input mesh (n rows, 3 columns)
MatrixXi newF; // incidence relations between faces and edges (f columns)
MatrixXi newR; // matrix indicating the partition of each vertex
MatrixXd Colorings; // the coloring
MatrixXd Proxies;
MatrixXi FaceAdjacencies; // face adjacency
int proxy_num; // number of proxies
int norme; //0 = norm L_2, 1 = norm L_2_1
HalfedgeDS* he;
int iterations;
vector<pair<int,double>> global_error_points; // contains the global_distortion_error according to the number of iterations
double error;
double precedent_error;
double treshold;

void debug_regions_vides(MatrixXi Regions, int p){
  cout<<"Regions vides"<<endl;
  bool trouve_j;
  for (int j=0 ; j<p ; j++){
    trouve_j = false;
    for (int i=0 ; i<Regions.rows() ; i++){
      if (Regions(i,0)==j){
        trouve_j = true;
      }
    }
    if (trouve_j == false){
      cout<<j<<endl;
    }
  }
  cout<<"fin"<<endl;
};

MatrixXd covariance(MatrixXd M) {
  Vector3d mean = M.colwise().mean();
  for (int j=0 ; j<M.rows(); j++) { 
    M.row(j) -= mean;
  }
  return M.transpose()*M;
}

pair<Vector3d,Vector3d> compute_ellipse_vectors(int c){
  vector<Vector3d> Radius;
  for(int i = 0; i < Faces.rows(); i++) {
    if (Regions(i,0)==c) {
      for (int j=0;j<3;j++) {
        Vector3d q = Vertices.row(Faces(i,j));
        Vector3d q2 = Proxies.row(Regions(i,0));
        Radius.push_back(q-q2);
      }
    }
  }
  int k = Radius.size();
  MatrixXd M;
  M.setZero(k,3);
  for (int j=0 ; j<k; j++) { 
    M.row(j) = Radius[j];
  }  
  EigenSolver<MatrixXd> eig(covariance(M)/k);
  MatrixXd ev = eig.eigenvectors().real();
  MatrixXd eva = eig.eigenvalues().real();
  
  Vector3d e1,e2;
  for (int l=0; l<3;l++) {
    if (eva(l) == eva.maxCoeff()) {
      e1 = ev.col(l)*pow(eva(l),0.5);
      eva(l) = -eva(l);
      break;
    }
  }
  for (int l=0; l<3;l++) {
    if (eva(l) == eva.maxCoeff()) {
      e2 = ev.col(l)*pow(eva(l),0.5);
    }
  }
  Vector3d n = Proxies.row(proxy_num+c);
  if ( n.dot(e1.cross(e2))<0) {
    return make_pair(-e1,e2);
  }
  return make_pair(e1,e2);
}

void draw_tangent(igl::opengl::glfw::Viewer &viewer) {
   for (int i =0; i<proxy_num;i++) {
    viewer.append_mesh();
    viewer.data(0).add_points(Proxies.row(i), Eigen::RowVector3d(1, 0, 0));
    viewer.data(0).add_edges(
        Proxies.row(i),
        Proxies.row(i) + Proxies.row(i+proxy_num)/10.0,
        Eigen::RowVector3d(1, 0, 0));
  }
}

void color_scheme(igl::opengl::glfw::Viewer &viewer, const MatrixXd& Vertices, const MatrixXi& Faces) {
  viewer.data().clear();
  int f = Faces.rows();
  MatrixXd nC(f,1);
  for (int i=0; i<f; i++){
    Vector3d c =(Vertices.row(Faces(i,0)) + Vertices.row(Faces(i,1)) + Vertices.row(Faces(i,2))) / 30.0;
    nC(i,0)=c(1);
  }
  igl::jet(nC,true,Colorings);
  viewer.data().set_mesh(Vertices, Faces);
  viewer.data().set_colors(Colorings);
}

void draw_anchors(igl::opengl::glfw::Viewer &viewer) {
  vector<vector<int>> anchors = anchor_points(*he, Regions, Vertices, Proxies,treshold);
  for(size_t i = 0; i < anchors.size(); i++) {
    for(size_t j = 0; j < anchors[i].size(); j++) {
      viewer.data(0).add_points(Vertices.row(anchors[i][j]), Eigen::RowVector3d(1,1,0));
    }
  }
    
}

void triangle_proxy(Vector3d x, Vector3d n, MatrixXd& newV, int k, Vector3d m1, Vector3d m2) {

  int M=20;
  newV.row(M*k) = x;
  for (int i=1; i<M; i++) {
    double t = i*2*M_PI/(M-1);
    newV.row(M*k+i) = x + sin(t)*m1 + cos(t)*m2;
      // result.row(i) <<1,2,3;
  }
}

void draw_prox(igl::opengl::glfw::Viewer &viewer) {
  int M=20;

  MatrixXd newV;
  newV.setZero(M*proxy_num,3);
  MatrixXi newF;
  MatrixXi newR0;
  newR0.setZero((M-1)*proxy_num,3);
  newF.setZero((M-1)*proxy_num,3);

  for(int i = 0; i < proxy_num; i++) {
    pair<Vector3d,Vector3d> vec = compute_ellipse_vectors(i);  
    triangle_proxy(Proxies.row(i),Proxies.row(i+proxy_num), newV, i, vec.first, vec.second);
    for (int j=1; j<M-1; j++) {
      newF.row((M-1)*i+j-1) << M*i+j,M*i,M*i+j+1;
      newR0((M-1)*i+j-1)=i;
    }
    newF.row((M-1)*(i+1)-1) << M*i+M-1,M*i,M*i+1;
    newR0((M-1)*(i+1)-1)=i;
  }
  viewer.data().clear();
  viewer.data().set_mesh(newV, newF);
  igl::jet(newR0,true,Colorings);
  viewer.data(0).set_colors(Colorings);

}

void one_iter(igl::opengl::glfw::Viewer &viewer) {
  proxy_fitting(Regions, Proxies, Vertices,  Faces, FaceAdjacencies, norme);
  Proxies = new_proxies(Regions, Faces, Vertices, proxy_num, norme);
  iterations += 1;
  error = global_distortion_error(Regions, Proxies, Vertices, Faces, norme);
  cout<<"Global Error : "<<error<<endl;
  global_error_points.push_back(make_pair(iterations,error));
  precedent_error = error;
  igl::jet(Regions,true,Colorings);
  viewer.data(0).set_colors(Colorings);
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
  cout << "pressed Key: " << key << " " << (unsigned int)key << endl;
  if (key=='1') {
    viewer.data().clear();
  }
  if (key=='2') {
    color_scheme(viewer, Vertices, Faces);
  }
  if (key=='3') {
    // debug_regions_vides(Regions,p);
    one_iter(viewer);
  }
  if (key=='4') {
    draw_anchors(viewer);
  }
  if (key=='5') {
    chrono::time_point<chrono::high_resolution_clock> begin = chrono::high_resolution_clock::now();

    treshold = 0.4;
    vector<vector<int>> Anchors = anchor_points(*he, Regions, Vertices, Proxies,treshold);

     //MatrixXi Cr = color_region(Regions, 6, Anchors,Vertices,*he);
     // viewer.append_mesh();
     // for(int j = 0; j < Cr.rows(); j++) {
     //   int i = Cr(j,0);
     //   if (i>-1) viewer.data(0).add_points(Vertices.row(j), Eigen::RowVector3d(i%3/2.0,i/9.0, i%2));
     // }
     // return true;

     pair<MatrixXi,MatrixXi> new_F_and_R = triangulation(Regions, Anchors, Vertices, Faces, *he);
     newF = new_F_and_R.first;
     newR = new_F_and_R.second;

     map<int,int> index = renumber(newF); //modifies Faces
     newV = new_V(*he, Vertices, Proxies, Regions, index);
     viewer.data().clear();
     // igl::jet(newR,true,Colorings);
     viewer.data().set_mesh(newV, newF);
     // viewer.data().set_colors(Colorings);
     color_scheme(viewer, newV, newF);
     cout <<"faces : "<<newF.rows() << endl;

    cout << "triangulation time : " << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - begin).count() << endl;
  }
  if (key=='6') {
    color_scheme(viewer, newV, newF);
  }
  if (key=='7') {
    draw_prox(viewer);
  }
  if (key=='8') {
    for (int i=0;i<10;i++) one_iter(viewer);
    cout << "    Done" <<endl;
  }
  if (key=='9') {
    for (int i=0;i<100;i++) one_iter(viewer);
    cout << "    Done" <<endl;
  }
  if (key == 'S' || (unsigned int)key == 83) {
    chrono::time_point<chrono::high_resolution_clock> begin = chrono::high_resolution_clock::now();

    vector<double> errors;
    while (fabs(error - precedent_error) > 0.9){
      precedent_error = error;
      proxy_fitting(Regions, Proxies, Vertices,  Faces, FaceAdjacencies, norme);
      Proxies = new_proxies(Regions, Faces, Vertices, proxy_num, norme);
      iterations += 1;
      error = global_distortion_error(Regions, Proxies, Vertices, Faces, norme);
      cout << error << endl;

      if (vector_contains(errors,error)){
        cout<<"cycle !"<<endl;
        break;
      }
     
      // error = global_distortion_error(Regions, Proxies, Vertices, Faces, norme);
      global_error_points.push_back(make_pair(iterations,error));
      errors.push_back(error);

      igl::jet(Regions,true,Colorings);
      viewer.data(0).set_colors(Colorings);
    }

    cout << "partition time : " << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - begin).count() << endl;
  }
  if (key == 'A' || (unsigned int)key == 65 ) {

    if(proxy_num >= 20) {
      viewer.data().clear();
      color_scheme(viewer, Vertices, Faces);
      // time point
      chrono::time_point<chrono::high_resolution_clock> begin = chrono::high_resolution_clock::now();

      // partition initialization
      proxy_num -= 10;
      initial_partition(proxy_num, Regions, Vertices, Faces, FaceAdjacencies, norme);
      Proxies = new_proxies(Regions, Faces, Vertices, proxy_num, norme);
      iterations = 1;
      error = global_distortion_error(Regions, Proxies, Vertices, Faces, norme); // E(TX,PY) = ∑(0~x)∑(0~y) E(tx,py))
      precedent_error = error - 1;                              // what?
      global_error_points.clear();
      global_error_points.push_back(make_pair(iterations, error));

      igl::jet(Regions, true, Colorings);
      viewer.data(0).set_colors(Colorings);

      cout << "partition time : " << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - begin).count() << endl;
    }

  }
  if (key == 'D' || (unsigned int)key == 68 ) {

    if (proxy_num <= Faces.rows() - 10) {
      viewer.data().clear();
      color_scheme(viewer, Vertices, Faces);
      chrono::time_point<chrono::high_resolution_clock> begin = chrono::high_resolution_clock::now();

      // partition initialization
      proxy_num += 10;
      initial_partition(proxy_num, Regions, Vertices, Faces, FaceAdjacencies, norme);
      Proxies = new_proxies(Regions, Faces, Vertices, proxy_num, norme);
      iterations = 1;
      error = global_distortion_error(Regions, Proxies, Vertices, Faces, norme); // E(TX,PY) = ∑(0~x)∑(0~y) E(tx,py))
      precedent_error = error - 1;                              // what?
      global_error_points.clear();
      global_error_points.push_back(make_pair(iterations, error));

      igl::jet(Regions, true, Colorings);
      viewer.data(0).set_colors(Colorings);

      cout << "partition time : " << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - begin).count() << endl;
    }
  }

  return false;
}

// ------------ main program ----------------
int main(int argc, char *argv[])
{
  chrono::time_point<chrono::high_resolution_clock> total_time = chrono::high_resolution_clock::now();

  chrono::time_point<chrono::high_resolution_clock> file_read_time = chrono::high_resolution_clock::now();
  string file = "../data/dragon.off";
  igl::readOFF(file, Vertices, Faces); // Load an input mesh in OFF format
  //  print the number of mesh elements
  cout << "Vertices: " << Vertices.rows() << endl;
  cout << "Faces:    " << Faces.rows() << endl;
  cout << "file read time : " << (chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - file_read_time).count())/1000.0 << " sec " << endl;

  chrono::time_point<chrono::high_resolution_clock> halfedge_build_time = chrono::high_resolution_clock::now();
  HalfedgeBuilder* builder = new HalfedgeBuilder();  
  HalfedgeDS he2 = builder->createMesh(Vertices.rows(), Faces); 
  he = &he2;
  cout << "halfedge build time : " << (chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - halfedge_build_time).count())/1000.0 << " sec " << endl;

  // Face adjacency
  chrono::time_point<chrono::high_resolution_clock> face_adjacency_build_time = chrono::high_resolution_clock::now();
  FaceAdjacencies = face_adjacency(Faces, Vertices.rows());
  cout << "face adjacency build time : " << (chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - face_adjacency_build_time).count())/1000.0 << " sec " << endl;

  // calculating normals,center and area of face of mesh
  chrono::time_point<chrono::high_resolution_clock> face_attribute_build_time = chrono::high_resolution_clock::now();
  initialize_normals_areas(Faces, Vertices);
  cout << "total attributes of face build time : " << (chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - face_attribute_build_time).count())/1000.0 << " sec " << endl;

  //coloring 
  // Partition_faces.setZero(Faces.rows(),1);
  // MatrixXd C;
  // tcolor(Partition_faces);
  // igl::jet(Partition_faces,true,C);

  // coloring adjacency 
  // MatrixXd Cf;
  // fcolor(Cf,FaceAdjacencies);
  // MatrixXd C;
  // igl::jet(Cf,true,C);

  // coloring distance 
  // MatrixXd Cf;
  // distance_color(Cf,Faces,Vertices,0);
  // MatrixXd C;
  // igl::jet(Cf,true,C);

	/*
	  p(int) : proxy number
	  Regions(MatrixXi&) : matrix indicating the partition of each vertex
	  Vertices(MatrixXd) : mesh vertices
	  Faces(MatrixXi) : mesh face
	  FaceAdjacencies(MatrixXi) : face adjacency
	  norme(int) : 0 is L_2, 1 is L_2_1
	*/
  chrono::time_point<chrono::high_resolution_clock> initial_partition_time = chrono::high_resolution_clock::now();
  proxy_num = 1000;
  norme = 1;
	initial_partition(proxy_num, Regions, Vertices, Faces, FaceAdjacencies, norme);
  cout << "initial partition time : " << (chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - initial_partition_time).count())/1000.0 << " sec " << endl;

  // Proxies is MatrixXd, its (2*k,3) 0~k is barycenter of each proxy in Proxies and k~2k is normal of each proxy in Proxies.
  Proxies = new_proxies(Regions, Faces, Vertices, proxy_num, norme);
  iterations = 1;
  error = global_distortion_error(Regions, Proxies, Vertices, Faces, norme); // E(TX,PY) = ∑(0~x)∑(0~y) E(tx,py))
  precedent_error = error - 1 ; // what?
  global_error_points.push_back(make_pair(iterations,error));
  igl::jet(Regions,true,Colorings);

  cout << "total time : " << (chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - total_time).count())/1000.0 << " sec " << endl;

  //showing normals
  // viewer.append_mesh();
  // for (int j=0;j<Faces.rows();j++) {
  //   Vector3d center = triangle_center(Faces.row(j),Vertices);
  //   Vector3d norm = triangle_normal(Faces.row(j),Vertices);
  //   viewer.data(0).add_edges(
  //       center.transpose(),
  //       center.transpose()+norm.transpose()/10.0,
  //       Eigen::RowVector3d(1, 0, 0));
  // }

  // launch viewer
  igl::opengl::glfw::Viewer viewer; // create the 3d viewer
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().set_mesh(Vertices, Faces); // load a face-based representation of the input 3d shape
  // viewer.data().set_colors(Colorings);
  color_scheme(viewer, Vertices, Faces);
  viewer.launch(); // run the editor
  
}

