#include "partitioning.h"
#include "distance.h"
#include "proxies.h"
#include "anchors.h"

#define MAXFLOAT 3.402823E38

//******************Partitionning******************
MatrixXi face_adjacency(const MatrixXi& Faces, int n) { // n: number of vertices
  // face adjacency matrix, O(m log(m)) performances
 
  int m = Faces.rows();
  MatrixXi Ad;
  Ad.setZero(m,3);
  map<int,int> edge;
  
  for (int i=0;i<m;i++) {
    int r[3] = {Faces(i,0),Faces(i,1),Faces(i,2)};
    sort(r, r + 3);
    // edge mapping: e -> e.min + n*e.max
    int e[3] = {r[0] + n*r[1],r[0] + n*r[2],r[1] + n*r[2]}; // order edges by vertices
    for (int l=0; l<3; l++){
      if (edge.count(e[l])!=0) { // if edge contains(e[l])
        int j = edge[e[l]];
        Ad(i,l) = j; // faces i and j are neighbourgs

        int r2[3] = {Faces(j,0),Faces(j,1),Faces(j,2)};
        sort(r2, r2 + 3);
        int e2[3] = {r2[0] + n*r2[1],r2[0] + n*r2[2],r2[1] + n*r2[2]};
        int k = (e[l]==e2[1]) + 2*(e[l]==e2[2]); // check corresponding index for j
        Ad(j,k) = i; // faces i and j are neighbourgs
      }
      else {
        edge[e[l]] = i;
      }
    }
  }
  return Ad;
}



vector<int> uniform_proxies(int k, int n) {
  vector<int> Proxies;
  for (int i=0;i<k;i++) {
    Proxies.push_back((i*n)/k);
  }
  return Proxies;
}

/*
    return proxies vertor that it save face(triangle) index.
*/
vector<int> random_proxies(int k, int n) {
  vector<int> Proxies;
  while (Proxies.size()<k){
    int x = rand() % n;  
    if (!vector_contains (Proxies,x)) Proxies.push_back(x);
  }
  return Proxies;
}

void tcolor(MatrixXi &R) {
  // color according to the face number
  int n = R.rows();
  for (int i=0; i<n;i++) {
    R(i,0) = i;
  }
}

void fcolor(MatrixXd &Cf, MatrixXi Ad) {
  // closeness to the original face in terms of raw distance
  Cf.setZero(Ad.rows(), 1);
  double color = 1.0;  
  priority_queue<pair<double, int>> q;
  pair<double, int> face;
  q.push(make_pair(color,0));
  while (q.size()!=0) {
    face = q.top();
    q.pop();
    if (Cf(face.second,0)==0) {
      Cf(face.second,0) = face.first;
      color = face.first * 0.95;
      q.push(make_pair(color,Ad(face.second,0)));
      q.push(make_pair(color,Ad(face.second,1)));
      q.push(make_pair(color,Ad(face.second,2)));
    }
  }
  
}


void distance_color(MatrixXd &Cf, MatrixXi F, MatrixXd V, int norme) {
  // closeness to the original face with the l2 or L21 distance
  Cf.setZero(F.rows(),1);
  Vector3d c = get_center(0);
  Vector3d n = get_normal(0);
  for (int i=0;i<F.rows();i++) {
    Cf(i) = distance(i,c,n,V, norme);
  }
}


int find_triangles_region(const vector<int>& SeedTriangles, MatrixXi &Regions, const MatrixXd Vertices, const MatrixXi Faces, const MatrixXi FaceAdjacencies, int norme) {
  // updates R and find for each region the furthest triangle and its region
  int face_num = Faces.rows();
  Regions = -MatrixXi::Ones(face_num, 1); // face_num rows 1 col, all of values are -1

  // Error and tag
  priority_queue< pair<double, int> > queue_priority;
  
  // number of proxies
  int proxy_num = SeedTriangles.size();
  // don't forget delete[]
  Vector3d* Proxies_center = new Vector3d[proxy_num]; // save barycenter of triangle
  Vector3d* Proxies_normal = new Vector3d[proxy_num]; // save normal of triangle
  //Vector3d Proxies_center[proxy_num];
  //Vector3d Proxies_normal[proxy_num];

  // the furthest triangle and distance of each proxy
  VectorXi furthest_triangle(proxy_num);
  VectorXd furthest_distance(proxy_num);
  furthest_distance.setZero(proxy_num);

  for (int proxy_index = 0; proxy_index < proxy_num; ++proxy_index) {
    // To do : optimize
    Proxies_center[proxy_index] = get_center(SeedTriangles[proxy_index]);
    Proxies_normal[proxy_index] = get_normal(SeedTriangles[proxy_index]);
    // SeedTriangles[proxy_index] is index of F, proxy_index is index of proxies.
    Regions(SeedTriangles[proxy_index]) = proxy_index;
    for (int ad = 0; ad < 3; ++ad) {
        // ad_face_index is index of F indicating a adjacency triangle.
        int ad_face_index = FaceAdjacencies(SeedTriangles[proxy_index],ad);
        // L_2_1 Error metrics between adjacency triangle and proxy
        double error_metric = distance(ad_face_index, Proxies_center[proxy_index], Proxies_normal[proxy_index], Vertices, norme);
        // "-d" because of deflaut priority_queue
        queue_priority.push(make_pair(-error_metric, ad_face_index + face_num * proxy_index));
        if (error_metric>furthest_distance(proxy_index) && !vector_contains(SeedTriangles,ad_face_index)) {

          furthest_distance(proxy_index) = error_metric;
          furthest_triangle(proxy_index) = ad_face_index;
        }
    }
  }

  pair<double, int> item;
  int face_index;
  int proxy_index;

  cout << "[JasonChen] the partition algorithm boostarup first time number of triangles is  " << queue_priority.size() << endl; // input of p = 20, now is 20 * 3 = 60

  // number of cycles = (face number - proxy number) * 3
  while (queue_priority.size()!=0) {
    item = queue_priority.top();
    queue_priority.pop();
    proxy_index = item.second / face_num; // index of seed Triangles
    face_index = item.second % face_num; // adjacency triangle of seed Triangles

    if (Regions(face_index) == -1) {

      Regions(face_index) = proxy_index;
      for (int ad = 0; ad < 3; ++ad) {
        int ad_face_index = FaceAdjacencies(face_index,ad);
        // if(Regions(ad_face_index) == -1) {
          double error_metric = distance(ad_face_index, Proxies_center[proxy_index], Proxies_normal[proxy_index], Vertices, norme);
          queue_priority.push(make_pair(-error_metric, ad_face_index + face_num * proxy_index)); // Number of cycles dependency
          if (error_metric > furthest_distance(proxy_index) && !vector_contains(SeedTriangles,ad_face_index)) {
            furthest_distance(proxy_index) = error_metric;
            furthest_triangle(proxy_index) = ad_face_index;
          }
        // }
      }
      
    }
  }

  double max = furthest_distance(0);
  int maxtri = furthest_triangle(0);
  for (int proxy_index = 1; proxy_index < proxy_num; ++proxy_index) {
    if (max < furthest_distance(proxy_index)) {
      maxtri = furthest_triangle(proxy_index);
      max = furthest_distance(proxy_index);
    } 
  }
  return maxtri;
}

void initial_partition(int proxy_num, MatrixXi &Regions, const MatrixXd& Vertices, const MatrixXi& Faces, const MatrixXi& FaceAdjacencies, int norme) {
  // the Triangles are proxies that it is index of face vector.
  vector<int> SeedTriangles = random_proxies(proxy_num, Faces.rows());
  find_triangles_region(SeedTriangles, Regions, Vertices, Faces, FaceAdjacencies, norme);
}


void initial_partition2(int p, MatrixXi &R, MatrixXd V, MatrixXi F, MatrixXi Ad, int norme) {
  //using furthest triangle
  int tri = rand() % F.rows();
  vector<int> Triangles;
  for (int i=0;i<p;i++) {
      cout << i<<endl;
      Triangles.push_back(tri);
      tri = find_triangles_region(Triangles,R,V,F,Ad,norme);
  }
}


VectorXi find_best_triangles(const MatrixXi& Regions, const MatrixXd& Proxies, const MatrixXd& Vertices, const MatrixXi& Faces, int norme) {
  int proxy_num = Proxies.rows()/2;
  VectorXi FittedTriangles(proxy_num);
  VectorXd distances(proxy_num);
  double normal_distance;
  
  for (int proxy_index = 0; proxy_index < proxy_num; ++proxy_index) distances(proxy_index) = MAXFLOAT;
  
  //look at each face i
  for (int face_index = 0; face_index < Faces.rows(); ++face_index) {
    int proxy_index = Regions(face_index,0);
    normal_distance = distance(face_index, Proxies.row(proxy_index), Proxies.row(proxy_index + proxy_num),  Vertices, norme);
    if (normal_distance < distances(proxy_index)) {
      distances(proxy_index) = normal_distance;
      FittedTriangles(proxy_index) = face_index;
    }
  }
  return FittedTriangles;
}


void proxy_fitting(MatrixXi &Regions, const MatrixXd& Proxies, const MatrixXd& Vertices, const MatrixXi& Faces, const MatrixXi& FaceAdjacencies, int norme) {
  int face_num = Faces.rows();
  // double error = 0; ???
  priority_queue<pair<double, int>> queue_priority; // distance, proxy

  // initialize proxies
  int proxy_num = Proxies.rows()/2;
  Vector3d* Proxies_center = new Vector3d[proxy_num];
  Vector3d* Proxies_normal = new Vector3d[proxy_num];
  //Vector3d Proxies_center[proxy_num];
  //Vector3d Proxies_normal[proxy_num];
  VectorXi FittedTriangles = find_best_triangles(Regions, Proxies, Vertices, Faces, norme);

  // reset R
  Regions = -MatrixXi::Ones(face_num, 1);
  for (int proxy_index = 0; proxy_index < proxy_num; ++proxy_index) {
    Proxies_center[proxy_index] = Proxies.row(proxy_index);
    Proxies_normal[proxy_index] = Proxies.row(proxy_index + proxy_num);
    Regions(FittedTriangles(proxy_index)) = proxy_index;
    // error += distance(FittedTriangles(proxy_index), Proxies_center[proxy_index], Proxies_normal[proxy_index], Vertices, norme); ???
    for (int ad = 0; ad < 3; ++ad) {
        int ad_face_index = FaceAdjacencies(FittedTriangles(proxy_index),ad);
        double error_metric = distance(ad_face_index, Proxies_center[proxy_index], Proxies_normal[proxy_index], Vertices, norme);
        queue_priority.push(make_pair(-error_metric, ad_face_index + face_num * proxy_index));
    }
  }

  pair<double, int> item;
  int face_index;
  int proxy_index;

  // number of cycles = (face number - proxy number) * 3
  while (queue_priority.size()!=0) {
    item = queue_priority.top();
    queue_priority.pop();
    proxy_index = item.second/face_num; // index of seed Triangles
    face_index = item.second%face_num; // adjacency triangle of seed Triangles

    if (Regions(face_index) == -1) {

      Regions(face_index) = proxy_index;
      for (int ad = 0; ad < 3; ++ad) {
        int ad_face_index = FaceAdjacencies(face_index,ad);
        // if(Regions(ad_face_index) == -1) {
          double error_metric = distance(ad_face_index, Proxies_center[proxy_index], Proxies_normal[proxy_index], Vertices, norme);
          queue_priority.push(make_pair(-error_metric, ad_face_index + face_num * proxy_index)); // Number of cycles dependency
        // }
      }
      
    }
  }

}

