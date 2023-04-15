#include "anchors.h"

/*
    return number of key in vertor v.
*/
bool vector_contains(const vector<int>& v, int key) {
    return count(v.begin(), v.end(), key);
}

bool vector_contains(const vector<double>& v, double key) {
    return count(v.begin(), v.end(), key);
}


int find_first_edge(HalfedgeDS& he, int anchor_vertex_index, int around_proxy_index, const MatrixXi& Regions) {
    int edge = he.getOpposite(he.getEdge(anchor_vertex_index));
    while (Regions(he.getFace(edge), 0) != around_proxy_index || Regions(he.getFace(he.getOpposite(edge)), 0) == around_proxy_index) {
        edge = he.getNext(he.getOpposite(edge));
    }
    return edge;
}
int find_next_edge(HalfedgeDS& he, int edge, int around_proxy_index, const MatrixXi& R) {
    edge = he.getNext(edge);
    while (R(he.getFace(he.getOpposite(edge)), 0) == around_proxy_index) {
        edge = he.getNext(he.getOpposite(edge));
    }
    return edge;
}

list <int> explore_boundary(HalfedgeDS& he, int v, int r, const MatrixXi& R) {
    list <int> bound;
    int edge = find_first_edge(he,v,r,R);
    int new_v = he.getTarget(edge);
    bound.push_back(new_v);
    while (new_v != v) {
        edge = find_next_edge(he,edge,r,R);
        new_v = he.getTarget(edge);
        bound.push_back(new_v);
    }
    return bound;
}


int add_anchors_on_edge(HalfedgeDS& he, int anchor_vertex_index, int around_proxy_index, const MatrixXi& Regions, const MatrixXd& Vertices, VectorXi& Anchors, const MatrixXd& Proxies, double treshold) {
    // modifies anchors to add anchors on edge
    int k = 0; //number of new anchors

    // find all anchors of the region R
    list<int>  anchors_list;
    int fedge = find_first_edge(he, anchor_vertex_index, around_proxy_index, Regions);
    int edge = fedge;
    int new_anchor_index = he.getTarget(edge);
    while (new_anchor_index != anchor_vertex_index) {
        if (Anchors(new_anchor_index) == 1) anchors_list.push_back(new_anchor_index);
        edge = find_next_edge(he, edge, around_proxy_index, Regions);
        new_anchor_index = he.getTarget(edge);
    }
    anchors_list.push_back(anchor_vertex_index);

    // go through each pair of anchor
    edge = fedge;
    int r2;
    while (!anchors_list.empty()) {
        new_anchor_index = anchors_list.front();
        anchors_list.pop_front();
        int new_v = he.getTarget(edge);
        double max_d = 0;
        int max_v = -1;
        while (new_v!= new_anchor_index) {
            edge = find_next_edge(he, edge, around_proxy_index, Regions);
            r2 = Regions(he.getFace(he.getOpposite(edge)), 0);
            new_v = he.getTarget(edge);
            double dist = distance_projection(Vertices, Proxies, anchor_vertex_index, new_anchor_index, new_v, around_proxy_index, r2);
            if (dist > max_d) {
                max_d = dist;
                max_v = new_v;
            }
        }
        if (max_d > treshold) { // treshold
            Anchors(max_v) = 1;
            k++;
        }
        anchor_vertex_index = new_anchor_index;
    }
    return k;
}

vector<int> find_anchors_on_region(HalfedgeDS& he, const int anchor, const int r, const MatrixXi& R, const VectorXi& anchors) {
    // modifies anchors to add anchors on edge
    int k = 0; //number of new anchors

    // find all anchors of the region R
    vector<int>  anchors_vector;
    anchors_vector.push_back(anchor);
    int fedge = find_first_edge(he,anchor,r,R);
    int edge = fedge;
    int new_anchor = he.getTarget(edge);
    while (new_anchor != anchor) {
        if (anchors(new_anchor)==1) {
            anchors_vector.push_back(new_anchor);
        }
        edge = find_next_edge(he,edge,r,R);
        new_anchor = he.getTarget(edge);
    }
    return anchors_vector;
}

// find proxies around a vertex by its index
vector<int> find_vertex_proxies(HalfedgeDS& he, int vertex_index, const MatrixXi& Regions) {
    vector<int> vertex_proxies;
    int proxy_index;
    int edge = he.getEdge(vertex_index);
    proxy_index = Regions(he.getFace(edge), 0);
    vertex_proxies.push_back(proxy_index);
    int pedge = he.getOpposite(he.getNext(edge));
    while (pedge != edge) {
        proxy_index = Regions(he.getFace(pedge), 0);
        if (!vector_contains(vertex_proxies, proxy_index)) {
            vertex_proxies.push_back(proxy_index);
        }
        pedge = he.getOpposite(he.getNext(pedge));
    }
    return vertex_proxies;
}


vector<vector<int>> anchor_points(HalfedgeDS& he, const MatrixXi& Regions, const MatrixXd& Vertices, const MatrixXd& Proxies, double treshold) { 
    int vertex_num = Vertices.rows();
    int proxy_num = Proxies.rows()/2;
    vector<vector<int>> vertex_proxies(vertex_num); //list of proxies
    VectorXi Anchors;
    VectorXi Seen; 
    Seen.setZero(proxy_num);
    Anchors.setZero(vertex_num);
    int anchor_num = 0; // number of Anchors

    // find for each vertex its list of proxies
    for (int vertex_index = 0; vertex_index < vertex_num; ++vertex_index) {
        vertex_proxies[vertex_index] = find_vertex_proxies(he, vertex_index, Regions);
        // add anchor if vertex has 3+ proxies
        if (vertex_proxies[vertex_index].size()>2) {
            Anchors(vertex_index) = 1;
            ++anchor_num; 
        }
    }

    // add Anchors between existing ones
    int proxy_index, kv;
    for (int vertex_index = 0; vertex_index < vertex_num; ++vertex_index) {
        if (vertex_proxies[vertex_index].size() > 2) {
            // to show the max deviation
            for(size_t index = 0; index < vertex_proxies[vertex_index].size(); index++) {
                proxy_index = vertex_proxies[vertex_index][index];
                if (Seen(proxy_index) == 0) {
                    Seen(proxy_index) = 1;
                    kv = add_anchors_on_edge(he, vertex_index, proxy_index, Regions, Vertices, Anchors, Proxies, treshold);
                    while (kv > 0) {
                        anchor_num += kv;
                        kv = add_anchors_on_edge(he, vertex_index, proxy_index, Regions, Vertices, Anchors, Proxies, treshold);
                    }
                }
            }
        }
    }
     Seen.setZero(proxy_num);


     // return list of polygons
     vector<vector<int>> proxy_polygons;
     vector<int> polygon_anchor_vertices;
     for (int index = 0; index < proxy_num; ++index) proxy_polygons.push_back(polygon_anchor_vertices);
     for (int vertex_index = 0; vertex_index < vertex_num; ++vertex_index) {
         if (vertex_proxies[vertex_index].size() > 2) {
             for(size_t index = 0; index < vertex_proxies[vertex_index].size(); ++index) {
                 proxy_index = vertex_proxies[vertex_index][index];
                 // if (Seen(proxy_index)==0) {
                 //     Seen(proxy_index)=1;
                 //     proxy_polygons[proxy_index] = find_anchors_on_region(he, vertex_index, proxy_index, Regions, Anchors);
                 // }
                 if (!vector_contains(proxy_polygons[proxy_index], vertex_index)) {
                     vector<int> new_anchors = find_anchors_on_region(he, vertex_index, proxy_index, Regions, Anchors);
                     proxy_polygons[proxy_index].insert(proxy_polygons[proxy_index].end(), new_anchors.begin(), new_anchors.end() );
                 }
             }
         }
     }
     return proxy_polygons;

    //return vertex_proxies;
}

