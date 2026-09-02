#pragma once

#include "implicit_point.h"
#include <cstring>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>


// Uncommenting the following macro definition makes the code use modified parts of hxt_SeqDel (Copyright (C) 2018 Célestin Marot).
// hxt_SeqDel is a sequential Delaunay triangulator hosted at https://git.immc.ucl.ac.be/hextreme/hxt_seqdel as of 2020.
// hxt_SeqDel is GPL licensed, meaning that if you uncomment the following line you accept the terms of the GPL license for
// the whole code which uses this library.
// If you need to use this code under the less restrictive LGPL license, please comment the following line.
// This will make the code slightly slower.
//#define USE_MAROTS_METHOD

#pragma intrinsic(fabs)

#define INFINITE_VERTEX UINT32_MAX // Vertex at infinity
#define DT_UNKNOWN  0              // Delaunay tetrahedron (DT) is marked as unknown
#define DT_OUT  1                  // Delaunay tetrahedron (DT) is marked as internal
#define DT_IN  2                   // Delaunay tetrahedron (DT) is marked as external

// Macros to set/unset/check a generic bit in the mark mask of tetrahedra
#define MARKBIT(m, twoPowBit) m |= ((uint32_t)(twoPowBit))
#define UNMARKBIT(m, twoPowBit) m &= (~((uint32_t)(twoPowBit)))
#define ISMARKEDBIT(m, twoPowBit) m & ((uint32_t)(twoPowBit)) 

// Prototype conversion for basic geometric predicates
inline int orient3d(const double* a, const double* b, const double* c, const double* d) {
    return orient3d(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], d[0], d[1], d[2]);
}

inline int insphere(const double* a, const double* b, const double* c, const double* d, const double* e) {
    return inSphere(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], d[0], d[1], d[2], e[0], e[1], e[2]);
}

// The most basic 3D point type that can be used to instantiate TetMesh_t
// without using implicit points.
class basicVec3d {
public:
    double coord[3]; // The three Cartesian coordinates

    basicVec3d() {}
    basicVec3d(double x, double y, double z) : coord{ x, y, z } {}

    double X() const { return coord[0]; }
    double Y() const { return coord[1]; }
    double Z() const { return coord[2]; }

    // For compatibility with implicit points
    bool getExactXYZCoordinates(bigrational& _x, bigrational& _y, bigrational& _z) const { 
        _x = bigfloat(X()); _y = bigfloat(Y()); _z = bigfloat(Z()); return true; 
    }

    bool getApproxXYZCoordinates(double& x, double& y, double& z, bool dummy) const {
        x = coord[0];
        y = coord[1];
        z = coord[2];
        return true;
    }

    //
    // EXACT predicates -----------------------------------------------------------------------
    //

    inline bool operator<(const basicVec3d& v) const {
        return (coord[0] < v.coord[0] || 
            (coord[0] == v.coord[0] && coord[1] < v.coord[1]) || 
            (coord[0] == v.coord[0] && coord[1] == v.coord[1] && coord[2] < v.coord[2]));
    }

    // Prototype conversion for basic geometric predicates
    static int orient3D(const basicVec3d& ap, const basicVec3d& bp, const basicVec3d& cp, const basicVec3d& dp) {
        return orient3d(ap.coord, bp.coord, cp.coord, dp.coord);
    }

    static int inSphere(const basicVec3d& ap, const basicVec3d& bp, const basicVec3d& cp, const basicVec3d& dp, const basicVec3d& ep) {
        return insphere(ap.coord, bp.coord, cp.coord, dp.coord, ep.coord);
    }

    // Derived predicates. Implemented based on model functions from Indirect Predicates library.
    static int orient2Dxy(const basicVec3d& a, const basicVec3d& b, const basicVec3d& c) { return orient2d(a.X(), a.Y(), b.X(), b.Y(), c.X(), c.Y()); }
    static int orient2Dyz(const basicVec3d& a, const basicVec3d& b, const basicVec3d& c) { return orient2d(a.Y(), a.Z(), b.Y(), b.Z(), c.Y(), c.Z()); }
    static int orient2Dzx(const basicVec3d& a, const basicVec3d& b, const basicVec3d& c) { return orient2d(a.Z(), a.X(), b.Z(), b.X(), c.Z(), c.X()); }
    static int inGabrielSphere(const basicVec3d& q, const basicVec3d& a, const basicVec3d& b, const basicVec3d& c);
    static int incirclexy(const basicVec3d& a, const basicVec3d& b, const basicVec3d& c, const basicVec3d& d);
    static int dotProductSign3D(const basicVec3d& a, const basicVec3d& b, const basicVec3d& c);
    static int lessThanOnX(const basicVec3d& a, const basicVec3d& b) { return (a.X() > b.X()) - (a.X() < b.X()); }
    static int lessThanOnY(const basicVec3d& a, const basicVec3d& b) { return (a.Y() > b.Y()) - (a.Y() < b.Y()); }
    static int lessThanOnZ(const basicVec3d& a, const basicVec3d& b) { return (a.Z() > b.Z()) - (a.Z() < b.Z()); }
    static int lessThan(const basicVec3d& a, const basicVec3d& b);
    static bool coincident(const basicVec3d& a, const basicVec3d& b) { return a.X() == b.X() && a.Y() == b.Y() && a.Z() == b.Z(); }
    static bool misaligned(const basicVec3d& A, const basicVec3d& B, const basicVec3d& C);
    static bool pointInInnerSegment(const basicVec3d& p, const basicVec3d& v1, const basicVec3d& v2);
    static bool pointInSegment(const basicVec3d& p, const basicVec3d& v1, const basicVec3d& v2);
    static bool pointInInnerTriangle(const basicVec3d& P, const basicVec3d& A, const basicVec3d& B, const basicVec3d& C);
    static bool pointInTriangle(const basicVec3d& P, const basicVec3d& A, const basicVec3d& B, const basicVec3d& C);
    static bool pointInTriangle(const basicVec3d& P, const basicVec3d& A, const basicVec3d& B, const basicVec3d& C, int& oAB, int& oBC, int& oCA);
    static bool innerSegmentsCross(const basicVec3d& A, const basicVec3d& B, const basicVec3d& P, const basicVec3d& Q);
    static bool segmentsCross(const basicVec3d& A, const basicVec3d& B, const basicVec3d& P, const basicVec3d& Q);
    static bool innerSegmentCrossesInnerTriangle(const basicVec3d& s1, const basicVec3d& s2, const basicVec3d& v1, const basicVec3d& v2, const basicVec3d& v3);
    static bool lineCrossesInnerTriangle(const basicVec3d& s1, const basicVec3d& s2, const basicVec3d& v1, const basicVec3d& v2, const basicVec3d& v3);
    static bool lineCrossesTriangle(const basicVec3d& s1, const basicVec3d& s2, const basicVec3d& v1, const basicVec3d& v2, const basicVec3d& v3);
    static bool innerSegmentCrossesTriangle(const basicVec3d& s1, const basicVec3d& s2, const basicVec3d& v1, const basicVec3d& v2, const basicVec3d& v3);
    static int orient2D(const basicVec3d& a, const basicVec3d& b, const basicVec3d& c, int n_max);
    static bool misaligned(const basicVec3d& A, const basicVec3d& B, const basicVec3d& C, int n_max);
    static bool pointInInnerSegment(const basicVec3d& p, const basicVec3d& v1, const basicVec3d& v2, int n_max);
    static bool pointInSegment(const basicVec3d& p, const basicVec3d& v1, const basicVec3d& v2, int n_max);
    static bool pointInInnerTriangle(const basicVec3d& P, const basicVec3d& A, const basicVec3d& B, const basicVec3d& C, int n_max);
    static bool pointInTriangle(const basicVec3d& P, const basicVec3d& A, const basicVec3d& B, const basicVec3d& C, int n_max);
    static bool innerSegmentsCross(const basicVec3d& A, const basicVec3d& B, const basicVec3d& P, const basicVec3d& Q, int n_max);
    static bool segmentsCross(const basicVec3d& A, const basicVec3d& B, const basicVec3d& P, const basicVec3d& Q, int n_max);

    // INEXACT operations -------------------------------------------------------------

    inline basicVec3d operator+(const basicVec3d& v) const { return basicVec3d(coord[0] + v.coord[0], coord[1] + v.coord[1], coord[2] + v.coord[2]); }
    inline basicVec3d operator-(const basicVec3d& v) const { return basicVec3d(coord[0] - v.coord[0], coord[1] - v.coord[1], coord[2] - v.coord[2]); }
    inline basicVec3d operator*(const double d) const { return basicVec3d(coord[0] * d, coord[1] * d, coord[2] * d); }
    inline void operator+=(const basicVec3d& v) { coord[0] += v.coord[0]; coord[1] += v.coord[1]; coord[2] += v.coord[2]; }
    inline void operator*=(const double d) { coord[0] *= d; coord[1] *= d; coord[2] *= d; }

    inline double dot(const basicVec3d& p) const { return (coord[0] * p.coord[0] + coord[1] * p.coord[1] + coord[2] * p.coord[2]); }
    inline basicVec3d cross(const basicVec3d& p) const { return basicVec3d(coord[1] * p.coord[2] - coord[2] * p.coord[1], coord[2] * p.coord[0] - coord[0] * p.coord[2], coord[0] * p.coord[1] - coord[1] * p.coord[0]); }
    inline double tripleProd(const basicVec3d& v2, const basicVec3d& v3) const {
        return ((v2.coord[0] * v3.coord[1] * coord[2]) - (v3.coord[0] * v2.coord[1] * coord[2])) +
            ((v3.coord[0] * coord[1] * v2.coord[2]) - (coord[0] * v3.coord[1] * v2.coord[2])) +
            ((coord[0] * v2.coord[1] * v3.coord[2]) - (v2.coord[0] * coord[1] * v3.coord[2]));
    }

    inline double operator*(const basicVec3d& d) const { return dot(d); }
    inline basicVec3d operator&(const basicVec3d& d) const { return cross(d); }

    // Squared length
    inline double sq_length() const { return dot(*this); }

    // Squared distance
    inline double dist_sq(const basicVec3d& v) const { return ((*this) - v).sq_length(); }

    static inline bool inSmallestSphere(const basicVec3d& pv, const basicVec3d& qv, const basicVec3d& rv) {
        return ((rv - pv).sq_length() + (rv - qv).sq_length()) <= (pv - qv).sq_length();
    }

    static inline bool hasLargerSphere(const basicVec3d& pv, const basicVec3d& qv, const basicVec3d& rv, const basicVec3d& sv) {
        const basicVec3d pms = pv - sv, qms = qv - sv, pmr = pv - rv, qmr = qv - rv;
        const double lens = pms.sq_length() * qms.sq_length();
        if (lens == 0) return true;
        const double lenr = pmr.sq_length() * qmr.sq_length();
        if (lenr == 0) return false;
        const double dots = pms.dot(qms);
        const double dotr = pmr.dot(qmr);

        return (dots * dots) * lenr < (dotr * dotr) * lens;
    }

    // TRUE if p is closer to q than to r
    static bool isCloserThan(const basicVec3d& pv, const basicVec3d& qv, const basicVec3d& rv) {
        return pv.dist_sq(qv) < pv.dist_sq(rv);
    }

    // TRUE if distance p-q is at most half the distance p-r
    static bool isAtMostTwiceDistanceThan(const basicVec3d& pv, const basicVec3d& qv, const basicVec3d& rv) {
        return pv.dist_sq(qv) * 4 < pv.dist_sq(rv);
    }
};

inline std::ostream& operator<<(std::ostream& os, const basicVec3d& p)
{
    return os << p.coord[0] << " " << p.coord[1] << " " << p.coord[2];
}


// Tetrahedron facet - used to reconstruct connectivity from indexed arrays

class tetFacet
{
public:
    uint32_t v[3]; // The three facet vertices
    uint64_t t; // Its incident tet's opposite corner
    tetFacet(uint32_t a, uint32_t b, uint32_t c, uint64_t d)
        : t(d)
    {
        v[0] = std::min(a, std::min(b, c));
        v[2] = std::max(a, std::max(b, c));
        v[1] = (a != v[0] && a != v[2]) ? (a) : ((b != v[0] && b != v[2]) ? (b) : (c));
    }

    bool operator==(const tetFacet& f) const
    {
        return v[0] == f.v[0] && v[1] == f.v[1] && v[2] == f.v[2];
    }

    bool operator<(const tetFacet& f) const
    {
        if (v[0] < f.v[0])
            return true;
        else if (v[0] > f.v[0])
            return false;
        else if (v[1] < f.v[1])
            return true;
        else if (v[1] > f.v[1])
            return false;
        else
            return (v[2] < f.v[2]);
    }
};

// Tetrahedral mesh data structure

template<class pointType>
class TetMesh_t {
public:
  // General purpose fields
  std::vector<pointType> vertices; // Vertices
  std::vector<uint64_t> inc_tet; // One tet incident upon each vertex
  std::vector<uint32_t> tet_node; // Tet corners
  std::vector<uint64_t> tet_neigh; // Tet opposites
  mutable std::vector<uint32_t> mark_tetrahedra; // Marks on tets
  mutable std::vector<uint8_t> marked_vertex; // Marks on vertices

  std::vector<uint64_t> Del_deleted; // Unlinked tets to be removed from the structure

  const bool has_outer_vertices; // This is TRUE if mesh vertices must survive after destruction

  double max_coord_x, max_coord_y, max_coord_z , max_coord;
  double o3d_static_filter; // Static filter for orient3d
  double isp_static_filter; // Static filter for insphere
  std::vector<double> tet_subdet; // Tet sub-determinants to speedup insphere tests

  // Constructor and destructor
  TetMesh_t() : has_outer_vertices(false), o3d_static_filter(DBL_MAX), isp_static_filter(DBL_MAX) {};
  TetMesh_t(bool h) : has_outer_vertices(h), o3d_static_filter(DBL_MAX), isp_static_filter(DBL_MAX) {};
  ~TetMesh_t() { if (!has_outer_vertices) flushVertices(); };


  /////// Global functions ///////

  // Number of vertices (infinite vertex is not counted)
  uint32_t numVertices() const { return (uint32_t)vertices.size(); }

  // Number of tetrahedra including ghosts
  uint32_t numTets() const { return (uint32_t)(tet_node.size() >> 2); }

  // Number of non-ghost tetrahedra
  uint32_t countNonGhostTets() const {
      return numTets() - (uint32_t)std::count(tet_node.begin(), tet_node.end(), INFINITE_VERTEX);
  }

  // Fill the vertex vector with newly-created objects (specific behaviour depends on pointType)
  void init_vertices(const double* coords, uint32_t num_v);

  // Init the connectivity structure starting from a serialized tet array
  void init_tets(const uint32_t* tet_idx, size_t num_t);

  // Init the static filters for geometric predicates
  void init_static_filters();
  void update_static_filters();

  // Destroy vertices
  void flushVertices() { if constexpr (std::is_pointer_v<pointType>) for (const pointType& p : vertices) delete p; }

  // Init the mesh with a tet connecting four non coplanar points in vertices
  void init(uint32_t& unswap_k, uint32_t& unswap_l);

  // Create a Delaunay tetrahedrization by incremental insertion
  void tetrahedrize();

  // Save the mesh to a .tet file
  // If inner_only is set, only tets tagged as DT_IN are saved
  bool saveTET(const char* filename, bool inner_only = false) const;

  // Save the mesh to a .vtu file (ParaView compatible - XML)
  bool saveVTU(const char* filename, bool inner_only = false) const;

  // Save the mesh to a .mesh file (MEDIT format)
  // If inner_only is set, only tets tagged as DT_IN are saved
  bool saveMEDIT(const char* filename, bool inner_only = false) const;

  // As above, but uses a binary format to avoid rounding
  bool saveBinaryTET(const char* filename, bool inner_only = false) const;

  // Save the interface between DT_IN and DT_OUT as an OFF file
  bool saveBoundaryToOFF(const char* filename) const;

  // As above, but saves rational coordinates and distinguishes between inner and outer tets
  bool saveRationalTET(const char* filename, bool inner_only = false);

  // Marks internal tets as DT_IN and ghost tets as DT_OUT
  void markInnerTets(uint64_t single_start = UINT64_MAX) {
      for (size_t i = 0; i < numTets(); i++)
          mark_tetrahedra[i] = (isGhost(i)) ? DT_OUT : DT_IN;
  }

  // Clear deleted tets after insertions
  void removeDelTets();

  // Resize the whole structure to contain 'new_size' tets
  void resizeTets(uint64_t new_size);
  void reserveTets(uint64_t new_capacity);

  // Return TRUE if at least one tet becomes flat or inverted after having
  // snapped its vertices to their closest floating-point representable positions.
  // Init num_flipped and num_flattened with the overall number of flips or flattings.
  bool hasBadSnappedOrientations(size_t& num_flipped, size_t& num_flattened) const;

  // Check whether the structure is coherent (use for debugging purposes)
  void checkMesh(bool checkDelaunay =true) const;

  /////// Local (element-based) functions ///////

  // TRUE if tet is ghost
  bool isGhost(uint64_t t) const { return tet_node[(t << 2) + 3] == INFINITE_VERTEX; }

  // TRUE if t has vertex v
  bool tetHasVertex(uint64_t t, uint32_t v) const;

  // Init 'ov' with the two vertices of tet which are not in 'v'
  void oppositeTetEdge(const uint64_t tet, const uint32_t v[2], uint32_t ov[2]) const;

  // Let t and n be face-adjacent tets.
  // This function returns the corner in t which is opposite to n
  uint64_t getCornerFromOppositeTet(uint64_t t, uint64_t n) const;

  // Return the i'th tet in neighbors 'n'
  inline uint64_t getIthNeighbor(const uint64_t* n, const uint64_t i) const { return n[i] & (~3); }

  // Fill v with the three vertices of t	
  void getFaceVertices(uint64_t t, uint32_t v[3]) const;

  // Same as above, but vertices are sorted based on their indexes
  void getFaceSortedVertices(uint64_t t, uint32_t v[3]) const;

  // Fill 'nt' with the two tets that share the vertices v1,v2,v3
  bool getTetsFromFaceVertices(uint32_t v1, uint32_t v2, uint32_t v3, uint64_t* nt) const;

  // Return the corner of t which is opposite to its face with vertices v1,v2,v3
  uint64_t tetOppositeCorner(uint64_t t, uint32_t v1, uint32_t v2, uint32_t v3) const;

  // Return the corner corresponding to vertex 'v' in the tet whose base corner is tb
  uint64_t tetCornerAtVertex(uint64_t tb, uint32_t v) const {
      return ((tet_node[tb] == v) * (tb)) + ((tet_node[tb + 1] == v) * (tb + 1)) + ((tet_node[tb + 2] == v) * (tb + 2)) + ((tet_node[tb + 3] == v) * (tb + 3));

      //while (tet_node[tb] != v) tb++;
      //return tb;
  }

  // Set the adjacency between the two corners c1 and c2
  void setMutualNeighbors(const uint64_t c1, const uint64_t c2) { tet_neigh[c1] = c2; tet_neigh[c2] = c1; }

  // Direct pointer to nodes and neighs
  uint32_t* getTetNodes(uint64_t tet) { return tet_node.data() + tet; }
  uint64_t* getTetNeighs(uint64_t tet) { return tet_neigh.data() + tet; }
  const uint32_t* getTetNodes(uint64_t tet) const { return tet_node.data() + tet; }
  const uint64_t* getTetNeighs(uint64_t tet) const { return tet_neigh.data() + tet; }

  // tetNi is a sum modulo 4 - used to traverse the nodes of a tet
  static size_t tetN1(const size_t i) { return (i + 1) & 3; }
  static size_t tetN2(const size_t i) { return (i + 2) & 3; }
  static size_t tetN3(const size_t i) { return (i + 3) & 3; }

  // tetONi - as above, but results in a coherent orientation
  static size_t tetON1(const size_t i) { return tetN1(i); }
  static size_t tetON2(const size_t i) { return (i & 2) ^ 3; }
  static size_t tetON3(const size_t i) { return (i + 3) & 2; }

  // Push a new isolated vertex in the structure
  void pushVertex(const pointType& p) {
      vertices.push_back(p);
      inc_tet.push_back(UINT64_MAX);
      marked_vertex.push_back(0);

      if constexpr (std::is_same_v<pointType, basicVec3d>) {
          double c;
          c = fabs(p.coord[0]); if (c > max_coord_x) { max_coord_x = c; if (c > max_coord) max_coord = c; }
          c = fabs(p.coord[1]); if (c > max_coord_y) { max_coord_y = c; if (c > max_coord) max_coord = c; }
          c = fabs(p.coord[2]); if (c > max_coord_z) { max_coord_z = c; if (c > max_coord) max_coord = c; }
          update_static_filters();
      }
  }

  // Pop an isolated vertex from the structure
  pointType popVertex() {
      pointType p = vertices.back();
      vertices.pop_back();
      inc_tet.pop_back();
      marked_vertex.pop_back();
      return p;
  }

  // Inserts an isolated vertex which is already in the vertices array.
  // ct is a hint for the algorithm to start searching the tet containing vi
  void insertExistingVertex(const uint32_t vi, uint64_t& ct);

  // Same as above, but only split the tet(s) containing the point
  // without forcing the Delaunay condition. Differently from the
  // above version, in this case 'ct' is assumed to contain the point.
  void insertExistingVertexNonDelaunay(const uint32_t vi, uint64_t ct);

  // fill 'cavityCorners' with tets bounding v_id's cavity
  void getUnconstrainedCavity(const uint32_t v_id, const uint64_t tet, std::vector<uint64_t>& cavityCorners);
  // fill 'adjacencies' with consecutive pairs of edge-adjacent tets in cavityCorners
  void getCavityConnectivity(const std::vector<uint64_t>& cavityCorners, std::vector<uint64_t>& adjacencies) const;
  // replace tetrahedra bounded by cavityCorners with new tets around vi
  void retetrahedrizeCavity(const uint32_t v_id, const std::vector<uint64_t>& cavityCorners, const std::vector<uint64_t>& adjacencies, const uint32_t mark);

  // Starting from 'tet', move by adjacencies until a tet is found that
  // contains vi. Return that tet.
  uint64_t searchTetrahedron(uint64_t tet, const uint32_t v_id);

  // Incident tetrahedra at a vertex
  void VT(uint32_t v, std::vector<uint64_t>& vt) const;

  // Same as VT, but this one includes ghost tets as well
  void VTfull(uint32_t v, std::vector<uint64_t>& vt) const;

  // Adjacent vertices
  void VV(uint32_t v, std::vector<uint32_t>& vv) const;

  // Incident tetrahedra at an edge
  void ET(uint32_t v1, uint32_t v2, std::vector<uint64_t>& et) const; // Non-ghost tets incident at v1-v2
  void ETfull(uint32_t v1, uint32_t v2, std::vector<uint64_t>& et) const; // Same as ET, but includes ghost tets
  void ETfast(uint32_t v1, uint32_t v2, uint64_t st_t, std::vector<uint64_t>& et) const; // Same as ET, but starts from a known incident tet st_t (faster)
  void ETfullfast(uint32_t v1, uint32_t v2, uint64_t st_t, std::vector<uint64_t>& et) const; // Combination of the previous two

  // Incident tetrahedra at an edge represented as ordered sequence of corners
  void ETcorners(uint32_t v1, uint32_t v2, std::vector<uint64_t>& et) const;

  // TRUE if v1 and v2 are connected by an edge
  bool hasEdge(uint32_t v1, uint32_t v2) const;

  // Swap the position of t1 and t2 in the structure and update all relations accordingly
  void swapTets(const uint64_t t1, const uint64_t t2);

  // Mark/unmark/check one single bit in tet mask
  inline void mark_Tet_1(const uint64_t t) const { mark_tetrahedra[t] |= ((uint32_t)2); }
  inline void unmark_Tet_1(const uint64_t t) const { mark_tetrahedra[t] &= (~((uint32_t)2)); }
  inline uint32_t is_marked_Tet_1(const uint64_t t) const { return mark_tetrahedra[t] & ((uint32_t)2); }
  inline void mark_Tet_2(const uint64_t t) const { mark_tetrahedra[t] |= ((uint32_t)4); }
  inline void unmark_Tet_2(const uint64_t t) const { mark_tetrahedra[t] &= (~((uint32_t)4)); }
  inline uint32_t is_marked_Tet_2(const uint64_t t) const { return mark_tetrahedra[t] & ((uint32_t)4); }
  inline void mark_Tet_31(const uint64_t t) const { mark_tetrahedra[t] |= ((uint32_t)2147483648); }
  inline void unmark_Tet_31(const uint64_t t) const { mark_tetrahedra[t] &= (~((uint32_t)2147483648)); }
  inline uint32_t is_marked_Tet_31(const uint64_t t) const { return mark_tetrahedra[t] & ((uint32_t)2147483648); }

  // Thes two functions mark/check one particular bit stating that a tet must be deleted.
  // Differently from above, here a tet is identified by its first corner.
  void markToDelete(uint64_t c) { mark_tetrahedra[c >> 2] |= ((uint32_t)1073741824); }
  void unmarkToDelete(uint64_t c) { mark_tetrahedra[c >> 2] &= (~(uint32_t)1073741824); }
  bool isToDelete(uint64_t c) const { return mark_tetrahedra[c >> 2] & ((uint32_t)1073741824); }
  bool isToDeleteSmall(uint64_t c) const { return mark_tetrahedra[c] & ((uint32_t)1073741824); }
  void moveDeletedToTail(uint64_t t, uint64_t l);
  void undeleteTets(size_t num);

  // Marks a tet (identified by its first corner) as 'removed' and add it to the queue
  // for eventual deletion.
  void pushAndMarkDeletedTets(uint64_t c) { Del_deleted.push_back(c); markToDelete(c); }

  // Predicates operating on vertex indexes
  int vOrient3D(uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4) const {
      if constexpr (std::is_pointer_v<pointType>) 
          return -std::remove_pointer_t<pointType>::orient3D(*vertices[v1], *vertices[v2], *vertices[v3], *vertices[v4]);
      else return -pointType::orient3D(vertices[v1], vertices[v2], vertices[v3], vertices[v4]);
  }

  int vInSphere(uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4, uint32_t v5) const {
      if constexpr (std::is_pointer_v<pointType>)
          return -std::remove_pointer_t<pointType>::inSphere(*vertices[v1], *vertices[v2], *vertices[v3], *vertices[v4], *vertices[v5]);
      else return -pointType::inSphere(vertices[v1], vertices[v2], vertices[v3], vertices[v4], vertices[v5]);
  }

  // Use the order of the five cospherical points in 'indices' to
  // return a nonzero though coherent inSphere result.
  int symbolicPerturbation(uint32_t indices[5]) const;

  // This is as vInSphere(v[0], v[1], v[2], v[3], v_id) but is guaranteed to
  // return a nonzero value by relying on the symbolic perturbation above.
  int vertexInTetSphere(const uint32_t v[4], uint32_t v_id) const;

  // Same as above, but the four vertices are the vertices of 'tet'.
  int vertexInTetSphere(uint64_t tet, uint32_t v_id) const;

  // Compute tet's sub-determinant
  void compute_subDet(const uint64_t tet);

#ifdef USE_MAROTS_METHOD

  void deleteInVertexTets(uint64_t tet, const uint32_t v_id);

  class DelTmp {
  public:
      uint32_t node[4];
      uint64_t bnd;

      DelTmp(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint64_t o) :
          node{ a, b, c, d }, bnd(o) {}
  };

  std::vector<DelTmp> Del_tmp;
  uint64_t numDelTmp() const { return Del_tmp.size(); }
  void flushDelTmp() { Del_tmp.clear(); }
  uint64_t* delTmpVec() const { return (uint64_t*)Del_tmp.data(); }
  void bnd_push(uint32_t v_id, uint32_t node1, uint32_t node2, uint32_t node3, uint64_t bnd) {
      Del_tmp.push_back(DelTmp(v_id, node1, node2, node3, bnd));
  }

  void deleteInSphereTets(uint64_t tet, const uint32_t v_id);
  void tetrahedrizeHole(uint64_t* tet);
#endif
};


/////////////////////////////////////////////////////////////////////////////
//
// F U N C T I O N   I M P L E M E N T A T I O N
//
/////////////////////////////////////////////////////////////////////////////

inline bool basicVec3d::misaligned(const basicVec3d& A, const basicVec3d& B, const basicVec3d& C) {
    return misaligned_t<basicVec3d>(A, B, C);
}

inline bool basicVec3d::innerSegmentsCross(const basicVec3d& A, const basicVec3d& B, const basicVec3d& P, const basicVec3d& Q) {
    return innerSegmentsCross_t<basicVec3d>(A, B, P, Q);
}

inline bool basicVec3d::segmentsCross(const basicVec3d& A, const basicVec3d& B, const basicVec3d& P, const basicVec3d& Q) {
    return segmentsCross_t<basicVec3d>(A, B, P, Q);
}

inline bool basicVec3d::innerSegmentCrossesInnerTriangle(const basicVec3d& s1, const basicVec3d& s2, const basicVec3d& v1, const basicVec3d& v2, const basicVec3d& v3) {
    return innerSegmentCrossesInnerTriangle_t<basicVec3d>(s1, s2, v1, v2, v3);
}

inline bool basicVec3d::pointInInnerSegment(const basicVec3d& p, const basicVec3d& v1, const basicVec3d& v2) {
    return pointInInnerSegment_t<basicVec3d>(p, v1, v2);
}

inline bool basicVec3d::pointInSegment(const basicVec3d& p, const basicVec3d& v1, const basicVec3d& v2) {
    return pointInSegment_t<basicVec3d>(p, v1, v2);
}

inline bool basicVec3d::pointInTriangle(const basicVec3d& p, const basicVec3d& a, const basicVec3d& b, const basicVec3d& c) {
    return pointInTriangle_t<basicVec3d>(p, a, b, c);
}

inline bool basicVec3d::pointInTriangle(const basicVec3d& p, const basicVec3d& a, const basicVec3d& b, const basicVec3d& c, int& o1, int& o2, int& o3) {
    return pointInTriangle_t<basicVec3d>(p, a, b, c, o1, o2, o3);
}

inline bool basicVec3d::pointInInnerTriangle(const basicVec3d& p, const basicVec3d& a, const basicVec3d& b, const basicVec3d& c) {
    return pointInInnerTriangle_t<basicVec3d>(p, a, b, c);
}

inline bool basicVec3d::lineCrossesInnerTriangle(const basicVec3d& s1, const basicVec3d& s2, const basicVec3d& v1, const basicVec3d& v2, const basicVec3d& v3) {
    return lineCrossesInnerTriangle_t<basicVec3d>(s1, s2, v1, v2, v3);
}

inline bool basicVec3d::lineCrossesTriangle(const basicVec3d& s1, const basicVec3d& s2, const basicVec3d& v1, const basicVec3d& v2, const basicVec3d& v3) {
    return lineCrossesTriangle_t<basicVec3d>(s1, s2, v1, v2, v3);
}

inline bool basicVec3d::innerSegmentCrossesTriangle(const basicVec3d& s1, const basicVec3d& s2, const basicVec3d& v1, const basicVec3d& v2, const basicVec3d& v3) {
    return innerSegmentCrossesTriangle_t<basicVec3d>(s1, s2, v1, v2, v3);
}

inline int basicVec3d::orient2D(const basicVec3d& a, const basicVec3d& b, const basicVec3d& c, int n_max) {
    return orient2D_t(a, b, c, n_max);
}

inline bool basicVec3d::misaligned(const basicVec3d& A, const basicVec3d& B, const basicVec3d& C, int n_max) {
    return misaligned_t(A, B, C, n_max);
}

inline bool basicVec3d::pointInInnerSegment(const basicVec3d& p, const basicVec3d& v1, const basicVec3d& v2, int xyz) {
    return pointInInnerSegment_t<basicVec3d>(p, v1, v2, xyz);
}

inline bool basicVec3d::pointInSegment(const basicVec3d& p, const basicVec3d& v1, const basicVec3d& v2, int xyz) {
    return pointInSegment_t<basicVec3d>(p, v1, v2, xyz);
}

inline bool basicVec3d::pointInInnerTriangle(const basicVec3d& p, const basicVec3d& a, const basicVec3d& b, const basicVec3d& c, int xyz) {
    return pointInInnerTriangle_t<basicVec3d>(p, a, b, c, xyz);
}

inline bool basicVec3d::pointInTriangle(const basicVec3d& p, const basicVec3d& a, const basicVec3d& b, const basicVec3d& c, int xyz) {
    return pointInTriangle_t<basicVec3d>(p, a, b, c, xyz);
}

inline bool basicVec3d::innerSegmentsCross(const basicVec3d& A, const basicVec3d& B, const basicVec3d& P, const basicVec3d& Q, int xyz) {
    return innerSegmentsCross_t<basicVec3d>(A, B, P, Q, xyz);
}

inline bool basicVec3d::segmentsCross(const basicVec3d& A, const basicVec3d& B, const basicVec3d& P, const basicVec3d& Q, int xyz) {
    return segmentsCross_t<basicVec3d>(A, B, P, Q, xyz);
}

inline int basicVec3d::inGabrielSphere(const basicVec3d& q, const basicVec3d& a, const basicVec3d& b, const basicVec3d& c) {
    return inGabrielSphere_ns(q.X(), q.Y(), q.Z(), a.X(), a.Y(), a.Z(), b.X(), b.Y(), b.Z(), c.X(), c.Y(), c.Z());
}

inline int basicVec3d::incirclexy(const basicVec3d& a, const basicVec3d& b, const basicVec3d& c, const basicVec3d& d) {
    return incircle_ns(a.X(), a.Y(), b.X(), b.Y(), c.X(), c.Y(), d.X(), d.Y());
}

inline int basicVec3d::dotProductSign3D(const basicVec3d& a, const basicVec3d& b, const basicVec3d& c) {
    return dotProductSign3D_ns(a.X(), a.Y(), a.Z(), b.X(), b.Y(), b.Z(), c.X(), c.Y(), c.Z());
}

inline int basicVec3d::lessThan(const basicVec3d& a, const basicVec3d& b) {
    int ret;
    if ((ret = ((a.X() > b.X()) - (a.X() < b.X())))) return ret;
    if ((ret = ((a.Y() > b.Y()) - (a.Y() < b.Y())))) return ret;
    return ((a.Z() > b.Z()) - (a.Z() < b.Z()));
}



template<class pointType> void TetMesh_t<pointType>::init_vertices(const double* coords, uint32_t num_v) {
    vertices.reserve(num_v);

    if constexpr (std::is_same_v<pointType, genericPoint *> || std::is_same_v<pointType, explicitPoint3D *>) {
        for (uint32_t i = 0; i < num_v; i++)
            vertices.push_back(new explicitPoint3D(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]));
    }
    else if constexpr (std::is_same_v<pointType, explicitPoint3D>) {
        for (uint32_t i = 0; i < num_v; i++)
            vertices.push_back(explicitPoint3D(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]));
    }
    else if constexpr (std::is_same_v<pointType, basicVec3d>) {
        for (uint32_t i = 0; i < num_v; i++)
            vertices.push_back(basicVec3d(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]));

        // Calculate static filters
        init_static_filters();
    }

    static_assert(std::is_same_v<pointType, genericPoint*> 
        || std::is_same_v<pointType, explicitPoint3D*>
        || std::is_same_v<pointType, explicitPoint3D>
        || std::is_same_v<pointType, basicVec3d>
        , "TetMesh_t can be instantiated with genericPoint, genericPoint *, explicitPoint3D, explicitPoint3D * or basicVec3d. Other types are not supported.");

    inc_tet.resize(num_v, UINT64_MAX);
    marked_vertex.resize(num_v, 0);
}

template<class pointType> void TetMesh_t<pointType>::init_tets(const uint32_t* tet_indices, size_t _num_tets) {
    resizeTets(_num_tets);
    uint32_t* tet_nodes = tet_node.data();

    for (size_t i = 0; i < _num_tets; i++) {
        uint32_t* nodes = tet_nodes + i * 4;
        nodes[0] = tet_indices[i * 4];
        nodes[1] = tet_indices[i * 4 + 1];
        nodes[2] = tet_indices[i * 4 + 2];
        nodes[3] = tet_indices[i * 4 + 3];
        inc_tet[nodes[0]] = inc_tet[nodes[1]] = inc_tet[nodes[2]] = inc_tet[nodes[3]] = i;
        if (vOrient3D(nodes[0], nodes[1], nodes[2], nodes[3]) < 0) std::swap(nodes[0], nodes[1]);
    }

    // Create tet half-facets and sort them lexycographically
    std::vector<tetFacet> facets;
    facets.reserve(_num_tets * 4);
    for (size_t i = 0; i < _num_tets * 4; i += 4) {
        uint32_t* nodes = tet_nodes + i;
        facets.push_back(tetFacet(nodes[0], nodes[1], nodes[2], i + 3));
        facets.push_back(tetFacet(nodes[0], nodes[1], nodes[3], i + 2));
        facets.push_back(tetFacet(nodes[0], nodes[2], nodes[3], i + 1));
        facets.push_back(tetFacet(nodes[3], nodes[1], nodes[2], i));
    }
    std::sort(facets.begin(), facets.end());

    // Make adjacencies
    std::fill(tet_neigh.begin(), tet_neigh.end(), UINT64_MAX);

    for (size_t i = 1; i < _num_tets * 4; i++) {
        const tetFacet& f1 = facets[i - 1];
        const tetFacet& f2 = facets[i];
        if (f1 == f2) {
            tet_neigh[f1.t] = f2.t;
            tet_neigh[f2.t] = f1.t;
        }
    }

    // Build a structure for ghost tets
    
    // First, count the ghosts (one for each facet having only one incident tet)
    size_t num_ghosts = 0;
    for (size_t i = 0; i < _num_tets * 4; i++)
        if (tet_neigh[facets[i].t] == UINT64_MAX) num_ghosts++;

    // Now, make room for ghosts 
    resizeTets(_num_tets + num_ghosts);
    facets.reserve((_num_tets + num_ghosts) * 4);
    uint64_t last = _num_tets * 4;

    // Build the ghost tets
    for (size_t i = 0; i < _num_tets * 4; i++) {
        const tetFacet& f = facets[i];
        if (tet_neigh[f.t] == UINT64_MAX) {
            facets.push_back(tetFacet(f.v[0], f.v[1], f.v[2], last + 3));
            facets.push_back(tetFacet(f.v[0], f.v[1], INFINITE_VERTEX, last + 2));
            facets.push_back(tetFacet(f.v[0], f.v[2], INFINITE_VERTEX, last + 1));
            facets.push_back(tetFacet(INFINITE_VERTEX, f.v[1], f.v[2], last));
            tet_node[last++] = f.v[0];
            tet_node[last++] = f.v[1];
            tet_node[last++] = f.v[2];
            // The following is ugly!
            if (vOrient3D(f.v[0], f.v[1], f.v[2], tet_node[f.t]) > 0) {
                std::swap(tet_node[last - 1], tet_node[last - 2]);
                std::swap(facets[facets.size() - 2].t, facets[facets.size() - 3].t);
            }
            tet_node[last++] = INFINITE_VERTEX;
        }
    }

            // Qui sotto si pu� fare meglio:
            // Le adiacenze dei tet di bordo con i loro ghost si fanno subito qui sopra
            // Quelle ghost-ghost si fanno riordinando e scorrendo solo la coda del vettore di facce

    // Re-sort facets
    std::sort(facets.begin(), facets.end());

    // Re-build adjacencies to include ghosts
    for (uint64_t i = 1; i < last; i++) {
        const tetFacet& f1 = facets[i - 1];
        const tetFacet& f2 = facets[i];
        if (f1 == f2) {
            tet_neigh[f1.t] = f2.t;
            tet_neigh[f2.t] = f1.t;
        }
    }

    if constexpr (std::is_same_v<pointType, basicVec3d>)
        for (size_t i = 0; i < _num_tets; i++) compute_subDet(i);
}

template<class pointType> void TetMesh_t<pointType>::init_static_filters() {
    if constexpr (std::is_same_v<pointType, basicVec3d>) {
        const uint32_t num_v = (uint32_t)vertices.size();

        // Calculate static filters
        double c;
        double minx = DBL_MAX, miny = DBL_MAX, minz = DBL_MAX;
        double maxx = -DBL_MAX, maxy = -DBL_MAX, maxz = -DBL_MAX;

        for (uint32_t i = 0; i < num_v; i++) {
            const double* cp = vertices[i].coord;
            c = *(cp++);
            if (c < minx) minx = c;
            if (c > maxx) maxx = c;
            c = *(cp++);
            if (c < miny) miny = c;
            if (c > maxy) maxy = c;
            c = *cp;
            if (c < minz) minz = c;
            if (c > maxz) maxz = c;
        }
        maxx = fabs(maxx);
        maxy = fabs(maxy);
        maxz = fabs(maxz);
        minx = fabs(minx);
        miny = fabs(miny);
        minz = fabs(minz);

        max_coord_x = std::max(maxx, minx);
        max_coord_y = std::max(maxy, miny);
        max_coord_z = std::max(maxz, minz);
        max_coord = std::max(max_coord_x, std::max(max_coord_y, max_coord_z));
        update_static_filters();
    }
}

template<class pointType> void TetMesh_t<pointType>::update_static_filters() {
    if constexpr (std::is_same_v<pointType, basicVec3d>) {
        o3d_static_filter = 5.1107127829973299e-15 * max_coord_x * max_coord_y * max_coord_z;
        isp_static_filter = 1.2466136531027298e-13 * max_coord_x * max_coord_y * max_coord_z * max_coord * max_coord;
    }
}

template<class pointType> void TetMesh_t<pointType>::init(uint32_t& unswap_k, uint32_t& unswap_l){
  const uint32_t n = numVertices();

  // Find non-coplanar vertices (we assume that no coincident vertices exist)
  int ori=0;
  uint32_t i=0, j=1, k=2, l=3;

  for (; ori == 0 && k < n - 1; k++)
      for (l = k + 1; ori == 0 && l < n; l++)
          ori = vOrient3D(i, j, k, l);

  l--; k--;

  if (ori == 0) {
      assert(0 && "TetMesh_t::init() - Input vertices do not define a volume");
      ip_error("TetMesh_t::init() - Input vertices do not define a volume.\n");
  }

  unswap_k = k;
  unswap_l = l;
  std::swap(vertices[k], vertices[2]); k=2;
  std::swap(vertices[l], vertices[3]); l=3;

  if(ori<0) std::swap(i, j); // Tets must have positive volume

  const uint32_t base_tet[] = { l, k, j, i, l, j, k, INFINITE_VERTEX, l, k, i, INFINITE_VERTEX, l, i, j, INFINITE_VERTEX, k, j, i, INFINITE_VERTEX };
  const uint64_t base_neigh[] = { 19, 15, 11, 7, 18, 10, 13, 3, 17, 14, 5, 2, 16, 6, 9, 1, 12, 8, 4, 0 };

  resizeTets(5);
  std::memcpy(getTetNodes(0), base_tet, 20 * sizeof(uint32_t));
  std::memcpy(getTetNeighs(0), base_neigh, 20 * sizeof(uint64_t));
  compute_subDet(0);
  compute_subDet(4);
  compute_subDet(8);
  compute_subDet(12);
  compute_subDet(16);

  // set the vertex-(one_of_the)incident-tetrahedron relation
  inc_tet[i] = inc_tet[j] = inc_tet[k] = inc_tet[l] = 0;
}


template<class pointType> void TetMesh_t<pointType>::tetrahedrize() {
    uint32_t uk, ul;
    init(uk, ul); // First tet is made of vertices 0, 1, uk, ul

    // Need to unswap immediately to keep correct indexing and
    // ensure symbolic perturbation is coherent
    if (ul != 3) {
        std::swap(vertices[ul], vertices[3]);
        std::swap(inc_tet[ul], inc_tet[3]);
        for (uint32_t& tn : tet_node) if (tn == 3) tn = ul; else if (tn == ul) tn = 3;
    }

    if (uk != 2) {
        std::swap(vertices[uk], vertices[2]);
        std::swap(inc_tet[uk], inc_tet[2]);
        for (uint32_t& tn : tet_node) if (tn == 2) tn = uk; else if (tn == uk) tn = 2;
    }

    uint64_t ct = 0;
    for (uint32_t i = 2; i < numVertices(); i++) if (i != uk && i != ul) insertExistingVertex(i, ct);

    removeDelTets();
}


template<class pointType> bool TetMesh_t<pointType>::saveTET(const char* filename, bool inner_only) const
{
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh_t::saveTET: Can't open file for writing.\n";
        return false;
    }

    f << numVertices() << " vertices\n";

    uint32_t ngnt = 0;
    for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN) ngnt++;

    if (inner_only) {
        f << ngnt << " tets\n";
        for (uint32_t i = 0; i < numVertices(); i++)
            if constexpr (std::is_pointer_v<pointType>) f << *vertices[i] << "\n";
            else f << vertices[i] << "\n";

        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    }
    else {
        f << ngnt << " inner tets\n";
        f << countNonGhostTets()-ngnt << " outer tets\n";
        for (uint32_t i = 0; i < numVertices(); i++)
            if constexpr (std::is_pointer_v<pointType>) f << *vertices[i] << "\n";
            else f << vertices[i] << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    }
    
    f.close();

    return true;
}

template<class pointType> bool TetMesh_t<pointType>::saveVTU(const char* filename, bool inner_only) const
{
    ofstream f(filename);
    if (!f) {
        std::cerr << "\nTetMesh_t::saveVTU: Can't open file for writing.\n";
        return false;
    }

    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    f << "  <UnstructuredGrid>\n";

    std::vector<uint32_t> tets_to_save;
    if (inner_only) {
        for (uint32_t i = 0; i < numTets(); i++) {
            if (mark_tetrahedra[i] == DT_IN) {
                tets_to_save.push_back(i);
            }
        }
    }
    else {
        for (uint32_t i = 0; i < numTets(); i++) {
            if (!isGhost(i)) {
                tets_to_save.push_back(i);
            }
        }
    }
    size_t num_t = tets_to_save.size();
    uint32_t num_v = numVertices();

    f << "    <Piece NumberOfPoints=\"" << num_v << "\" NumberOfCells=\"" << num_t << "\">\n";
    f << "      <Points>\n";
    f << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    for (uint32_t i = 0; i < num_v; i++) {
        if constexpr (std::is_pointer_v<pointType>) f << "          " << *vertices[i] << "\n";
        else f << "          " << vertices[i] << "\n";
    }
    f << "        </DataArray>\n";
    f << "      </Points>\n";
    f << "      <Cells>\n";
    f << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (uint32_t i : tets_to_save) {
        const uint32_t* nodes = tet_node.data() + i * 4;
        f << "          " << nodes[0] << " " << nodes[1] << " " << nodes[2] << " " << nodes[3] << "\n";
    }
    f << "        </DataArray>\n";
    f << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    f << "          ";
    for (uint32_t i = 1; i <= num_t; i++) {
        f << i * 4 << " ";
    }
    f << "\n        </DataArray>\n";
    f << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    f << "          ";
    for (uint32_t i = 0; i < num_t; i++) {
        f << "10 "; // VTK_TETRA
    }
    f << "\n        </DataArray>\n";
    f << "      </Cells>\n";
    f << "    </Piece>\n";
    f << "  </UnstructuredGrid>\n";
    f << "</VTKFile>\n";

    f.close();
    return true;
}

template<class pointType> bool TetMesh_t<pointType>::saveMEDIT(const char* filename, bool inner_only) const
{
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh_t::saveMEDIT: Can't open file for writing.\n";
        return false;
    }

    f << "MeshVersionFormatted 2\nDimension\n3\n";

    f << "Vertices\n" << numVertices() << "\n";

    uint32_t ngnt = 0;
    for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN) ngnt++;

    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    if (inner_only) {
        for (uint32_t i = 0; i < numVertices(); i++)
            if constexpr (std::is_pointer_v<pointType>) f << *vertices[i] << " 1\n";
            else f << vertices[i] << " 1\n";
        f << "Tetrahedra\n" << ngnt << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << tet_node[i * 4]+1 << " " << tet_node[i * 4 + 2] + 1 << " " << tet_node[i * 4 + 1] + 1 << " " << tet_node[i * 4 + 3] + 1 << " 1\n";
    }
    else {
        for (uint32_t i = 0; i < numVertices(); i++)
            if constexpr (std::is_pointer_v<pointType>) f << *vertices[i] << " 1\n";
            else f << vertices[i] << " 1\n";
        f << "Tetrahedra\n" << countNonGhostTets() << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << tet_node[i * 4] + 1 << " " << tet_node[i * 4 + 2] + 1 << " " << tet_node[i * 4 + 1] + 1 << " " << tet_node[i * 4 + 3] + 1 << " 1\n";
        for (uint32_t i = 0; i < numTets(); i++) if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
            f << tet_node[i * 4] + 1 << " " << tet_node[i * 4 + 2] + 1 << " " << tet_node[i * 4 + 1] + 1 << " " << tet_node[i * 4 + 3] + 1 << " 2\n";
    }

    f.close();

    return true;
}


template<class pointType> bool TetMesh_t<pointType>::saveBinaryTET(const char* filename, bool inner_only) const
{
    ofstream f(filename, ios::binary);

    if (!f) {
        std::cerr << "\nTetMesh_t::saveBinaryTET: Can't open file for writing.\n";
        return false;
    }

    uint32_t num_v = numVertices(), num_t = 0;

    for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN) num_t++;

    f << num_v << " vertices\n";

    if (inner_only) {
        f << num_t << " tets\n";
    }
    else {
        f << num_t << " inner tets\n";
        f << countNonGhostTets() - num_t << " outer tets\n";
    }

    double c[3];
    for (uint32_t i = 0; i < numVertices(); i++) {
        if constexpr (std::is_pointer_v<pointType>) vertices[i]->getApproxXYZCoordinates(c[0], c[1], c[2], true);
        else vertices[i].getApproxXYZCoordinates(c[0], c[1], c[2], true);
        f.write((const char*)(&c), sizeof(double) * 3);
    }

    const uint32_t* tnd = tet_node.data();

    if (inner_only) {
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f.write((const char*)(tnd + i * 4), sizeof(uint32_t) * 4);
    }
    else {
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f.write((const char*)(tnd + i * 4), sizeof(uint32_t) * 4);
        for (uint32_t i = 0; i < numTets(); i++) if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
            f.write((const char*)(tnd + i * 4), sizeof(uint32_t) * 4);
    }

    f.close();

    return true;
}

template<class pointType> bool TetMesh_t<pointType>::saveBoundaryToOFF(const char* filename) const {
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh_t::saveBoundaryToOFF: Can't open file for writing.\n";
        return false;
    }

    f << "OFF\n" << numVertices() << " ";

    size_t num_tris = 0;
    for (uint64_t i = 0; i < tet_node.size(); i++)
        if (i > tet_neigh[i] && mark_tetrahedra[tet_neigh[i] >> 2] != mark_tetrahedra[i >> 2]) num_tris++;

    f << num_tris << " 0\n";

    for (uint32_t i = 0; i < numVertices(); i++)
        if constexpr (std::is_pointer_v<pointType>) f << *vertices[i] << "\n";
        else f << vertices[i] << "\n";

    uint32_t fv[3];
    for (uint64_t i = 0; i < tet_node.size(); i++)
        if (i > tet_neigh[i] && mark_tetrahedra[tet_neigh[i] >> 2] != mark_tetrahedra[i >> 2]) {
            getFaceVertices(i, fv);
            if (fv[0] == INFINITE_VERTEX || fv[1] == INFINITE_VERTEX || fv[2] == INFINITE_VERTEX) ip_error("Attempting to save skin of invalid In/Out classification.\n");
            f << "3 " << fv[0] << " " << fv[1] << " " << fv[2] << "\n";
        }
    f.close();

    return true;
}

template<class pointType> bool TetMesh_t<pointType>::saveRationalTET(const char* filename, bool inner_only)
{
#ifdef USE_INDIRECT_PREDS
    ofstream f(filename);

    if (!f) {
        std::cerr << "\nTetMesh_t::saveRationalTET: Can't open file for writing.\n";
        return false;
    }

    f << numVertices() << " vertices\n";

    uint32_t ngnt = 0;
    for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN) ngnt++;

    if (inner_only) {
        f << ngnt << " tets\n";
        for (uint32_t i = 0; i < numVertices(); i++) {
            bigrational c[3];
            vertices[i]->getExactXYZCoordinates(c[0], c[1], c[2]);
            f << c[0] << " " << c[1] << " " << c[2] << "\n";
        }
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    }
    else {
        f << ngnt << " inner tets\n";
        f << countNonGhostTets() - ngnt << " outer tets\n";
        for (uint32_t i = 0; i < numVertices(); i++) {
            bigrational c[3];
            vertices[i]->getExactXYZCoordinates(c[0], c[1], c[2]);
            f << c[0] << " " << c[1] << " " << c[2] << "\n";
        }
        for (uint32_t i = 0; i < numTets(); i++) if (mark_tetrahedra[i] == DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
        for (uint32_t i = 0; i < numTets(); i++) if (!isGhost(i) && mark_tetrahedra[i] != DT_IN)
            f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    }

    f.close();
#endif

    return true;
}

// Swap t and l while assuming that l is a valid tet (not to be deleted)
template<class pointType> inline void TetMesh_t<pointType>::moveDeletedToTail(uint64_t t, uint64_t l) {
    const uint64_t t2 = t >> 2, l2 = l >> 2;
    uint32_t* tnt = tet_node.data() + t, *lnt = tet_node.data() + l;
    uint64_t* tnn = tet_neigh.data() + t, *lnn = tet_neigh.data() + l;
    const uint32_t* te = tnt + 4;

    do {
        *tnt++ = *lnt;

        uint64_t neigh = *lnn++;
        *tnn++ = neigh;
        tet_neigh[neigh] = t++;

        if (*lnt != INFINITE_VERTEX && inc_tet[*lnt] == l2)
            inc_tet[*lnt] = t2;
        lnt++;
    } while (tnt != te);

    mark_tetrahedra[t2] = mark_tetrahedra[l2];
    if constexpr (std::is_same_v<pointType, basicVec3d>) {
        double* tsd = tet_subdet.data();
        std::memcpy(tsd + t, tsd + l, 4 * sizeof(double));
    }
}

template<class pointType> void TetMesh_t<pointType>::removeDelTets() {
    uint64_t* dp = Del_deleted.data();
    uint64_t* de = dp + Del_deleted.size();
    uint64_t last_valid = numTets() - 1;
    while (isToDeleteSmall(last_valid)) last_valid--;

    while (dp != de) {
        if (*dp < (last_valid<<2)) {
            moveDeletedToTail(*dp, (last_valid << 2));
            last_valid--;
            while (isToDeleteSmall(last_valid)) last_valid--;
        }
        dp++;
    }

    resizeTets(++last_valid);
    Del_deleted.clear();
}

template<class pointType> bool TetMesh_t<pointType>::tetHasVertex(uint64_t t, uint32_t v) const {
    t <<= 2;
    return tet_node[t] == v || tet_node[t + 1] == v || tet_node[t + 2] == v || tet_node[t + 3] == v;
}

template<class pointType> void TetMesh_t<pointType>::oppositeTetEdge(const uint64_t tet, const uint32_t v[2], uint32_t ov[2]) const {
    int i = 0, j = 0;
    while (i < 4) {
        const uint32_t w = tet_node[tet + i];
        if (w != v[0] && w != v[1]) ov[j++] = w;
        i++;
    }
    assert(j == 2);
}

template<class pointType> uint64_t TetMesh_t<pointType>::getCornerFromOppositeTet(uint64_t t, uint64_t n) const {
    t <<= 2;
    for (int i = 0; i < 4; i++)
        if ((tet_neigh[t + i] >> 2) == n)
            return tet_neigh[t + i];
    assert(0);
    return UINT64_MAX;
}

template<class pointType> void TetMesh_t<pointType>::getFaceVertices(uint64_t t, uint32_t v[3]) const {
    uint64_t tv = t & 3;
    const uint32_t* Node = tet_node.data() + (t - tv);
    v[0] = Node[(++tv) & 3];
    v[1] = Node[(++tv) & 3];
    v[2] = Node[(++tv) & 3];
}

template<class pointType> void TetMesh_t<pointType>::getFaceSortedVertices(uint64_t t, uint32_t v[3]) const {
    getFaceVertices(t, v);
    if (v[0] > v[1]) std::swap(v[0], v[1]);
    if (v[1] > v[2]) std::swap(v[1], v[2]);
    if (v[0] > v[1]) std::swap(v[0], v[1]);
}

template<class pointType> bool TetMesh_t<pointType>::getTetsFromFaceVertices(uint32_t v1, uint32_t v2, uint32_t v3, uint64_t* nt) const {
    static thread_local std::vector<uint64_t> vt; // Static to avoid reallocation at each call
    VTfull(v1, vt);
    int i = 0;
    for (uint64_t t : vt) if (tetHasVertex(t, v2) && tetHasVertex(t, v3)) nt[i++] = t;
    vt.clear();
    return (i == 2);
}

template<class pointType> uint64_t TetMesh_t<pointType>::tetOppositeCorner(uint64_t t, uint32_t v1, uint32_t v2, uint32_t v3) const {
    const uint64_t tb = t << 2;
    const uint32_t* n = tet_node.data() + tb;
    for (int i = 0; i < 3; i++)
        if (n[i] != v1 && n[i] != v2 && n[i] != v3)
            return tet_neigh[tb + i];
    assert(n[3] != v1 && n[3] != v2 && n[3] != v3);
    return tet_neigh[tb + 3];
}

template<class pointType> void TetMesh_t<pointType>::resizeTets(uint64_t new_size) {
    mark_tetrahedra.resize(new_size, 0);
    new_size <<= 2;
    tet_node.resize(new_size);
    tet_neigh.resize(new_size);
    if constexpr (std::is_same_v<pointType, basicVec3d>) tet_subdet.resize(new_size);
}

template<class pointType> void TetMesh_t<pointType>::reserveTets(uint64_t new_capacity) {
    mark_tetrahedra.reserve(new_capacity);
    new_capacity <<= 2;
    tet_node.reserve(new_capacity);
    tet_neigh.reserve(new_capacity);
    if constexpr (std::is_same_v<pointType, basicVec3d>) tet_subdet.reserve(new_capacity);
}

template<class pointType> uint64_t TetMesh_t<pointType>::searchTetrahedron(uint64_t tet, const uint32_t v_id)
{
    if (tet_node[tet + 3] == INFINITE_VERTEX)
        tet = getIthNeighbor(getTetNeighs(tet), 3);

    uint64_t i, f0 = 4;
    do {
        const uint32_t* Node = getTetNodes(tet);
        if (Node[3] == INFINITE_VERTEX) return tet;

        const uint64_t* Neigh = getTetNeighs(tet);
        for (i = 0; i < 4; i++)
            if (i != f0 && vOrient3D(Node[tetON1(i)], Node[tetON2(i)], Node[tetON3(i)], v_id) < 0) {
                tet = getIthNeighbor(Neigh, i);
                f0 = Neigh[i] & 3;
                break;
            }
    } while (i != 4);

    return tet;
}

template<class pointType> int TetMesh_t<pointType>::symbolicPerturbation(uint32_t indices[5]) const {
    int swaps = 0;
    int n = 5;
    int count;
    do {
        count = 0;
        n--;
        for (int i = 0; i < n; i++) {
            if (indices[i] > indices[i + 1]) {
                std::swap(indices[i], indices[i + 1]);
                count++;
            }
        }
        swaps += count;
    } while (count);

    n = vOrient3D(indices[1], indices[2], indices[3], indices[4]);
    if (n) return (swaps % 2) ? (-n) : n;

    n = vOrient3D(indices[0], indices[2], indices[3], indices[4]);
    return (swaps % 2) ? (n) : (-n);
}

template<class pointType> int TetMesh_t<pointType>::vertexInTetSphere(const uint32_t Node[4], uint32_t v_id) const {
    int det = vInSphere(Node[0], Node[1], Node[2], Node[3], v_id);
    if (det) return det;
    uint32_t nn[5] = { Node[0],Node[1],Node[2],Node[3],v_id };
    det = symbolicPerturbation(nn);
    if (det == 0.0) {
        //std::cout << *vertices[Node[0]] << " (ID: " << Node[0] << ")\n";
        //std::cout << *vertices[Node[1]] << " (ID: " << Node[1] << ")\n";
        //std::cout << *vertices[Node[2]] << " (ID: " << Node[2] << ")\n";
        //std::cout << *vertices[Node[3]] << " (ID: " << Node[3] << ")\n";
        //std::cout << *vertices[v_id] << " (ID: " << v_id << ")\n";
        assert(0 && "Symbolic perturbation failed! Should not happen");
        ip_error("Symbolic perturbation failed! Should not happen.\n");
    }
    return det;
}

template<class pointType> int TetMesh_t<pointType>::vertexInTetSphere(uint64_t tet, uint32_t v_id) const
{
    const uint32_t* Node = getTetNodes(tet);

    if constexpr (std::is_same_v<pointType, basicVec3d>) {
        const double* SubDet = tet_subdet.data() + tet;

        const double* a = vertices[Node[0]].coord;
        const double* e = vertices[v_id].coord;

        const double aex = e[0] - a[0];
        const double aey = e[1] - a[1];
        const double aez = e[2] - a[2];

        if (Node[3] == INFINITE_VERTEX) {
            double det = aex * SubDet[0] + aey * SubDet[1] + aez * SubDet[2];
            if (fabs(det) > o3d_static_filter) {
                return (det > 0) - (det < 0);
            }

            const double* b = vertices[Node[1]].coord;
            const double* c = vertices[Node[2]].coord;

            det = -orient3d(a, b, c, e);
            if (det != 0.0) {
                return (det > 0) - (det < 0);
            }

            const uint32_t nnode[4] = { Node[0], Node[1], Node[2], tet_node[tet_neigh[tet + 3]] };
            return -vertexInTetSphere(nnode, v_id);
        }

        const double aer = aex * aex + aey * aey + aez * aez;
        double det = aex * SubDet[0] - aey * SubDet[1] + aez * SubDet[2] - aer * SubDet[3];
        if (fabs(det) > isp_static_filter) {
            return (det > 0) - (det < 0);
        }

        return vertexInTetSphere(Node, v_id);
    }
    else {
        if (Node[3] != INFINITE_VERTEX) return vertexInTetSphere(Node, v_id);
        else {
            const int det = vOrient3D(Node[0], Node[1], Node[2], v_id);
            if (det != 0) return (det > 0) - (det < 0);
            const uint32_t nn[4] = { Node[1], Node[0], Node[2], tet_node[tet_neigh[tet + 3]] };
            return vertexInTetSphere(nn, v_id);
        }
    }
}

#ifdef USE_MAROTS_METHOD
template<class pointType> void TetMesh_t<pointType>::deleteInSphereTets(uint64_t tet, const uint32_t v_id)
{
  pushAndMarkDeletedTets(tet);

  for(uint64_t t = Del_deleted.size() - 1; t < Del_deleted.size(); t++) {
    uint64_t tet = Del_deleted[t];
    uint64_t* Neigh = getTetNeighs(tet);
    uint32_t* Node = getTetNodes(tet);

    uint64_t neigh = getIthNeighbor(Neigh, 0);
    if(!isToDelete(neigh)){
      if(vertexInTetSphere(neigh, v_id)<0) bnd_push(v_id, Node[1], Node[2], Node[3], Neigh[0]);
      else pushAndMarkDeletedTets(neigh);
    }

    neigh = getIthNeighbor(Neigh, 1);
    if(!isToDelete(neigh)){
      if(vertexInTetSphere(neigh, v_id)<0) bnd_push(v_id, Node[2], Node[0], Node[3], Neigh[1]);
      else pushAndMarkDeletedTets(neigh);
    }

    neigh = getIthNeighbor(Neigh, 2);
    if(!isToDelete(neigh)){
      if(vertexInTetSphere(neigh, v_id)<0) bnd_push(v_id, Node[0], Node[1], Node[3], Neigh[2]);
      else pushAndMarkDeletedTets(neigh);
    }

    neigh = getIthNeighbor(Neigh, 3);
    if(!isToDelete(neigh)){
      if(vertexInTetSphere(neigh, v_id)<0){
        if(Node[1]<Node[2])
          bnd_push(v_id, Node[0], Node[2], Node[1], Neigh[3]);
        else
          bnd_push(v_id, Node[1], Node[0], Node[2], Neigh[3]);
      }
      else pushAndMarkDeletedTets(neigh);
    }
  }
}


template<class pointType> void TetMesh_t<pointType>::tetrahedrizeHole(uint64_t* tet){
  uint64_t clength = Del_deleted.size(); // Num tets removed
  uint64_t blength = numDelTmp(); // Num tets to insert

  uint64_t tn = numTets();

  if(blength > clength){
    for (uint64_t i = clength; i<blength; i++, tn++)
        Del_deleted.push_back(tn<<2);

    clength = blength;
    resizeTets(tn);
  }

  uint64_t start = clength - blength;

  for (uint64_t i=0; i<blength; i++)
  {
    const uint64_t tet = Del_deleted[i + start];
    uint32_t* Node = getTetNodes(tet);

    Node[0] = Del_tmp[i].node[0];
    Node[1] = Del_tmp[i].node[1];
    Node[2] = Del_tmp[i].node[2];
    Node[3] = Del_tmp[i].node[3];

    uint64_t bnd = Del_tmp[i].bnd;
    tet_neigh[tet] = bnd;
    tet_neigh[bnd] = tet;
    Del_tmp[i].bnd = tet;

    mark_tetrahedra[tet >> 2] = 0;
    compute_subDet(tet);

    if(tet_node[tet+3]!=INFINITE_VERTEX)
      for(uint32_t j=0; j<4; j++)
          inc_tet[tet_node[tet + j]] = tet>>2;
  }

  uint64_t tlength = 0;
  const uint64_t middle = blength * 3 / 2;

  uint64_t* Tmp = delTmpVec();
  const unsigned index[4] = { 2,3,1,2 };

  for (uint64_t i = 0; i < blength; i++)
  {
      uint64_t tet = Del_deleted[start + i];
      const uint32_t* Node = getTetNodes(tet);

      for (uint64_t j = 0; j < 3; j++)
      {
          uint64_t key = ((uint64_t)Node[index[j]] << 32) + Node[index[j + 1]];
          tet++;

          uint64_t k;
          for (k = 0; k < tlength; k++) if (Tmp[k] == key) break;

          if (k == tlength) {
              Tmp[tlength] = (key >> 32) + (key << 32);
              Tmp[middle + tlength] = tet;
              tlength++;
          }
          else {
              uint64_t pairValue = Tmp[middle + k];
              tet_neigh[tet] = pairValue;
              tet_neigh[pairValue] = tet;
              tlength--;
              if (k < tlength) {
                  Tmp[k] = Tmp[tlength];
                  Tmp[middle + k] = Tmp[middle + tlength];
              }
          }
      }
  }

  flushDelTmp();
  *tet = Del_deleted[start];
  Del_deleted.resize(start);
}

template<class pointType> void TetMesh_t<pointType>::insertExistingVertex(const uint32_t vi, uint64_t& ct)
{
    ct = searchTetrahedron(ct, vi);
    deleteInSphereTets(ct, vi);
    tetrahedrizeHole(&ct);
    uint64_t lt = ct;
    if (tet_node[lt + 3] == INFINITE_VERTEX) lt = tet_neigh[lt + 3];
    inc_tet[vi] = lt >> 2;
}

template<class pointType> void TetMesh_t<pointType>::insertExistingVertexNonDelaunay(const uint32_t vi, uint64_t ct)
{
    deleteInVertexTets(ct, vi);
    tetrahedrizeHole(&ct);
    uint64_t lt = ct;
    if (tet_node[lt + 3] == INFINITE_VERTEX) lt = tet_neigh[lt + 3];
    inc_tet[vi] = lt >> 2;
}

template<class pointType> void TetMesh_t<pointType>::deleteInVertexTets(uint64_t tet, const uint32_t v_id)
{
    pushAndMarkDeletedTets(tet);

    for (uint64_t t = Del_deleted.size() - 1; t < Del_deleted.size(); t++) {
        uint64_t tet = Del_deleted[t];
        uint64_t* Neigh = getTetNeighs(tet);
        uint32_t* Node = getTetNodes(tet);

        uint64_t neigh = getIthNeighbor(Neigh, 0);
        if (!isToDelete(neigh)) {
            if (vOrient3D(Node[1], Node[2], Node[3], v_id) != 0) {
                bnd_push(v_id, Node[1], Node[2], Node[3], Neigh[0]);
            }
            else pushAndMarkDeletedTets(neigh);
        }

        neigh = getIthNeighbor(Neigh, 1);
        if (!isToDelete(neigh)) {
            if (vOrient3D(Node[2], Node[0], Node[3], v_id) != 0) {
                bnd_push(v_id, Node[2], Node[0], Node[3], Neigh[1]);
            }
            else pushAndMarkDeletedTets(neigh);
        }

        neigh = getIthNeighbor(Neigh, 2);
        if (!isToDelete(neigh)) {
            if (vOrient3D(Node[0], Node[1], Node[3], v_id) != 0) {
                bnd_push(v_id, Node[0], Node[1], Node[3], Neigh[2]);
            }
            else pushAndMarkDeletedTets(neigh);
        }

        neigh = getIthNeighbor(Neigh, 3);
        if (!isToDelete(neigh)) {
            if (vOrient3D(Node[0], Node[2], Node[1], v_id) != 0) {
                if (Node[1] < Node[2])
                    bnd_push(v_id, Node[0], Node[2], Node[1], Neigh[3]);
                else
                    bnd_push(v_id, Node[1], Node[0], Node[2], Neigh[3]);
            }
            else pushAndMarkDeletedTets(neigh);
        }
    }
}

#else

// Expand by adjacencies to collect all tets whose circumsphere contains v_id
template<class pointType> void TetMesh_t<pointType>::getUnconstrainedCavity(const uint32_t v_id, const uint64_t tet, std::vector<uint64_t>& cavityCorners) {
    uint64_t ntet = searchTetrahedron(tet, v_id) >> 2;

    const uint64_t* tet_neigh_data = tet_neigh.data();

    size_t first = Del_deleted.size();
    pushAndMarkDeletedTets(ntet << 2);

    for (size_t i = first; i < Del_deleted.size(); i++) {
        const uint64_t* nb = tet_neigh_data + Del_deleted[i];

        ntet = (*nb) & (~3);
        if (!isToDelete(ntet)) {
            if (vertexInTetSphere(ntet, v_id) < 0) cavityCorners.push_back(*nb);
            else pushAndMarkDeletedTets(ntet);
        }

        ntet = (*(++nb)) & (~3);
        if (!isToDelete(ntet)) {
            if (vertexInTetSphere(ntet, v_id) < 0) cavityCorners.push_back(*nb);
            else pushAndMarkDeletedTets(ntet);
        }

        ntet = (*(++nb)) & (~3);
        if (!isToDelete(ntet)) {
            if (vertexInTetSphere(ntet, v_id) < 0) cavityCorners.push_back(*nb);
            else pushAndMarkDeletedTets(ntet);
        }

        ntet = (*(++nb)) & (~3);
        if (!isToDelete(ntet)) {
            if (vertexInTetSphere(ntet, v_id) < 0) cavityCorners.push_back(*nb);
            else pushAndMarkDeletedTets(ntet);
        }
    }
}

template<class pointType> void TetMesh_t<pointType>::undeleteTets(size_t num_deltets) {
    for (size_t i = Del_deleted.size() - num_deltets; i < Del_deleted.size(); i++)
        mark_tetrahedra[Del_deleted[i] >> 2] &= (~((uint32_t)1073741824)); // No longer to delete
    Del_deleted.resize(Del_deleted.size() - num_deltets);
}

template<class pointType> void TetMesh_t<pointType>::getCavityConnectivity(const std::vector<uint64_t>& cavityCorners, std::vector<uint64_t>& adjs) const {
    const uint32_t* tet_node_data = tet_node.data();
    const uint64_t* tet_neigh_data = tet_neigh.data();
    const uint32_t* tn;
    uint32_t v1, v2, v3, v[5];
    uint64_t b, c, c0, i, j = 0;

    adjs.resize(cavityCorners.size() * 3);

    for (const uint64_t t : cavityCorners) {
        c0 = tet_neigh_data[t];
        b = c0 & 3;
        tn = tet_node_data + c0 - b;
        v[0] = v[3] = tn[tetON1(b)];
        v[1] = v[4] = tn[tetON2(b)];
        v[2] = tn[tetON3(b)];

        for (i = 0; i < 3; i++) {
            v1 = v[i]; v2 = v[i + 1];
            if (v1 < v2) {
                c = c0;
                do {
                    v3 = tet_node_data[c];
                    c &= (~3);
                    tn = tet_node_data + c;
                    while (*tn == v1 || *tn == v2 || *tn == v3) tn++;
                    c = tet_neigh_data[tn - tet_node_data];
                } while (isToDelete(c));
                adjs[j++] = c;
                adjs[j++] = t;
            }
        }
    }
    assert(adjs.size() == j);
}

template<class pointType> void TetMesh_t<pointType>::retetrahedrizeCavity(const uint32_t v_id, const std::vector<uint64_t>& cavityCorners, const std::vector<uint64_t>& adjs, const uint32_t mark) {
    static const int fi[4][3] = { {2, 1, 3} ,{0, 2, 3} ,{1, 0, 3} ,{0, 1, 2} };
    uint32_t v1;
    uint32_t* tet_node_data = tet_node.data();
    uint64_t* tet_neigh_data = tet_neigh.data();

    // Resize the mesh to host the new tets
    uint64_t ntb, newpos = tet_node.size();
    if (cavityCorners.size() > Del_deleted.size()) {
        resizeTets(numTets() + (cavityCorners.size() - Del_deleted.size()));
        tet_node_data = tet_node.data();
        tet_neigh_data = tet_neigh.data();
    }

    // Create the new tets
    for (const uint64_t c : cavityCorners) {
        unmarkToDelete(c);
        if (Del_deleted.empty()) {
            ntb = newpos;
            newpos += 4;
        }
        else {
            ntb = Del_deleted.back();
            Del_deleted.pop_back();
        }
        const uint64_t cb = c & 3;
        const uint32_t* cr = tet_node_data + (c - cb);
        uint32_t* cn = tet_node_data + ntb;

        assert(cr[fi[cb][2]] == INFINITE_VERTEX || vOrient3D(v_id, cr[fi[cb][0]], cr[fi[cb][1]], cr[fi[cb][2]]) > 0);

        *cn++ = v_id;
        *cn++ = cr[fi[cb][0]];
        *cn++ = cr[fi[cb][1]];
        *cn++ = cr[fi[cb][2]];

        tet_neigh_data[ntb] = c; tet_neigh_data[c] = ntb;

        compute_subDet(ntb);

        ntb >>= 2;
        if ((*(--cn)) != INFINITE_VERTEX) {
            inc_tet[*cn] = ntb;
            inc_tet[*(--cn)] = ntb;
            inc_tet[*(--cn)] = ntb;
            inc_tet[v_id] = ntb;
        }
        mark_tetrahedra[ntb] = mark;
    }

    // Restore inner connectivity
    for (size_t i = 0; i < adjs.size(); ) {
        uint64_t c = tet_neigh_data[adjs[i++]] & (~3); c++;
        uint64_t n = tet_neigh_data[adjs[i++]] & (~3); n++;
        uint64_t k = c;

        v1 = tet_node_data[c];
        if (v1 == tet_node_data[n + 1] || v1 == tet_node_data[n + 2] || v1 == tet_node_data[n]) v1 = tet_node_data[++c];
        if (v1 == tet_node_data[n + 1] || v1 == tet_node_data[n + 2] || v1 == tet_node_data[n]) ++c;
        v1 = tet_node_data[n];
        if (v1 == tet_node_data[k + 1] || v1 == tet_node_data[k + 2] || v1 == tet_node_data[k]) v1 = tet_node_data[++n];
        if (v1 == tet_node_data[k + 1] || v1 == tet_node_data[k + 2] || v1 == tet_node_data[k]) ++n;

        tet_neigh_data[c] = n; tet_neigh_data[n] = c;
    }
}

// Collect all tets whose circumsphere contains v_id and replace them
// with a star of new tets originating at v_id

template<class pointType> void TetMesh_t<pointType>::insertExistingVertex(const uint32_t v_id, uint64_t& tet)
{
    static thread_local std::vector<uint64_t> cavityCorners, adjs; // Static to avoid reallocation on each call
    getUnconstrainedCavity(v_id, tet, cavityCorners);
    getCavityConnectivity(cavityCorners, adjs);
    retetrahedrizeCavity(v_id, cavityCorners, adjs, 0);

    tet = tet_neigh[cavityCorners.back()];

    cavityCorners.clear();
}

template<class pointType> void TetMesh_t<pointType>::insertExistingVertexNonDelaunay(const uint32_t v_id, uint64_t ntet)
{
    static thread_local std::vector<uint64_t> cavityCorners, adjs;

    ntet &= (~3ULL);

    assert(tet_node[ntet + 3] != INFINITE_VERTEX);

    const uint64_t* tet_neigh_data = tet_neigh.data();

    size_t first = Del_deleted.size();
    pushAndMarkDeletedTets(ntet);

    // fi[k] lists the three nodes of the tet face opposite node k, in the same order
    // retetrahedrizeCavity uses to build the new (v_id, face) tet.
    static const int fi[4][3] = { {2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2} };
    const uint32_t* tet_node_data = tet_node.data(); // pushAndMarkDeletedTets does not resize tet_node

    for (size_t i = first; i < Del_deleted.size(); i++) {
        const uint64_t* nb = tet_neigh_data + Del_deleted[i];

        for (int f = 0; f < 4; f++) {
            const uint64_t nbv = nb[f];
            const uint64_t nt = nbv & (~3ULL);
            if (isToDelete(nt)) continue;

            // Absorb the neighbour into the cavity iff v_id lies on the face SHARED with it
            // -- the one retetrahedrizeCavity builds the new (v_id, face) tet from, i.e. the
            // neighbour's face opposite node (nbv & 3). The previous code always tested the
            // neighbour's {1,2,3} face (opposite node 0); when the shared face was a different
            // one, a shared face containing v_id could be kept as a cavity boundary, producing
            // a zero-volume (v_id, face) tet -- the assert in retetrahedrizeCavity (debug), or
            // an endless loop (release). A ghost neighbour (node 3 == INFINITE_VERTEX) is
            // never absorbed.
            const uint64_t cb = nbv & 3;
            const uint32_t* cr = tet_node_data + nt;
            if (cr[3] == INFINITE_VERTEX
                || vOrient3D(cr[fi[cb][0]], cr[fi[cb][1]], cr[fi[cb][2]], v_id) != 0)
                cavityCorners.push_back(nbv);
            else
                pushAndMarkDeletedTets(nt);
        }
    }

    getCavityConnectivity(cavityCorners, adjs);
    retetrahedrizeCavity(v_id, cavityCorners, adjs, 0);

    cavityCorners.clear();
}

#endif

template<class pointType> void TetMesh_t<pointType>::VT(uint32_t v, std::vector<uint64_t>& vt) const {
    static thread_local std::vector<uint64_t> vt_queue; // Static to avoid reallocation at each call
    uint64_t t = inc_tet[v];

    vt_queue.push_back(tetCornerAtVertex(t << 2, v));
    mark_Tet_31(t);

    for (size_t i = 0; i < vt_queue.size(); i++) {
        t = vt_queue[i];
        const uint64_t sb = t & 3;
        const uint64_t* tg = tet_neigh.data() + t - sb;
        for (int j = 1; j < 4; j++) {
            const uint64_t tb = tg[(sb+j)&3];
            const uint64_t tbb = tb >> 2;
            if (tet_node[tb] != INFINITE_VERTEX && !is_marked_Tet_31(tbb)) {
                const uint64_t nt = tetCornerAtVertex(tb & (~3), v);
                vt_queue.push_back(nt); 
                mark_Tet_31(tbb); 
            }
        }
    }

    for (uint64_t t : vt_queue) {
        t >>= 2;
        unmark_Tet_31(t);
        vt.push_back(t);
    }

    vt_queue.clear();
}

template<class pointType> void TetMesh_t<pointType>::VV(uint32_t v, std::vector<uint32_t>& vv) const {
    static thread_local std::vector<uint64_t> vt_queue; // Static to avoid reallocation at each call
    uint64_t t = inc_tet[v];
    const uint64_t tb = t << 2;

    const uint64_t s = tetCornerAtVertex(tb, v);
    vt_queue.push_back(s);
    mark_Tet_31(t);

    const uint32_t* tn = tet_node.data() + tb;
    const uint64_t sb = s & 3;
    for (int j = 1; j < 4; j++) {
        const uint32_t w = tn[(sb + j) & 3];
        marked_vertex[w] |= 128;
        vv.push_back(w);
    }

    for (size_t i = 0; i < vt_queue.size(); i++) {
        t = vt_queue[i];
        const uint64_t sb = t & 3;
        const uint64_t* tg = tet_neigh.data() + t - sb;
        for (int j = 1; j < 4; j++) {
            const uint64_t tb = tg[(sb + j) & 3];
            const uint64_t tbb = tb >> 2;
            const uint32_t w = tet_node[tb];
            if (w != INFINITE_VERTEX && !is_marked_Tet_31(tbb)) {
                vt_queue.push_back(tetCornerAtVertex(tb & (~3), v));
                mark_Tet_31(tbb);
                if (!(marked_vertex[w] & 128)) {
                    marked_vertex[w] |= 128;
                    vv.push_back(w);
                }
            }
        }
    }

    for (uint64_t t : vt_queue) unmark_Tet_31(t>>2);
    vt_queue.clear();
    for (uint32_t w : vv) marked_vertex[w] &= 127;
}

template<class pointType> void TetMesh_t<pointType>::ET(uint32_t v1, uint32_t v2, std::vector<uint64_t>& et) const {
    VT(v1, et);
    for (size_t i = 0; i < et.size();)
        if (!tetHasVertex(et[i], v2)) {
            std::swap(et[i], et[et.size() - 1]);
            et.pop_back();
        }
        else i++;
}

template<class pointType> void TetMesh_t<pointType>::ETfull(uint32_t v1, uint32_t v2, std::vector<uint64_t>& et) const {
    VTfull(v1, et);
    for (size_t i = 0; i < et.size();)
        if (!tetHasVertex(et[i], v2)) {
            std::swap(et[i], et[et.size() - 1]);
            et.pop_back();
        }
        else i++;
}

template<class pointType> void TetMesh_t<pointType>::ETfast(uint32_t v1, uint32_t v2, uint64_t st_t, std::vector<uint64_t>& et) const {
    assert(tetHasVertex(st_t, v1) && tetHasVertex(st_t, v2));

    uint64_t prev_t = UINT64_MAX, next_t, cur_t = st_t;
    do {
        const uint64_t tb = cur_t * 4;
        if (tet_node[tb + 3] != INFINITE_VERTEX) et.push_back(cur_t);
        for (int i = 0; i < 4; i++) {
            next_t = tet_neigh[tb + i] >> 2;
            if (next_t != prev_t && tetHasVertex(next_t, v1) && tetHasVertex(next_t, v2)) break;
        }
        prev_t = cur_t;
        cur_t = next_t;
    } while (cur_t != st_t);
}

template<class pointType> void TetMesh_t<pointType>::ETfullfast(uint32_t v1, uint32_t v2, uint64_t st_t, std::vector<uint64_t>& et) const {
    assert(tetHasVertex(st_t, v1) && tetHasVertex(st_t, v2));

    uint64_t prev_t = UINT64_MAX, next_t, cur_t = st_t;
    do {
        const uint64_t tb = cur_t * 4;
        et.push_back(cur_t);
        for (int i = 0; i < 4; i++) {
            next_t = tet_neigh[tb + i] >> 2;
            if (next_t != prev_t && tetHasVertex(next_t, v1) && tetHasVertex(next_t, v2)) break;
        }
        prev_t = cur_t;
        cur_t = next_t;
    } while (cur_t != st_t);
}

template<class pointType> void TetMesh_t<pointType>::ETcorners(uint32_t v1, uint32_t v2, std::vector<uint64_t>& et) const {
    uint64_t t;
    VTfull(v1, et);
    for (uint64_t s : et) if (tetHasVertex(s, v2)) { t = (s<<2); break; }
    
    while (tet_node[t] == v1 || tet_node[t] == v2) t++;

    et.clear();

    uint64_t c0 = t;
    do {
        et.push_back(t); // Add tet
        uint64_t oc = tet_neigh[t] & (~3); // Get next base
        uint32_t cv = tet_node[t];
        t &= (~3);
        while (tet_node[t] == v1 || tet_node[t] == v2 || tet_node[t] == cv) t++;
        t = tetCornerAtVertex(oc, tet_node[t]); // Get corresp corner at opposite tet
    } while (t != c0);
}

template<class pointType> void TetMesh_t<pointType>::VTfull(uint32_t v, std::vector<uint64_t>& vt) const {
    static thread_local std::vector<uint64_t> vt_queue; // Static to avoid reallocation at each call
    uint64_t s, t = inc_tet[v];
    vt_queue.push_back(t);
    mark_Tet_31(t);

    while (!vt_queue.empty()) {
        t = vt_queue.back();
        vt_queue.pop_back();
        vt.push_back(t);
        t <<= 2;
        s = tet_neigh[t] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) { vt_queue.push_back(s); mark_Tet_31(s); }
        s = tet_neigh[t + 1] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) { vt_queue.push_back(s); mark_Tet_31(s); }
        s = tet_neigh[t + 2] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) { vt_queue.push_back(s); mark_Tet_31(s); }
        s = tet_neigh[t + 3] >> 2;
        if (!is_marked_Tet_31(s) && tetHasVertex(s, v)) { vt_queue.push_back(s); mark_Tet_31(s); }
    }

    for (uint64_t t : vt) unmark_Tet_31(t);
}


template<class pointType> bool TetMesh_t<pointType>::hasEdge(uint32_t v1, uint32_t v2) const {
    static thread_local std::vector<uint64_t> vt_queue; // Static to avoid reallocation at each call

    uint64_t t = inc_tet[v1];
    const uint64_t tb = t << 2;
    if (tet_node[tb] == v2 || tet_node[tb + 1] == v2 || tet_node[tb + 2] == v2 || tet_node[tb + 3] == v2) return true;

    vt_queue.push_back(tetCornerAtVertex(tb, v1));
    mark_Tet_31(t);

    for (size_t i = 0; i < vt_queue.size(); i++) {
        t = vt_queue[i];
        const uint64_t sb = t & 3;
        const uint64_t* tg = tet_neigh.data() + t - sb;
        for (int j = 1; j < 4; j++) {
            const uint64_t tb = tg[(sb + j) & 3];
            const uint64_t tbb = tb >> 2;
            const uint32_t w = tet_node[tb];
            if (w != INFINITE_VERTEX && !is_marked_Tet_31(tbb)) {
                vt_queue.push_back(tetCornerAtVertex(tbb << 2, v1));
                mark_Tet_31(tbb);
                if (w == v2) {
                    for (uint64_t t : vt_queue) unmark_Tet_31(t >> 2);
                    vt_queue.clear();
                    return true;
                }
            }
        }
    }

    for (uint64_t t : vt_queue) unmark_Tet_31(t >> 2);
    vt_queue.clear();
    return false;
}


template<class pointType> void TetMesh_t<pointType>::swapTets(const uint64_t t1, const uint64_t t2) 
{
    if (t1 == t2) return;

    const uint64_t t1_id = t1<<2;
    const uint64_t t2_id = t2<<2;

    // update VT base relation
    for (int i = 0; i < 3; i++) if (inc_tet[tet_node[t1_id + i]] == t1) inc_tet[tet_node[t1_id + i]] = t2;
    if (tet_node[t1_id + 3] != INFINITE_VERTEX && inc_tet[tet_node[t1_id + 3]] == t1) inc_tet[tet_node[t1_id + 3]] = t2;

    for (int i = 0; i < 3; i++) if (inc_tet[tet_node[t2_id + i]] == t2) inc_tet[tet_node[t2_id + i]] = t1;
    if (tet_node[t2_id + 3] != INFINITE_VERTEX && inc_tet[tet_node[t2_id + 3]] == t2) inc_tet[tet_node[t2_id + 3]] = t1;

    // Update nodes and marks
    for (int i = 0; i < 4; i++) std::swap(tet_node[t1_id + i], tet_node[t2_id + i]);
    std::swap(mark_tetrahedra[t1], mark_tetrahedra[t2]);
    if constexpr (std::is_same_v<pointType, basicVec3d>) {
        double* tsd = tet_subdet.data();
        std::swap(tsd[t1_id], tsd[t2_id]);
        std::swap(tsd[t1_id+1], tsd[t2_id+1]);
        std::swap(tsd[t1_id+2], tsd[t2_id+2]);
        std::swap(tsd[t1_id+3], tsd[t2_id+3]);
    }

    // update neigh-neigh relations
    const uint64_t ng1[] = { tet_neigh[t1_id + 0], tet_neigh[t1_id + 1], tet_neigh[t1_id + 2], tet_neigh[t1_id + 3] };
    const uint64_t ng2[] = { tet_neigh[t2_id + 0], tet_neigh[t2_id + 1], tet_neigh[t2_id + 2], tet_neigh[t2_id + 3] };

    for (int i = 0; i < 4; i++) if ((ng2[i] >> 2) != t1) tet_neigh[ng2[i]] = t1_id + i;
    for (int i = 0; i < 4; i++) if ((ng1[i] >> 2) != t2) tet_neigh[ng1[i]] = t2_id + i;

    for (int i = 0; i < 4; i++)
        if ((ng2[i] >> 2) != t1) tet_neigh[t1_id + i] = tet_neigh[t2_id + i];
        else tet_neigh[t1_id + i] = (tet_neigh[t2_id + i] & 3) + (t2 << 2);

    for (int i = 0; i < 4; i++)
        if ((ng1[i] >> 2) != t2) tet_neigh[t2_id + i] = ng1[i];
        else tet_neigh[t2_id + i] = (ng1[i] & 3) + (t1 << 2);
}

template<class pointType> bool TetMesh_t<pointType>::hasBadSnappedOrientations(size_t& num_flipped, size_t& num_flattened) const {
    const uint32_t* tn = tet_node.data();
    const uint32_t* end = tn + tet_node.size();
    num_flipped = num_flattened = 0;
    double v[12];
    while (tn < end) {
        if (tn[3] != INFINITE_VERTEX) {
            for (int i = 0; i < 4; i++) {
                const pointType& p = vertices[tn[i]];
                if constexpr (std::is_pointer_v<pointType>) p->getApproxXYZCoordinates(v[i * 3], v[i * 3 + 1], v[i * 3 + 2], true);
                else p.getApproxXYZCoordinates(v[i * 3], v[i * 3 + 1], v[i * 3 + 2], true);
            }
            const int o = orient3d(v, v+3, v+6, v+9);
            if (o > 0) num_flipped++;
            else if (o == 0) num_flattened++;
        }
        tn += 4;
    }

    return (num_flipped || num_flattened);
}

template<class pointType> void TetMesh_t<pointType>::compute_subDet(const uint64_t tet)
{
    if constexpr (std::is_same_v<pointType, basicVec3d>) {
        double ab[4], ac[4];

        uint32_t* Node = tet_node.data() + tet;
        double* SubDet = tet_subdet.data() + tet;

        double* a = vertices[Node[0]].coord;
        double* b = vertices[Node[1]].coord;
        double* c = vertices[Node[2]].coord;

        if (Node[3] != INFINITE_VERTEX) {
            double* d = vertices[Node[3]].coord;
            double ad[4];

            unsigned i;
            for (i = 0; i < 3; i++) {
                ab[i] = b[i] - a[i];
                ac[i] = c[i] - a[i];
                ad[i] = d[i] - a[i];
            }

            ab[3] = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
            ac[3] = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
            ad[3] = ad[0] * ad[0] + ad[1] * ad[1] + ad[2] * ad[2];

            double cd12 = ac[2] * ad[1] - ac[1] * ad[2];
            double db12 = ad[2] * ab[1] - ad[1] * ab[2];
            double bc12 = ab[2] * ac[1] - ab[1] * ac[2];

            double cd30 = ac[0] * ad[3] - ac[3] * ad[0];
            double db30 = ad[0] * ab[3] - ad[3] * ab[0];
            double bc30 = ab[0] * ac[3] - ab[3] * ac[0];

            SubDet[0] = ab[3] * cd12 + ac[3] * db12 + ad[3] * bc12;
            SubDet[1] = ab[2] * cd30 + ac[2] * db30 + ad[2] * bc30;
            SubDet[2] = ab[1] * cd30 + ac[1] * db30 + ad[1] * bc30;
            SubDet[3] = ab[0] * cd12 + ac[0] * db12 + ad[0] * bc12;
        }
        else {
            unsigned i;
            for (i = 0; i < 3; i++) {
                ab[i] = b[i] - a[i];
                ac[i] = c[i] - a[i];
            }

            SubDet[0] = ac[1] * ab[2] - ac[2] * ab[1];
            SubDet[1] = ac[2] * ab[0] - ac[0] * ab[2];
            SubDet[2] = ac[0] * ab[1] - ac[1] * ab[0];
            SubDet[3] = SubDet[0] * SubDet[0] + SubDet[1] * SubDet[1] + SubDet[2] * SubDet[2];
        }
    }
}

template<class pointType> void TetMesh_t<pointType>::checkMesh(bool checkDelaunay) const {
    size_t i;
    const uint32_t num_vertices = (uint32_t)vertices.size();
    // Check tet nodes	
    for (i = 0; i < numTets(); i++) if (!isToDelete(i<<2)) {
        const uint32_t* tn = tet_node.data() + i * 4;
        if (tn[0] >= num_vertices) assert(0 && "Wrong tet node!\n");
        if (tn[1] >= num_vertices) assert(0 && "Wrong tet node!\n");
        if (tn[2] >= num_vertices) assert(0 && "Wrong tet node!\n");
        if (tn[3] != INFINITE_VERTEX && tet_node[i * 4 + 3] >= num_vertices) assert(0 && "Wrong tet node!\n");
        if (tn[0] == tn[1] || tn[0] == tn[2] || tn[0] == tn[3]
            || tn[1] == tn[2] || tn[1] == tn[3] || tn[2] == tn[3]) 
            assert(0 && "Wrong tet node indexes!\n");
    }

    // Check neighbors	
    for (i = 0; i < numTets() * 4; i++) if (!isToDelete(i))
        if (tet_neigh[i] >= tet_neigh.size() || tet_neigh[tet_neigh[i]] != i)
            assert(0 && "Wrong neighbor!\n");

    // Check neighbor-node coherence
    for (i = 0; i < numTets() * 4; i++) if (!isToDelete(i)) {
        if (tetHasVertex(tet_neigh[i] >> 2, tet_node[i]))
            assert(0 && "Incoherent neighbor!\n");
        else {
            uint32_t v[3];
            getFaceVertices(i, v);
            if (!tetHasVertex(tet_neigh[i] >> 2, v[0])) assert(0 && "Incoherent face at neighbors!\n");
            if (!tetHasVertex(tet_neigh[i] >> 2, v[1])) assert(0 && "Incoherent face at neighbors!\n");
            if (!tetHasVertex(tet_neigh[i] >> 2, v[2])) assert(0 && "Incoherent face at neighbors!\n");
        }
    }

    // Check vt*	
    for (i = 0; i < num_vertices; i++) if (inc_tet[i]!=UINT64_MAX) {
        if (inc_tet[i] >= numTets())
            assert(0 && "Wrong vt* (out of range)!\n");
        if (isGhost(inc_tet[i]))
            assert(0 && "Wrong vt* (ghost tet)!\n");
        if (isToDeleteSmall(inc_tet[i]))
            assert(0 && "Wrong vt* (deleted tet)!\n");
        const uint32_t* tn = tet_node.data() + inc_tet[i] * 4;
        if (tn[0] != i && tn[1] != i && tn[2] != i && tn[3] != i)
            assert(0 && "Wrong vt*!\n");
    }

    // Check marks
    //for (i = 0; i < numTets(); i++) if (!isToDelete(i<<2))
    //    if (mark_tetrahedra[i])
    //        assert(0 && "Marked tet\n");

    // Check geometry
    for (i = 0; i < numTets(); i++) if (!isToDelete(i<<2)) {
        const uint32_t* tn = tet_node.data() + i * 4;
        if (tn[3] != INFINITE_VERTEX && vOrient3D(tn[0], tn[1], tn[2], tn[3]) <= 0) assert(0 && "Inverted/degn tet\n");
    }

    if (checkDelaunay) {
        for (size_t i = 0; i < numTets(); i++) if (!isToDelete(i<<2)) {
            const uint32_t* n = tet_node.data() + (i * 4);
            if (n[3] == INFINITE_VERTEX) continue;
            for (int j = 0; j < 4; j++) {
                uint64_t oppc = tet_neigh[i * 4 + j];
                uint32_t ov = tet_node[oppc];
                if (ov != INFINITE_VERTEX && ((mark_tetrahedra[i] & 3) == (mark_tetrahedra[oppc>>2] & 3)) && vertexInTetSphere(n, ov) > 0) assert(0 && "Non delaunay\n");
            }
        }
    }

    printf("checkMesh passed\n");
}
