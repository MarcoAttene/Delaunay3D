// Regression test for insertExistingVertexNonDelaunay when the inserted vertex lies
// exactly on an edge shared by several tets.
//
// The cube [-1,2]^3 is tetrahedralized into 6 tets that all share the main diagonal edge
// 0--7 (from (-1,-1,-1) to (2,2,2)). Vertices are then inserted with the *non-Delaunay*
// insertion; one of them, (0,0,1), ends up lying exactly on an edge shared by several tets
// of the current mesh.
//
// Before the fix, the cavity was grown by testing the neighbour's {1,2,3} face (opposite
// node 0) instead of the face actually SHARED with it (opposite node nbv&3). A shared face
// that the new vertex lies on was then kept as a cavity boundary, so retetrahedrizeCavity
// built a zero-volume (v_id, face) tet. In a debug build that trips the orientation assert
// in retetrahedrizeCavity; in a release build (assert compiled out) the routine never
// terminates.
//
// With the fix the insertions complete and every resulting finite tet is non-degenerate.
// Exit code 0 = pass, 1 = failure.

#include "delaunay.h"

#include <cstdint>
#include <cstdio>
#include <vector>

int main()
{
    const double cube[8 * 3] = {
        -1, -1, -1,   2, -1, -1,   -1, 2, -1,   2, 2, -1,
        -1, -1,  2,   2, -1,  2,   -1, 2,  2,   2, 2,  2,
    };
    // 6 tets, all sharing edge 0--7 (the (-1,-1,-1)->(2,2,2) diagonal).
    const uint32_t cube_tets[6 * 4] = {
        0, 1, 3, 7,   0, 5, 1, 7,   0, 3, 2, 7,
        0, 2, 6, 7,   0, 4, 5, 7,   0, 6, 4, 7,
    };

    TetMesh_t<basicVec3d> mesh;
    mesh.init_vertices(cube, 8);
    mesh.init_tets(cube_tets, 6);

    // Insert four extra vertices. (0,0,0) is on the shared diagonal edge; after it is
    // inserted, (0,0,1) lands on an edge shared by several of the new tets -- the trigger.
    const double extra[4 * 3] = {
        0, 0, 0,
        0, 0, 1,
        0, 1, 0.5,
        1, 0, 0.5,
    };
    const uint32_t first_extra = mesh.numVertices();
    for (int i = 0; i < 4; i++)
        mesh.pushVertex(basicVec3d(extra[3 * i], extra[3 * i + 1], extra[3 * i + 2]));

    uint64_t tet = 0;
    for (uint32_t vi = first_extra; vi < mesh.numVertices(); vi++) {
        tet = mesh.searchTetrahedron(tet, vi);
        if (tet == UINT64_MAX) {
            printf("FAIL: vertex %u not located in any tet\n", vi);
            return 1;
        }
        mesh.insertExistingVertexNonDelaunay(vi, tet); // asserts (debug) / loops forever (release) before the fix
    }

    // Every finite (non-ghost) tet must be non-degenerate and consistently oriented: no
    // zero-volume sliver (the bug's symptom) and no flipped tet. Rather than assume the
    // mesh's orientation convention, require all finite tets to share the sign of the first.
    int degenerate = 0, flipped = 0, ref = 0;
    for (uint32_t t = 0; t < mesh.numTets(); t++) {
        if (mesh.isGhost(t)) continue;
        const uint32_t* n = mesh.tet_node.data() + 4 * t;
        const int s = basicVec3d::orient3D(mesh.vertices[n[0]], mesh.vertices[n[1]],
                                           mesh.vertices[n[2]], mesh.vertices[n[3]]);
        if (s == 0) { degenerate++; continue; }
        if (ref == 0) ref = s;
        else if (s != ref) flipped++;
    }
    if (degenerate || flipped) {
        printf("FAIL: %d degenerate and %d flipped finite tet(s)\n", degenerate, flipped);
        return 1;
    }

    printf("OK: %u vertices inserted on shared edges, all finite tets positively oriented\n",
           mesh.numVertices());
    return 0;
}
