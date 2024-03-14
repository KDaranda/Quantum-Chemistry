##################################
# ATOM CONNECTIVITY
# Determines bond structure,
# calculates structural parameters
##################################
"""Connectivity module.

Provides functions to create connectivity data (internal coordinates)
from the Cartesian (XYZ) coordinate set. It's possible to convert data
back to the Cartesian coordinates given 3 fixed coords for orientation.

The following data types are used for connectivity info:
coordinates - {  m: (x, y, z) }
distances   - { (m, n): r }
bonds       - { (m, n): r }
angles      - { (m, n, p): a }
order       - [ (m, n, p) ]
dihedrals   - { (m, n, p, q): d }
bondmap     - {  m: [n, ...] }
"""


import math
import math_vector as vec


def get_coords (struct):
    """Read in coordinates from input data structure."""
    return { num: xyz for (num, sym, xyz) in struct }


def get_distances (coords):
    """Calculate distances between points in coordinate set."""
    distances = {}
    for m in sorted(coords):
        for n in sorted(coords):
            if m < n:
                distances[(m, n)] = vec.distance(coords[m], coords[n])
    return distances


# GUESSWORK
RCUT = 1.65     # Includes all carbon bonds
FSIM = 1.12     # Bond similarity factor
RMAX = 2.50     # Maximum bond length


def get_bonds (coords, mode='R+C', separateRedun=False):
    """Gather and return bonds between points.

    OPTIONAL ARGUMENTS
        'mode' selects connectivity type:
    "MINI" - minimal connectivity (single bond for each atom),
    "RCUT" - MINI + all short bonds (determined by RCUT),
    "CLUS" - MINI + bonds of similar length (determined by FSIM),
    "R+C"  - RCUT and CLUS together [default].
    (Redundant modes are useful for chemical consideration.)
        'separateRedun' controls output of redundant bonds:
    False - all bonds are returned together [default],
    True  - primary and redundant bonds are returned separately.
    """
    # Set up data
    N = len(coords)         # Number of atoms
    trees = set()           # Connected atom groups
    bonds, redun = {}, {}   # Primary and redundant bonds
    r_max = 0               # Longest primary bond
    # Gather bonds
    distances = get_distances(coords)
    for bond in sorted(distances, key=distances.get):
        (m, n), r = bond, distances[bond]
        # Stop if outer radius is reached
        if r > RMAX:
            break
        # Atom groups containing this bond
        joined = set()
        for tree in trees:
            if m in tree or n in tree:
                joined.add(tree)
        # CONNECTIVITY:
        if len(joined) == 0:
            # Non-connected bond (new atom group)
            bonds[bond] = r_max = r
            trees.add(frozenset(bond))
        elif len(joined) == 2:
            # Bond joins two atom groups
            bonds[bond] = r_max = r
            trees = trees.difference(joined)
            trees.add(frozenset( atom for tree in joined for atom in tree ))
        elif len(joined) == 1:
            # Bond connects to single atom group
            tree = joined.pop()
            if not (m in tree and n in tree):
                # Add primary bond
                bonds[bond] = r_max = r
                trees.remove(tree)
                trees.add(tree.union(bond))
            else:
                # Minimal mode: skip redundant bonds
                if mode == "MINI" and len(tree) < N:
                    pass
                # Cut-off mode: add short redundant bonds
                elif (mode == "RCUT" or mode == "R+C") and r <= RCUT:
                    redun[bond] = r
                # Cluster mode: add bonds similar to the primary
                elif (mode == "CLUS" or mode == "R+C"):
                    # Look up shortest primary bonds for the current atoms
                    ordered_bonds = sorted(bonds, key=bonds.get)
                    rm = next(( bonds[b] for b in ordered_bonds if m in b ), 0)
                    rn = next(( bonds[b] for b in ordered_bonds if n in b ), 0)
                    # Add bond if it's similar to the primary bonds
                    if r <= max(rm, rn) * FSIM:
                        redun[bond] = r
                # Stop if none of above criteria are met
                else:
                    break
    # Return redundant bonds separately if requested
    if separateRedun:
        return (bonds, redun)
    bonds.update(redun)
    return bonds


def get_graph (bonds, directed=False):
    """Return connection (bonding) graph.

    OPTIONAL ARGUMENTS
        'directed' controls graph directionality:
    False - connections to lower numbers are included [default],
    True  - connections to lower numbers are excluded.
    """
    graph = {}
    for (m, n) in sorted(bonds):
        graph[m] = graph.get(m, []) + [n]
        if not directed:
            graph[n] = graph.get(n, []) + [m]
    return graph


def get_angles (bonds, coords=None, naturalOrder=False):
    """Return angles defined between provided bonds.

    OPTIONAL ARGUMENTS
        'coords' is used to provide coordinate data
    for calculation of angle values [default: None].
        'naturalOrder' controls output of key ordering:
    False - no custom order is provided [default],
    True  - list of keys, ordered by their appearance, is returned.
    """
    # Get directed bonding graph
    graph = get_graph(bonds, directed=True)
    # Gather angles
    angles = {}
    if naturalOrder:
        order = []
    for m in sorted(graph):
        for n in graph[m]:
            # Forward angles
            for q in graph.get(n, []):
                angles[(m, n, q)] = None
                if naturalOrder:
                    order += [(m, n, q)]
            # Corner angles
            for p in graph[m]:
                if n < p:
                    angles[(n, m, p)] = None
                    if naturalOrder:
                        order += [(n, m, p)]
            # Backward angles
            for k in sorted(graph):
                if n in graph[k] and m < k and k < n:
                    angles[(m, n, k)] = None
                    if naturalOrder:
                        order += [(m, n, k)]
    # Calculate values if coordinates are provided
    if coords:
        for (m, n, p) in angles:
            angles[(m, n, p)] = vec.angle(coords[m], coords[n], coords[p])
    return (angles, order) if naturalOrder else angles


# Brute-force version;
# flips and messes up natural ordering of dihedrals
def get_diheds (angles, coords=None):
    """Return dihedral angles defined between provided angles.

    OPTIONAL ARGUMENTS
        'coords' is used to provide coordinate data
    for calculation of angle values [default: None].
    """
    # Gather dihedral angles
    diheds = {}
    for (k, m, n) in sorted(angles):
        for (p, q, s) in sorted(angles):
            # Forward dihedrals
            if (m, n) == (p, q) and (s, n, m, k) not in diheds:
                diheds[(k, m, n, s)] = None
            # Front-to-back dihedrals
            elif (m, n) == (s, q) and (p, n, m, k) not in diheds:
                diheds[(k, m, n, p)] = None
            # Back-to-front dihedrals
            elif (k, m) == (q, p) and (n, m, k, s) not in diheds:
                diheds[(s, k, m, n)] = None
    # Calculate values if coordinates are provided
    if coords:
        for (m, n, p, q) in diheds:
            diheds[(m, n, p, q)] = vec.dihed(coords[m], coords[n],
                                             coords[p], coords[q])
    return diheds


# Ordered version;
# requires to store valence angle order
def get_diheds_ordered (ordered_angles, coords=None):
    """Return dihedral angles defined between angles in natural order.

    OPTIONAL ARGUMENTS
        'coords' is used to provide coordinate data
    for calculation of angle values [default: None].
    """
    # Gather dihedral angles
    diheds = {}
    for (i, (k, m, n)) in enumerate(ordered_angles):
        for (p, q, s) in ordered_angles[i+1:]:
            # Forward dihedrals
            if (m, n) == (p, q):
                diheds[(k, m, n, s)] = None
            # Front-to-back dihedrals
            elif (m, n) == (s, q):
                diheds[(k, m, n, p)] = None
            # Back-to-front dihedrals
            elif (k, m) == (q, p):
                diheds[(s, k, m, n)] = None
    # Calculate values if coordinates are provided
    if coords:
        for (m, n, p, q) in diheds:
            diheds[(m, n, p, q)] = vec.dihed(coords[m], coords[n],
                                             coords[p], coords[q])
    return diheds


def get_xyz (A, B, C, rBC, rCD, aBCD, dABCD):
    """Convert internal coordinates of a 3D point (D) to Cartesian.

    Function implements the SN-NERF algorithm.
    """
    # Step 1: calculate relative position of D in reference orientation
    r, theta, phi = rCD, math.pi - vec.to_rad(aBCD), vec.to_rad(dABCD)
    Dx = r * math.cos(theta)
    Dy = r * math.sin(theta) * math.cos(phi)
    Dz = r * math.sin(theta) * math.sin(phi)
    # Step 2: rotate D to the orientation defined by A, B, C
    uBC = vec.scale(1. / rBC, vec.vec(B, C))
    nABC = vec.unit(vec.cross(vec.vec(A, B), uBC))
    (Mx, My, Mz) = (uBC, vec.cross(nABC, uBC), nABC)
    Dr = vec.rotate(zip(*(Mx, My, Mz)), (Dx, Dy, Dz))
    # Step 3: get absolute position of D by adding C
    return vec.add(C, Dr)


def get_carts (bonds, angles, diheds, c_init):
    """Convert internal coordinate set to Cartesians.

    OPTIONAL ARGUMENTS
        'c_init' are Cartesian coordinates of three initial points,
    used to define system position and rotation.
    """
    # Gather atom numbers
    nums = set( n for bond in bonds for n in bond if n not in c_init )
    # Add initial coordinates
    coords = { n: c_init[n] for n in c_init }
    # Calculate all other coordinates
    for num in nums:
        for (m, n, p, q) in diheds:
            if num == q and all( (x in coords) for x in (m, n, p) ):
                A, B, C = coords[m], coords[n], coords[p]
                rBC = bonds[(n, p) if (n, p) in bonds else (p, n)]
                rCD = bonds[(p, q) if (p, q) in bonds else (q, p)]
                aBCD = angles[(n, p, q) if (n, p, q) in angles else (q, p, n)]
                dABCD = diheds[(m, n, p, q)]
                coords[q] = get_xyz(A, B, C, rBC, rCD, aBCD, dABCD)
                break
    return coords

'''
#######
# TESTS
#######

def params_ref (fname):
    """Reference values for test function."""
    data = {
        "test/3bR2.txt": {
            "bonds": {
                (6, 9): 1.0902054125264653, (23, 37): 1.094105844135292,
                (10, 11): 1.227834995354425, (5, 6): 1.3952842315438814,
                (14, 17): 1.530693671967713, (33, 36): 1.1019197307726183,
                (10, 12): 1.2284108064165666, (40, 41): 1.1036549547014232,
                (17, 19): 1.100632451876193, (1, 6): 1.394632781548247,
                (26, 39): 1.0907997850114381, (22, 23): 1.386229912776737,
                (1, 2): 1.3989341783279157, (24, 25): 1.3993689448998072,
                (44, 45): 1.1058010823063975, (29, 31): 1.097655640301183,
                (40, 42): 1.1002235100564794, (4, 13): 1.3561500441086896,
                (4, 5): 1.4112116934131465, (21, 33): 1.5460483211950395,
                (1, 10): 1.467693571259341, (13, 20): 1.4798129505738216,
                (20, 21): 1.5926953967513688, (27, 28): 1.4054285000333526,
                (17, 18): 1.096946225924954, (33, 35): 1.1004208948379706,
                (44, 46): 1.0971797739700635, (25, 26): 1.4046836672322354,
                (29, 30): 1.1013972348094943, (33, 34): 1.1005799616974679,
                (40, 43): 1.103610441640527, (25, 38): 1.093294801302009,
                (5, 14): 1.5089813856963248, (14, 15): 1.1038936741955723,
                (17, 20): 1.5319450826149088, (23, 24): 1.4121307260908953,
                (2, 7): 1.0888097528664042, (44, 47): 1.100654766450407,
                (24, 40): 1.5107663620209448, (26, 27): 1.3940100997313472,
                (14, 16): 1.1012255332764491, (21, 22): 1.5236269497724833,
                (21, 29): 1.53636345364663, (29, 32): 1.1002375443225887,
                (3, 8): 1.0908371093421785, (22, 27): 1.401750519539408,
                (20, 28): 1.4516855819091818, (2, 3): 1.3873117084710993,
                (28, 44): 1.4540307727451987, (3, 4): 1.4082418047029421
            },
        },
    }
    return data[fname] if fname in data else None


import init
#import logging
import read_base as inp
import pretty_print as pretty


def test (level="warning"):
    """Test of connectivity math."""
    if level is not None:
        init.start_log (level)
    # Name of trial data file
    fname = "test/3bR2.txt"
    (struct, title) = inp.fread_xyz (fname)
    # Process data file
    coords = get_coords (struct)
    d_mat = get_distances (coords)
    bonds = get_bonds (coords)
    b_map = get_graph (bonds, True)
    angles, order = get_angles (bonds, coords, naturalOrder=True)
    diheds = get_diheds_ordered (order, coords)
    c_back = get_carts (bonds, angles, diheds,
                        {x: coords[x] for x in (1, 2, 3)})
    # Convert reference and trial data to strings
    coords_ref = pretty.coordinates (coords)
    pretty.dist_matrix (d_mat)
    bonds_try = pretty.param_table (bonds)
    pretty.bondmap (b_map)
    pretty.param_table (angles)
    pretty.param_table (diheds)
    bonds_ref = pretty.param_table (params_ref(fname)["bonds"])
    coords_try = pretty.coordinates (c_back)
    # Fail immediately if parameter counts differ
    B0, A0, D0 = 50, 91, 136
    if len(bonds) != B0 or len(angles) != A0 or len(diheds) != D0:
        print (len(bonds), len(angles), len(diheds), (B0, A0, D0))
        return init.FAIL
    # Fail if strings are different
    for (trial, ref) in ((bonds_try, bonds_ref), (coords_try, coords_ref)):
        if trial != ref:
            print (pretty.mismatch (ref, trial))
            return init.FAIL
    return init.PASS


############
# LOCAL MAIN
############
import time


def main():
    test()

if __name__ == "__main__":
    start = time.perf_counter()
    test()
    stop  = time.perf_counter()
    print ('Test took %.6f ms' % (1000*(stop-start)))
'''