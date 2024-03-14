##############################
# 3D VECTOR MATH
# Does what it says on the tin
##############################

import math


def add (A, B):
    """Add two 3D vectors."""
    (Ax, Ay, Az) = A
    (Bx, By, Bz) = B
    return (Ax + Bx, Ay + By, Az + Bz)


def sub (A, B):
    """Subtract two 3D vectors."""
    (Ax, Ay, Az) = A
    (Bx, By, Bz) = B
    return (Ax - Bx, Ay - By, Az - Bz)


def vec (A, B):
    """Return 3D vector pointing from A to B."""
    return sub(B, A)


def scale (c, A):
    """Multiply 3D vector by a scalar."""
    (Ax, Ay, Az) = A
    return (c * Ax, c * Ay, c * Az)


def dot (A, B):
    """Return dot product of two 3D vectors."""
    (Ax, Ay, Az) = A
    (Bx, By, Bz) = B
    return Ax * Bx + Ay * By + Az * Bz


def cross (A, B):
    """Return cross product of two 3D vectors."""
    (Ax, Ay, Az) = A
    (Bx, By, Bz) = B
    Rx = Ay * Bz - Az * By
    Ry = Az * Bx - Ax * Bz
    Rz = Ax * By - Ay * Bx
    return (Rx, Ry, Rz)


def length (A):
    """Return length (magnitude) of 3D vector."""
    (Ax, Ay, Az) = A
    return math.sqrt (Ax * Ax + Ay * Ay + Az * Az)


def unit (A):
    """Return normalized (unit) 3D vector."""
    (Ax, Ay, Az), rA = A, length(A)
    return (Ax / rA, Ay / rA, Az / rA)


def rotate (M, A):
    """Apply rotation matrix to 3D vector."""
    (Mx, My, Mz) = M
    return (dot(Mx, A), dot(My, A), dot(Mz, A))


def to_rad (deg):
    """Convert degrees to radians."""
    return deg * (math.pi / 180.)


def to_deg (rad):
    """Convert radians to degrees."""
    return rad * (180. / math.pi)


def distance (A, B):
    """Return distance between two 3D points."""
    return length(vec(A, B))


def angle_r (A, B, C):
    """Return angle between three 3D points (in radians)."""
    U, V, V_ = vec(A, B), vec(B, C), vec(C, B)
    UxV = cross (U, V)
    return math.atan2 (length(UxV), dot(U, V_))


def angle_d (A, B, C):
    """Return angle between three 3D points (in degrees)."""
    return to_deg (angle_r(A, B, C))


def dihed_r (A, B, C, D):
    """Return dihedral angle between four 3D points (in radians)."""
    U, V, Z = vec(A, B), vec(B, C), vec(C, D)
    Usv = scale (length(V), U)
    UxV, VxZ = cross (U, V), cross (V, Z)
    return math.atan2 (dot(Usv, VxZ), dot(UxV, VxZ))


def dihed_d (A, B, C, D):
    """Return dihedral angle between four 3D points (in degrees)."""
    return to_deg (dihed_r(A, B, C, D))


def di (dihed, PI=180.):
    """Normalize dihedral angle to [-PI; PI]."""
    while abs(dihed) > PI:
        dihed = dihed - math.copysign(2*PI, dihed)
    return dihed


# ALIASES FOR COMMON FUNCTIONS
bond = distance
angle = angle_d
dihed = dihed_d

'''
#######
# TESTS
#######

import init
import logging


def test (level=None):
    """Test vector math module."""
    if level is not None:
        init.start_log(level)
    ACC = init.ACC
    # Reference values
    b0, a0, d0 = 1.45168, 112.184, -169.525
    A, B = (0.707782, 0.261779, 7.643724), (1.077923, -0.359432, 6.384962)
    C, D = (2.267849, 0.341463, 5.721879), (2.512836, -0.128903, 4.285996)
    # Calculated values
    b = bond (A, B)
    a = angle(A, B, C)
    d = dihed(A, B, C, D)
    # Test result accuracy
    EPS = 5.0 * 10**(-ACC)
    logging.info ("Desired accuracy - %.1e", EPS)
    for (trial, ref) in ((b, b0), (a, a0), (d, d0)):
        delta = (trial - ref) / ref
        # Output formatting
        fstring = "%" + str(ACC+8) + "." + str(ACC+1) + "e "
        fstring += fstring
        fstring += "%9.2e"
        logging.debug (fstring, ref, trial, delta)
        # Fail if accuracy not reached
        if delta >= EPS:
            return init.FAIL
    return init.PASS


############
# LOCAL MAIN
############
def main ():
    test()

if __name__ == "__main__":
    main()
'''