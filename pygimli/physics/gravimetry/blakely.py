import numpy as np

def magneticDipole(Q, M, P=None, x=None, y=0., z=0., alpha=0, cylinder=False):
    """Compute magnetic field.
    
    Parameters
    ----------
    Q : array 3x1
        position of magnetic dipole
    M : array 3x1   
        magnetization vector
    P : array 3xN
        measuring positions
    x : array N
        positions as profile
    y/z : array | float
        y and z positions
    alpha : profile direction (degree)
        profile direction
    cylinder : bool [False]
        use line/cylinder instead of point/sphere
    """
    if x is not None:
        if isinstance(z, (int, float)):
            z = np.ones_like(x) * z
        if isinstance(y, (int, float)):
            y = np.ones_like(x) * y
        
        P = np.column_stack([x * np.cos(np.deg2rad(alpha)), 
                             y + x * np.sin(np.deg2rad(alpha)), z])
    
    P -= Q
    r = np.sqrt(np.sum(P**2, axis=1))  # distance as scalar
    M0 = np.linalg.norm(M)
    M = np.array(M, dtype=float) / M0  # unit vector
    R = P / np.reshape(r, [-1, 1])  # norm vectors
    my0 = 4 * np.pi * 1e-7
    if cylinder:
        fak = np.reshape(my0 * M0 / 2 / np.pi / r**2, [-1, 1])
        return (np.reshape(R.dot(M), [-1, 1]) * R * 2 - M) * fak
    else:
        fak = np.reshape(my0 * M0 / 4 / np.pi / r**3, [-1, 1])
        return (np.reshape(R.dot(M), [-1, 1]) * R * 3 - M) * fak
