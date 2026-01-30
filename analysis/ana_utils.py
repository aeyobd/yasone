import shapely
import numpy as np


def is_in_poly(x, y, x_poly, y_poly):
    N = len(x)
    poly = shapely.Polygon(np.vstack([x_poly, y_poly]).T)
    
    in_poly = np.full(N, False)
    for i in range(N):
        if np.isfinite(x[i]) and np.isfinite(y[i]):
            point = shapely.Point(x[i], y[i])
            in_poly[i] = point.within(poly)

    return in_poly
