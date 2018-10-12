import numpy as np
import scipy.ndimage as nd
import math
from scipy.cluster.vq import vq

from gmap_tools import read_ply_mesh, array_unique
from array_dict import array_dict

class Point:
	def __init__(x,y,z):
        self.x = x
        self.y = y
        self.z = z




def display(self, color = [190,205,205], add = False):
    from openalea.plantgl.all import Scene, Shape, Material, FaceSet, Viewer
        from random import randint
        s = Scene()
        for facedart in self.elements(2):
                lastdart = facedart
                positions = []
                for dart in self.orderedorbit(facedart,[0,1]):
                        if self.alpha(0, dart) != lastdart:
                                positions.append(self.get_position(dart))
                        lastdart = dart
                if color is None:
                        mat = Material((randint(0,255),randint(0,255),randint(0,255)))
                else:
                        mat = Material(tuple(color),diffuse=0.25)
                s.add(Shape(FaceSet(positions, [range(len(positions))]) , mat, facedart ))
        if add : 
            Viewer.add(s)
        else : 
            Viewer.display(s)

def bernstein(i, n, u):
    """ Bernstein Polynomials
   Initial conditions : $B_0^0 = 1$, $B_{-1}^n = B_{n+1}^n = 0$
   $B_i^{n+1}(u) =uB_{i-1}^n(u) + (1-u)B_{i}^n(u)$ 
    """
    if i == 0 and n == 0:
    	return 1;

    if n == -1 or i == n+1:
    	return 0;

    f1 = u*bernstein(i-1,n,u)
    f2 = (1-u)*bernstein(i,n,u)
    return f1 + f2

def spline(control_points, nb_points=10):
    """ Evaluate the spline curve defined by its control points. 
    return a list of nb_points.
    """
	res = []
	if nb_points == 0:
		return res

	for i in xrange(nb_points):
		u = i/nb_points
		sumx = 0
		sumy = 0
		sumz = 0
		for point in control_points:
			sumx += bernstein(i,n,u) * point.x
			sumy += bernstein(i,n,u) * point.y
			sumz += bernstein(i,n,u) * point.z
		res.append(Point(sumx,sumy,sumz))
	return res

def plot_spline_crv(ctrls, pts):
    """ 
    Parameters
    ==========
      - ctrl: control points
      - pts : evaluated points on the curve
    """
    scene = Scene()
    crv = Shape(geometry=Polyline(pts), appearance=Material((125,12,12)))
    scene.add(crv)

    ctrl = Shape(geometry=Polyline(ctrls), appearance=Material((125,12,12)))
    scene.add(ctrl)
    # To complete: Draw the control points and the line between each ones.

    Viewer.display(scene)


if __name__ == '__main__':
    P1 = Point(0,0,0)
    P2 = Point(5,5,0)
    P3 = Point(10,5,0)
    P4 = Point(15,0,0)

    ctrls = []
    ctrls.append(P1)
    ctrls.append(P2)
    ctrls.append(P3)
    ctrls.append(P4)

    res = spline(ctrls,5)
    plot_spline_crv(ctrls,res)
