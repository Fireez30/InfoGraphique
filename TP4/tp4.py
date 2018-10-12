import math
from openalea.plantgl.all import *


def bernstein(i, n, u):
    """ Bernstein Polynomials
   Initial conditions : $B_0^0 = 1$, $B_{-1}^n = B_{n+1}^n = 0$
   $B_i^{n+1}(u) =uB_{i-1}^n(u) + (1-u)B_{i}^n(u)$ 
    """
    if i == 0 and n == 0:
    	return 1

    if i < 0 or i > n:
    	return 0



    f1 = u*bernstein(i-1,n-1,u)
    f2 = (1-u)*bernstein(i,n-1,u)
    return f1 + f2

def casteljau(control_points, nb_points):
	points = []
	for i in range(nb_points):
		u = i/float(nb_points-1)
		tmp1 = control_points

		while(len(tmp1) > 1):
			tmp2 = []
			for j in range(len(tmp1)-1):
				Ptmp1 = tmp1[j]
				Ptmp2 = tmp1[j+1]

				tmp = [Ptmp1[0]*(1-u)+Ptmp2[0]*u,Ptmp1[1]*(1-u)+Ptmp2[1]*u,0]
				tmp2.append(tmp)

			tmp1 = tmp2

		points.append(tmp1[0])

	return points

def spline(control_points, nb_points=10):
	res = []
	if nb_points == 0:
		return res

	for i in range(nb_points):
		u = i/float(nb_points-1)
		sumx = 0
		sumy = 0
		sumz = 0
		for id in range(len(control_points)):
			sumx += bernstein(id,len(control_points)-1,u) * control_points[id][0]
			sumy += bernstein(id,len(control_points)-1,u) * control_points[id][1]
			sumz += bernstein(id,len(control_points)-1,u) * control_points[id][2]
		print(sumx)
		print(" ")
		print(sumy)
		print(" ")
		print(sumz)
		res.append([sumx,sumy,sumz])
	return res


def plot_spline_surface(ctrl_net, points):
    """
    Parameters
    ==========
      - ctrl_net : the net of control points (list of list)
      - points : a set of evaluated points (list of list)
    """
    scene = Scene()
    n = len(points)
    m = len(points[0])

    # Compute a mesh (i.e. TriangleSet) for the set of points
    pointList = [pt for rank in points for pt in rank]
    indexList = []

    for i in range(n-1):
        for j in range(m-1):
            ii = i*m+j
            i1 = (ii, ii+1, ii+m)
            i2 = (ii+1,ii+m+1,ii+m)
            indexList.append(i1)
            indexList.append(i2)

    
    surf = Shape(TriangleSet(pointList, indexList), appearance=Material((12,125,12)))
    scene.add(surf)

    
    # plot the control net
    n = len(ctrl_net)
    m = len(ctrl_net[0])
    for pts in ctrl_net:
        crv = Shape(geometry=Polyline(pts), appearance=Material((125,12,12)))
        scene.add(crv)
        for pt in pts:
            scene.add(Shape(Translated(Vector3(pt),Sphere(radius=0.1))))
            
    for i in range(m):
        pts = [ctrl_net[j][i] for j in range(n)]
        crv = Shape(geometry=Polyline(pts), appearance=Material((12,12,125)))
        scene.add(crv)
        
    Viewer.display(scene)


def surface_bezier(ctrls,ctrls2,n,m):
	res = []
	for i in range(n):
		u = i/float(n -1)

		pu = [ctrls[0][0]*(1-u)+ctrls[-1][0]*u,ctrls[0][1]*(1-u)+ctrls[-1][1]*u,ctrls[0][2]*(1-u)+ctrls[-1][2]*u]
		vect = [pu[0]-ctrls2[0][0],pu[1]-ctrls2[0][1],pu[2]-ctrls2[0][2]]

		points = []
		for id in range(len(ctrls2)):
			newP = [ctrls2[id][0]+vect[0],ctrls2[id][1]+vect[1],ctrls2[id][2]+vect[2]]
			points.append(newP)
		res.append(points)

	return res

def plot_spline_crv(ctrls, pts):
    """ 
    Parameters
    ==========
      - ctrl: control points
      - pts : evaluated points on the curve
    """
    from openalea.plantgl.all import Scene, Shape, Material, FaceSet, Viewer, Polyline
    from random import randint
    scene = Scene()
    crv = Shape(geometry=Polyline(pts), appearance=Material((12,12,125)))
    scene.add(crv)

    ctrl = Shape(geometry=Polyline(ctrls), appearance=Material((12,125,12)))
    scene.add(ctrl)
    # To complete: Draw the control points and the line between each ones.

    Viewer.display(scene)


if __name__ == '__main__':
    P1 = [0,0,0]
    P2 = [5,5,0]
    P3 = [10,5,0]
    P4 = [15,0,0]

    ctrls = []
    ctrls.append(P1)
    ctrls.append(P2)
    ctrls.append(P3)
    ctrls.append(P4)

    P21 = [1,1,0]
    P22 = [4,6,0]
    P23 = [8,0,0]
    P24 = [15,6,0]

    ctrls2 = []
    ctrls2.append(P21)
    ctrls2.append(P22)
    ctrls2.append(P23)
    ctrls2.append(P24)

    res = casteljau(ctrls,4)
    res2 = casteljau(ctrls2,4)

    net = surface_bezier(ctrls,ctrls2,4,4)
    plot_spline_surface([ctrls,ctrls2],net)