import math
from openalea.plantgl.all import *

#Calcul polynôme Bernstein

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

#Calcul casteljau

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

#Fonction qui calcul n points sur une courbe à partir des points de controles 

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
		res.append([sumx,sumy,sumz])
	return res

#Affichage des courbes

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

"""
def surface_bezier(ctrls,ctrls2,n,m):
    res = []
    controles = []

    for i in range(n) :
	   u = i/float(n -1)

	   points = []
	   for id in range(len(ctrls2)):
	       newP = [ctrls2[id][0]+vect[0],ctrls2[id][1]+vect[1],ctrls2[id][2]+vect[2]]
	       points.append(newP)
	   res.append(points)

    return res
"""
"""
def real_bezier_sufarce(ctrl1,n,m):
    res = []
    for list in ctrl1:
            res.append(casteljau(ctrl1,n))
    for index in range(len(ctrl1)):
        points = []
        for list in ctrl1:
            points.append(list[index])
        res.append(casteljau(points))
    return res
"""
#nbCourbeHori = n
#nbCourbeVerti = m
def real_bezier_sufarce(ctrl1,n,m):
    res = []

    listeCourbesHori = []
    listeCourbesVerti = []
    for i in range(n):
    	listeCourbesHori.append(spline(ctrl1[i], m))
 	
 	listPtsCtrCourbeVerti = [[] for i in range(m)]

 	for k in range(m):
 		for elmts in listeCourbesHori:
 			listPtsCtrCourbeVerti[k].append(elmts[k])

    for j in range(m):
    	listeCourbesVerti.append(spline(listPtsCtrCourbeVerti[j], m))
    res = [ctrl1,listPtsCtrCourbeVerti]

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
"""
def basis(i, k, u, knots):
    i: ith basis
        k: degree
        u : parameter in [u0,u(n+k+1)]
        knots : knot vector
   
    if k == 0 :
        return 1

    nu = (i/k)

    if u[i] <= nu and u[i+1] > nu:
        return 0

    f1 = (nu-u[i])/(u[i]+k-u[i])*basis(i,k-1,u,knots)
    f2 = ((u[i+k+1]-nu)/(u[i+k+1])-(u[i+1]))*basis(i+1,k-1,u,knots)

    return f1 + f2
"""
def basis(i, k, u, knots):
    """ i: ith basis
        k: degree
        u : parameter in [u0,u(n+k+1)]
        knots : knot vector
    """
    if(k == 0):
    	if (knots[i] <= u) and (u < knots[i+1]):
    		return 1
    	else: return 0
    else:
    	return (((u-knots[i])/(knots[i+k]-knots[i])) * basis(i,k-1,u,knots)) + (((knots[i+k+1]-u)/(knots[i+k+1]-knots[i+1])) * basis(i+1,k-1,u,knots))



def knot_vector(k, n, u_min=0., u_max=1.):
    """ Uniform knot vector.
    """
    m = k+n+2
    knots = [u_min]*(k+1)
    n_internals = m-2*k-1

    pas = (u_max - u_min) / n
    value = u_min
    print(len(knots))
    for i in range(n-1):
    	value += pas
    	knots.append(value)
	print(len(knots),m)
 	for i in range(k):
 		knots.append(u_max)
    # complete
    #uniforme knot vector means, have k same values at each end (beggining and end)  of the vector
    # ti+1 - ti = constant pour tout i
 
 	

    assert(len(knots) == m)
    return knots

def spline_b(control_points, nb_points=10):
	v = knot_vector(len(control_points)-1,len(control_points)-1)
	res = []
	if nb_points == 0:
		return res
	for i in range(nb_points):
		u = i/float(nb_points-1)
		sumx = 0
		sumy = 0
		sumz = 0
		for id in range(len(control_points)):
			sumx += basis(id,len(control_points)-1,u,v) * control_points[id][0]
			sumy += basis(id,len(control_points)-1,u,v) * control_points[id][1]
			sumz += basis(id,len(control_points)-1,u,v) * control_points[id][2]
		res.append([sumx,sumy,sumz])
	return res

def rcalc(i,k,u,knots,control_points_weight):
    f1 = basis(i,k,u,knots)*control_point_weight[i]
    sum = 0
    for j in range(len(control_points_weight)):
        sum += basis(j,k,u,knots)*control_point_weight[j]

    if sum == 0:
        sum = 1

    return f1/sum


def nurbs(control_points,control_points_weight, nb_points=10):
    flag = True
    for i in control_points_weight:
        if (i != 1):
            flag = false

    if (flag):
        result = spline_b(control_points,nb_points)
        return result

    v = knot_vector(len(control_points)-1,len(control_points)-1)
    res = []
    if nb_points == 0:
        return res

    for i in range(nb_points):
        u = i/float(nb_points-1)
        sumx = 0
        sumy = 0
        sumz = 0
        for id in range(len(control_points)):
            sumx += rcalc(id,len(control_points)-1,u,v) * control_points_weight[id]
            sumy += basis(id,len(control_points)-1,u,v) * control_points_weight[id]
            sumz += basis(id,len(control_points)-1,u,v) * control_points_weight[id]
        res.append([sumx,sumy,sumz])
    return res

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

    res = spline_b(ctrls,4)
    #res2 = casteljau(ctrls2,4)

    plot_spline_crv(ctrls, res)

    P1 = [0.0,0.0,0.0]
    P2 = [5.0,0.0,10.0]
    P3 = [10.0,0.0,0.0]
    P4 = [0.0,5.0,0.0]
    P5 = [5.0,5.0,0.0]
    P6 = [10.0,5.0,0.0]
    P7 = [0.0,10,0.0]
    P8 = [5.0,10.0,10.0]
    P9 = [10.0,10.0,0.0]
    ptC1 =[P1,P2,P3]
    ptC2 =[P4,P5,P6]
    ptC3 =[P7,P8,P9]
    listPtsCtrCourbeHori = [ptC1,ptC2,ptC3]

    listsCourbes = real_bezier_sufarce(listPtsCtrCourbeHori,3,10)
    pts = listsCourbes[1]

    ctrl_net = listsCourbes[0]
    #plot_spline_surface(ctrl_net, pts)
