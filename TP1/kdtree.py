# coding: utf8
#import openalea.plantgl.all as pgl
import numpy as np 
from numpy.linalg import norm

class KDNode:
    def __init__(self, pivot = None, axis = None, left_child = None, right_child = None):
        self.pivot       = pivot
        self.axis        = axis
        self.left_child  = left_child
        self.right_child = right_child
 
 

def create_kdtree(point_list, minbucketsize=3, depth=0):
    if (len(point_list) > 2*minbucketsize):
        axis = depth % 3
        point_list.sort(key=lambda point_list: point_list[axis])
        med = len(point_list) // 2
        n = KDNode()
        n.pivot = point_list[med]
        n.axis = axis
        n.left_child = create_kdtree(point_list[:med],minbucketsize,depth+1)
        n.right_child = create_kdtree(point_list[med+1:],minbucketsize,depth+1)
        return n
    else:
        return point_list
 
def generate_random_point(size=[1,1,1], distribution='uniform'):
    from random import uniform, gauss
    if distribution == 'uniform':
        return np.array([uniform(-size[0],size[0]), uniform(-size[1],size[1]), uniform(-size[2],size[2])]) 
    elif distribution == 'gaussian':
        return np.array([gauss(0,size[0]/3.), gauss(0,size[1]/3.), gauss(0,size[1]/3.)]) 
 
 
def generate_random_pointlist(size=[1,1,1], nb = 100, distribution='uniform'):
    return [generate_random_point(size, distribution=distribution) for i in xrange(nb)]
 
def generate_random_kdtree(nbpoints=100, size=[1,1,1], minbucketsize=2):
    points = generate_random_pointlist(nb = nbpoints, size=size, distribution='uniform')
    mkdtree = create_kdtree(points, minbucketsize)
    view_kdtree(mkdtree)
    return mkdtree

def brute_force_closest(point, pointlist):
    """ Find the closest points of 'point' in 'pointlist' using a brute force approach """
    import sys
    pid, d = -1, sys.maxint
    for i, p in enumerate(pointlist):
        nd = norm(point-p) 
        if nd < d:
            d = nd
            pid = i
    return pointlist[pid]
 


def kdtree_representation(ax, kdtree, bbox=[[-1., 1.],[-1., 1.],[-1., 1.]]):
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import numpy as np
    if isinstance(kdtree, KDNode):
        dim = kdtree.axis
        plane_bbox = [b for i,b in enumerate(bbox) if i != dim]
        plane_points = []
        plane_points += [np.insert([plane_bbox[0][0],plane_bbox[1][0]],dim,kdtree.pivot[dim])]
        plane_points += [np.insert([plane_bbox[0][0],plane_bbox[1][1]],dim,kdtree.pivot[dim])]
        plane_points += [np.insert([plane_bbox[0][1],plane_bbox[1][1]],dim,kdtree.pivot[dim])]
        plane_points += [np.insert([plane_bbox[0][1],plane_bbox[1][0]],dim,kdtree.pivot[dim])]
        
        plane_points = [[v.tolist() for v in plane_points]]

        left_bbox = np.copy(bbox).astype(float)
        right_bbox = np.copy(bbox).astype(float)
        left_bbox[dim,1] = kdtree.pivot[dim]
        right_bbox[dim,0] = kdtree.pivot[dim]

        poly = Poly3DCollection( plane_points )

        poly.set_edgecolor(list(np.insert([0,0],dim,1)))
        poly.set_facecolors([0,0,0,0.2])
        ax.add_collection3d(poly)
        ax.scatter([kdtree.pivot[0]], [kdtree.pivot[1]], [kdtree.pivot[2]], c=(180/255.,150/255.,50/255.) )
        kdtree_representation(ax,kdtree.left_child, bbox=left_bbox)
        kdtree_representation(ax,kdtree.right_child, bbox=right_bbox)

    else:
        assert (type(kdtree) == list) or (isinstance(kdtree,np.ndarray))
        kdtree = np.array(kdtree)
        ax.scatter(kdtree[:,0],kdtree[:,1],kdtree[:,2], c=(150/255.,150/255.,150/255.) )



def view_kdtree(kdtree, bbox=[[-1., 1.],[-1., 1.],[-1., 1.]]):
    import  matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    kdtree_representation(ax, kdtree, bbox)
    plt.show()


 
def kdtree_closest(point, kdtree):
    axis = kdtree.axis
    if point < kdtree.pivot:
        best = kdtree_closest(point,kdtree.left_child)
        opp = kdtree.right_child
    else:
        best = kdtree_closest(point,kdtree.right_child)
        opp = kdtree.left_child
    # De la même manière qu’a la création, 
    # on détermine la subdivision dans laquelle le point appartient, 
    # on trouve le point le plus proche dans la subdivision 
    # (appelé meilleur candidat)
    # On remonte dans la hiérarchie en testant à chaque niveau,
    # Si un meilleur candidat n’est pas dans la subdivision opposée ou au pivot :
    # On teste cela en comparant la distance de p au meilleur candidat courant et 
    # la distance de p au plan de subdivision (donné par la coordonnée n du pivot)
        
    # Si la distance est trop grande:
        # On teste si le pivot ne fournit pas un meilleur candidat.
        # On teste les subdivisions opposées avec la même procédure
    # Si la distance au candidat est plus petite, on remonte
    return None
 
def test_kdtree(nbtest=100, nbpoints=100, size=[1,1,1], minbucketsize=2):
    import time
 
    points = generate_random_pointlist(nb = nbpoints, size=size, distribution='uniform')
    mkdtree = create_kdtree(points, minbucketsize)
    #view_kdtree(mkdtree, radius=0.03, bbox=[[-float(s),float(s)] for s in size])
    kdtime, bftime = 0,0
    for i in xrange(nbtest):
        testpoint = generate_random_point()
        t = time.time()
        kpoint = kdtree_closest(testpoint, mkdtree)
        kdtime += time.time()-t
        t = time.time()
        bfpoint = brute_force_closest(testpoint, points)
        bftime += time.time()-t
        if norm(kpoint - bfpoint) > 1e-5: 
            raise ValueError('Invalid closest point')
    if nbtest == 1:
        view_kdtree(mkdtree, bbox=[[-float(s),float(s)] for s in size])
    print 'Comparative execution time : KD-Tree [', kdtime,'], BruteForce [', bftime,']'
 
    return kdtime, bftime
 
 
def plot_execution_time(nbpoints_min=10, nbpoints_max=5000):
    import matplotlib.pyplot as plt
    
    # Réaliser une boucle sur le nombre de points
    # Stocker les temps d'exécution renvoyés par test_kdtree dans deux listes
    
    # Visualiser les courbes à l'aide de la fonction plot (http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot)
    plt.figure("Execution Time")
    #plt.plot(...,color='r',label="KD-Tree")
    #plt.plot(...,color='b',label="Brute Force")
    plt.legend()
    plt.show()

if __name__ == '__main__':
    #test_kdtree(nbtest=1, nbpoints= 100)
    generate_random_kdtree()