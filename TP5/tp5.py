# coding: utf8

import numpy as np

class GMap:
    def __init__(self, degree=2):
        """ 
        Constructor 
        """

        self.maxid = 0
        self.alphas = { 0 : {}, 1 : {}, 2 : {} }
        self.positions = {}

    def darts(self): 
        """ 
        Return a list of id representing the darts of the structure 
        """
        return self.alphas[0].keys()

    def alpha(self, degree, dart):
        """ Return the application of the alpha_deg on dart """
        return self.alphas[degree][dart]

    def alpha_composed(self, list_of_alpha_value, dart):
        """ 
        Return the application of a composition of alphas on dart 
        """
        for alpha in list_of_alpha_value:
            dart = self.alpha(alpha, dart)
        return dart

    def is_free(self, degree, dart):
        """ 
        Test if dart is free for alpha_degree (if it is a fixed point) 
        """
        return self.alpha(degree,dart) == dart



    def add_dart(self):
        """ 
        Create a new dart and return its id. 
        Set its alpha_i to itself (fixed points) 
        """
        dart = self.maxid
        self.maxid += 1
        for degree in self.alphas.keys():
            self.alphas[degree][dart] = dart
        return dart













    def is_valid(self):
        """ 
        Test the validity of the structure. 
        Check if there is pending dart for alpha_0 and alpha_1 (fixed point) 
        """
        for dart, alpha_0_of_dart in self.alphas[0].items():
             if dart == alpha_0_of_dart : return False # no fixed point
             if dart != self.alpha(0,alpha_0_of_dart) : return False # alpha_0 is an involution

        for dart, alpha_1_of_dart in self.alphas[1].items():
             if dart == alpha_1_of_dart : return False # no fixed point
             if dart != self.alpha(1,alpha_1_of_dart) : return False # alpha_1 is an involution

        for dart in self.darts(): # alpha_0 alpha_2 is an involution
            if self.alpha_composed([0,2,0,2],dart) != dart: return False

        return True

    def link_darts(self,degree, dart1, dart2): 
        """ 
        Link the two darts with a relation alpha_degree
        """
        assert self.is_free(degree,dart1) and self.is_free(degree,dart2)
        self.alphas[degree][dart1] = dart2
        self.alphas[degree][dart2] = dart1

    def print_alphas(self, withposition = True):
        """ 
        Print for each dart, the value of the different alpha applications.
        """ 
        try:
            from colorama import Style, Fore
        except:
            print "Try to install colorama (pip install colorama) for a better-looking display!"
            for d in self.darts():
                print d," | ",self.alpha(0,d),self.alpha(1,d),self.alpha(2,d)  ,  self.positions.get(d,'') if withposition else ''
        else:
            print "d     α0  α1  α2" +(' positions' if withposition else '')
            for d in self.darts():
                print d," | ",Fore.MAGENTA+str(self.alpha(0,d))," ",Fore.GREEN+str(self.alpha(1,d))," ",Fore.BLUE+str(self.alpha(2,d))," ",Fore.BLACK+str(self.positions.get(d,'') if withposition else '')," ",Style.RESET_ALL 


    def orbit(self, dart, list_of_alpha_value):
        """ 
        Return the orbit of dart using a list of alpha relation.
        Example of use : gmap.orbit(0,[0,1]).
        In Python, you can use the set structure to process only once all darts of the orbit.  
        """
        orbit = []
        marked = set([])
        toprocess = [dart]

        while len(toprocess)>0:
            d = toprocess.pop(0)
            if not d in marked:
                orbit.append(d)
                marked.add(d)
                for degree in list_of_alpha_value:
                    toprocess.append(self.alpha(degree,d))

        return orbit


    def orderedorbit(self, dart, list_of_alpha_value):
        """
        Return the ordered orbit of dart using a list of alpha relations by applying
        repeatingly the alpha relations of the list to dart.
        Example of use. gmap.orderedorbit(0,[0,1]).
        Warning: No fixed point for the given alpha should be contained.
        """
        orbit = []
        current_dart = dart
        current_alpha_index = 0
        n_alpha = len(list_of_alpha_value)
        while (current_dart != dart) or orbit==[]:
            orbit.append(current_dart)
            current_alpha = list_of_alpha_value[current_alpha_index]
            current_dart = self.alpha(current_alpha,current_dart)
            current_alpha_index = (current_alpha_index+1) % n_alpha
        return orbit


    def sew_dart(self, degree, dart1, dart2, merge_attribute = True):
        """
        Sew two elements of degree 'degree' that start at dart1 and dart2.
        Determine first the orbits of dart to sew and heck if they are compatible.
        Sew pairs of corresponding darts, and if they have different embedding 
        positions, merge them. 
        """
        if degree == 1:
            self.link_darts(1, dart1, dart2)
        else:
            alpha_list = [0]
            orbit1 = self.orbit(dart1, alpha_list)
            orbit2 = self.orbit(dart2, alpha_list)
            if len(orbit1) != len(orbit2):
                raise ValueError('Incompatible orbits', orbit1, orbit2)
            for d1,d2 in zip(orbit1, orbit2):
                self.link_darts(degree, d1, d2)
                if merge_attribute:
                    d1e = self.get_embedding_dart(d1, self.positions)
                    d2e = self.get_embedding_dart(d2, self.positions)
                    if d1e in self.positions and d2e in self.positions:
                            pos = (self.positions[d1e] + self.positions[d2e]) / 2.
                            del self.positions[d2e]
                            self.positions[d1e] = pos



    def elements(self, degree):
        """ 
        Return one dart per element of degree. For this, consider all darts as initial set S. 
        Take the first dart d, remove from the set all darts of the orbit starting from d and 
        corresponding to element of degree degree. Take then next element from set S and do the 
        same until S is empty. 
        Return all darts d that were used. """
        
        elements = []
        darts = set(self.darts())

        list_of_alpha_value = range(3)
        list_of_alpha_value.remove(degree)

        while len(darts) > 0:
            dart = darts.pop()
            elementi = self.orbit(dart, list_of_alpha_value)
            darts -= set(elementi)
            elements.append(dart)

        return elements


    def incident_cells(self, dart, degree, incidentdegree):
        """
        Return all the element of degree incidentdegree
        that are incident to the element dart of degree degree.
        (Typically all edges around a point)
        For this iterate over all the dart of the orbit of (dart, degree).
        For each dart d of this orbit, get all the darts coresponding
        to the orbit of the element (d, incidentdegree) and remove them
        from the original set.
        """
        results = []

        alphas = range(3)
        alphas.remove(degree) 

        incidentalphas = range(3)
        incidentalphas.remove(incidentdegree) 

        marked = set()

        for d in self.orbit(dart, alphas):
            if not d in marked:
                results.append(d)
                marked |= set(self.orbit(d, incidentalphas))

        return results
        

    def insert_edge(self, dart):
        """ 
        Insert an edge at the point represented by dart.
        Return a dart corresponding to the dandling edge end.
        """

        dart1 = self.alpha(1, dart)
        newdarts = [self.add_dart() for i in xrange(4)]
        
        self.link_darts(0, newdarts[0], newdarts[1])
        self.link_darts(0, newdarts[3], newdarts[2])
        
        self.link_darts(2, newdarts[0], newdarts[3])
        self.link_darts(2, newdarts[1], newdarts[2])

        self.alphas[1][dart] = newdarts[0]
        self.alphas[1][newdarts[0]] = dart

        self.alphas[1][dart1] = newdarts[3]
        self.alphas[1][newdarts[3]] = dart1

        return newdarts[1]


    def split_face(self, dart1, dart2=None):
        """
        Split face by inserting an edge between dart1 and dart2 
        """

        if dart2 is None:
            dart2 = self.alpha_composed([0,1,0],dart1)

        dedge = self.insert_edge(dart1)

        dart2a1 = self.alpha(1,dart2)
        dedgea2 = self.alpha(2, dedge)

        self.alphas[1][dart2] = dedge
        self.alphas[1][dedge] = dart2

        self.alphas[1][dart2a1] = dedgea2
        self.alphas[1][dedgea2] = dart2a1


    def split_edge(self, dart):
        """ 
        Operator to split an edge. 
        Return a dart corresponding to the new points
        """
        orbit1 = self.orbit(dart,[2])
        orbit2 = self.orbit(self.alpha(0,dart),[2])

        newdart1 = [self.add_dart() for i in orbit1]
        newdart2 = [self.add_dart() for i in orbit2]
        
        for d, nd in zip(orbit1+orbit2, newdart1+newdart2):
            self.alphas[0][d] = nd
            self.alphas[0][nd] = d
        
        for nd1, nd2 in zip(newdart1, newdart2):
            self.link_darts(1, nd1, nd2)
        
        for nd in newdart1+newdart2:
            if self.is_free(2, nd) and not self.is_free(2, self.alpha(0, nd)):
                self.link_darts(2,nd, self.alpha(0,self.alpha(2,self.alpha(0,nd))))

        return newdart1[0]

    

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


    def element_center(self, dart, degree):
        list_of_alpha_value = range(3)
        list_of_alpha_value.remove(degree)

        element_positions = [self.get_position(d) for d in self.orbit(dart,list_of_alpha_value)]
        return np.mean(element_positions,axis=0)


    def dart_display(self, radius=0.1, coef=0.8, add=False):
        import openalea.plantgl.all as pgl

        sphere = pgl.Sphere(radius,slices=16,stacks=16)
        coal = pgl.Material(ambient=(8,10,13),diffuse=3.,specular=(89,89,89),shininess=0.3)
        purple = pgl.Material(ambient=(72,28,72),diffuse=2.,specular=(89,89,89),shininess=0.3)
        green = pgl.Material(ambient=(0,88,9),diffuse=2.,specular=(89,89,89),shininess=0.3)
        blue = pgl.Material(ambient=(9,0,88),diffuse=2.,specular=(89,89,89),shininess=0.3)

        font = pgl.Font(size=10)

        s = pgl.Scene()

        dart_points = {}

        for dart in self.darts():
            dart_point = self.get_position(dart)
            dart_face_center = self.element_center(dart,2)
            dart_edge_center = self.element_center(dart,1)

            dart_face_point = dart_face_center + coef*(dart_point-dart_face_center)
            dart_face_edge_center = dart_face_center + coef*(dart_edge_center-dart_face_center)

            dart_edge_point = dart_face_edge_center + coef*(dart_face_point-dart_face_edge_center)
            dart_middle_edge_point = dart_face_edge_center + 0.33*(dart_edge_point-dart_face_edge_center)

            dart_points[dart] = [dart_edge_point,dart_middle_edge_point]

            s += pgl.Shape(pgl.Translated(dart_points[dart][0],sphere),coal)
            # s += pgl.Shape(pgl.Translated(np.mean(dart_points[dart],axis=0), pgl.Text(str(dart),fontstyle=font)), coal, id=dart)
            s += pgl.Shape(pgl.Polyline(dart_points[dart],width=2),coal)

        for dart in self.darts():
            alpha_0_points = []
            alpha_0_points += [dart_points[dart][1]]
            alpha_0_points += [dart_points[self.alpha(0,dart)][1]]
            s += pgl.Shape(pgl.Polyline(alpha_0_points,width=5),purple)

            alpha_1_points = []
            alpha_1_points += [0.66*dart_points[dart][0] + 0.33*dart_points[dart][1]]
            alpha_1_points += [0.66*dart_points[self.alpha(1,dart)][0] + 0.33*dart_points[self.alpha(1,dart)][1]]
            s += pgl.Shape(pgl.Polyline(alpha_1_points,width=5),green)

            alpha_2_points = []
            alpha_2_points += [0.33*dart_points[dart][0] + 0.66*dart_points[dart][1]]
            alpha_2_points += [0.33*dart_points[self.alpha(2,dart)][0] + 0.66*dart_points[self.alpha(2,dart)][1]]
            s += pgl.Shape(pgl.Polyline(alpha_2_points,width=5),blue)

        if add : 
            pgl.Viewer.add(s)
        else : 
            pgl.Viewer.display(s)

    def eulercharacteristic(self):
        """
        Compute the Euler-Poincare characteristic of the subdivision
        """
        return sum([pow(-1,i)*len(self.elements(i)) for i in xrange(3)])


    def dual(self):
        """
        Compute the dual of the object
        Create a new GMap object with the same darts but reversed alpha relations
        Update the positions of the dual 0-cells as the centers of the 2-cells
        """
        dual_gmap = GMap()
        for alpha_value in [0,1,2]:
            dual_gmap.alphas[alpha_value] = self.alphas[2-alpha_value]

        for face_dart in self.elements(2):
            dual_gmap.set_position(face_dart,self.element_center(face_dart,2))

        return dual_gmap

    def get_embedding_dart(self, dart, propertydict, degree=0):
        """ 
        Check if a dart of the orbit representing the vertex has already been 
        associated with a value in propertydict. If yes, return this dart, else
        return the dart passed as argument.
        """
        alphas = range(3)
        alphas.remove(degree)
        for d in self.orbit(dart, alphas):
            if propertydict.has_key(d):
                return d
        return dart

    def get_position(self, dart):
        """
        Retrieve the coordinates associated to the vertex &lt;alpha_1, alpha_2&gt;(dart) 
        """
        return self.positions.get(self.get_embedding_dart(dart,self.positions,0))


    def set_position(self, dart, position) :
        """
        Associate coordinates with the vertex &lt;alpha_1,alpha_2&gt;(dart)
        """
        self.positions[self.get_embedding_dart(dart,self.positions)] = np.array(position)

def meanpoint(point_list):
    mean = [0,0,0]
    for elmt in point_list:
        mean[0] += elmt[0]
        mean[1] += elmt[1]
        mean[2] += elmt[2]
    mean[0] = mean[0] / len(point_list)
    mean[1] = mean[1] / len(point_list)
    mean[2] = mean[2] / len(point_list)

    return mean

def catmullclark(gmap):
    # Compute a dictionary of face points equal to the face center
    # (use functions elements and element_center)
    facepoints = dict([   ])
    faces = gmap.elements(2)
    for p in faces:
        facepoints[p] = gmap.element_center(p,2)
    # Create a local function to access to the face point from any dart of the face
    # (use function get_embedding_dart with the right degree)
    def get_facepoint(fdart):
        return facepoints[gmap.get_embedding_dart(fdart,facepoints,2)]

    # Define a function to compute the position of an edge points:
    # For this, list the positions of the edge vertices (use incident_cells and get_position)
    # and the face points of the edge faces (use incident_cells and get_facepoint) and
    # return the mean
    def compute_edgepoint(edart):
        corners = dict([])
        for p in incident_cells(edart,0,1):
            corners[p] = gmap.get_position(p)
        for f in incident_cells(edart,0,2):
            corners[f] = get_position(get_facepoint(f))
        newPos = mean(corners)

    # Create a local function to access the edgepoint from any dart of the edge
    def get_edgepoint(edart):
        return edgepoints[gmap.get_embedding_dart(edart,edgepoints,1)]

    edgepoints = dict([  ])
    p2 = gmap.elements(1)
    for brin in p2:
        edgepoints[brin] = compute_edgepoint(brin)
    # Define a function to compute the new position of a vertex:


    def compute_vertexpoint(vdart):
        # For this, compute the mean of edgepoints of incident edges (E),
        iedges = incident_cells(edart,1)
        ep2 = dict([])
        for key in iedges:
            ep2[key] = get_edgepoint(key)
        E = meanpoint(ep2)
        # the mean of facepoints of incident faces (F) and the current
        ifaces = incident_cells(edart,2)
        ep3 = dict([])
        for key in iface:
            ep3[key] = get_facepoint(key)
        F = meanpoint(ep3)
        # position of the vertex (V). (use incident_cells, get_position,
        # get_edgepoint and get_facepoint)
        V = gmap.get_position(edart)
        # Use the valence of the vertex (k, number of incident edges) to 
        # compute the new position:
        k = len (incident_cells(edart,1))
        # V* &lt;- ((k-3)V + 2E + F)/k
        Vstar[0] = ((k-3)*V[0]+2*E[0]+F[0])/k
        Vstar[1] = ((k-3)*V[1]+2*E[1]+F[1])/k
        Vstar[2] = ((k-3)*V[2]+2*E[2]+F[2])/k

        return Vstar

    # Set the new position to the vertex points

    # Create new vertices in the topological structure corresponding
    # to edge points:
    # For this, go through all the edges and split them
    edges = gmap.elements(1)
    InsDarts = []
    for e in edges:
        tmp = gmap.split_edge(e)
        InsDarts.append(tmp)
    # (use split_edge) and set the position to the new vertex.
    # Doing so, fill a list with all the new inserted darts 
    # (vertex orbits of the new darts created by the split)
    edgepoint_darts = []

    # Finally, create the new vertices in the topological structure 
    # corresponding to face points, and connect them to the edge point 
    # vertices inserted previously:
    # Iterate over all the faces
        # Store the postion of the face point corresponding to the face
        # Iterate over all the incident vertices of the face
            # Check if the vertex corresponds to an edge point (inserted 
            # at the previous step)
                # If it is an edge point store it in a list to process
        # Iterate over all the edge points to process
            # Create a new edge from the vertex (use insert_edge)
            # Store the darts at the other end of the edge (vertex orbits
            # of the new darts created by the insertion) in a list to process
        # Iterate over all the end darts to process
            # If the dart d is free by alpha_1:
                # Find its next dart in the new (quad) face :
                # alpha_0(alpha_1(alpha_0(alpha_1(alpha_0(alpha_1(alpha_0(d)))))))
                # Link those two darts by alpha_1
        # Select one of the end darts and set is position to the face point


def add_topo_square(gmap):
    darts = [gmap.add_dart() for i in xrange(8)]
    for i in xrange(4):
        gmap.link_darts(0, darts[2*i], darts[2*i+1])
    for i in xrange(4):
        gmap.link_darts(1, darts[2*i+1], darts[(2*i+2) % 8])
    return darts
 
def topo_square():
    gmap = GMap()
    add_topo_square(gmap)
    return gmap

def square(xsize = 5, ysize  = 5, height = 0):
    gmap = GMap()
    darts = add_topo_square(gmap)
 
    for darti, position in zip([0,2,4,6],[ [xsize, ysize, height], [xsize, -ysize, height] , [-xsize, -ysize, height], [-xsize, ysize, height]]):
        gmap.set_position(darts[darti], position)
    return gmap

def cube(xsize = 5, ysize  = 5 , zsize = 5):
    g = GMap()
    squares = [add_topo_square(g) for i in xrange(6)]
 
    # sew top square to lateral squares
    g.sew_dart(2, squares[0][0], squares[1][1] )
    g.sew_dart(2, squares[0][2], squares[4][1] )
    g.sew_dart(2, squares[0][4], squares[3][1] )
    g.sew_dart(2, squares[0][6], squares[2][1] )
 
    # sew bottom square to lateral squares
    g.sew_dart(2, squares[5][0], squares[1][5] )
    g.sew_dart(2, squares[5][2], squares[2][5] )
    g.sew_dart(2, squares[5][4], squares[3][5] )
    g.sew_dart(2, squares[5][6], squares[4][5] )
 
    # sew lateral squares between each other
    g.sew_dart(2, squares[1][2], squares[2][7] )
    g.sew_dart(2, squares[2][2], squares[3][7] )
    g.sew_dart(2, squares[3][2], squares[4][7] )
    g.sew_dart(2, squares[4][2], squares[1][7] )
 
    for darti, position in zip([0,2,4,6],[ [xsize, ysize, zsize], [xsize, -ysize, zsize] , [-xsize, -ysize, zsize], [-xsize, ysize, zsize]]):
        dart = squares[0][darti]
        g.set_position(dart, position)
     
    for darti, position in zip([0,2,4,6],[ [xsize, -ysize, -zsize], [xsize, ysize, -zsize] , [-xsize, +ysize, -zsize], [-xsize, -ysize, -zsize]]):
        dart = squares[5][darti]
        g.set_position(dart, position)
 
    return g
 

def holeshape(xsize = 5, ysize = 5, zsize = 5, internalratio = 0.5):
    assert 0 < internalratio < 1
 
    g = GMap()
    squares = [add_topo_square(g) for i in xrange(16)]
 
    # sew upper squares between each other
    g.sew_dart(2, squares[0][2], squares[1][1] )
    g.sew_dart(2, squares[1][4], squares[2][3] )
    g.sew_dart(2, squares[2][6], squares[3][5] )
    g.sew_dart(2, squares[3][0], squares[0][7] )
 
    # sew upper squares with external lateral
    g.sew_dart(2, squares[0][0], squares[8][1] )
    g.sew_dart(2, squares[1][2], squares[9][1] )
    g.sew_dart(2, squares[2][4], squares[10][1] )
    g.sew_dart(2, squares[3][6], squares[11][1] )
 
    # # sew upper squares with internal lateral
    g.sew_dart(2, squares[0][5], squares[12][0] )
    g.sew_dart(2, squares[1][7], squares[13][0] )
    g.sew_dart(2, squares[2][1], squares[14][0] )
    g.sew_dart(2, squares[3][3], squares[15][0] )
 
    # sew lower squares between each other
    g.sew_dart(2, squares[4][6], squares[5][1] )
    g.sew_dart(2, squares[5][4], squares[6][7] )
    g.sew_dart(2, squares[6][2], squares[7][5] )
    g.sew_dart(2, squares[7][0], squares[4][3] )
 
    # sew lower squares with external lateral
    g.sew_dart(2, squares[4][0], squares[8][5] )
    g.sew_dart(2, squares[5][6], squares[9][5] )
    g.sew_dart(2, squares[6][4], squares[10][5] )
    g.sew_dart(2, squares[7][2], squares[11][5] )
 
    # sew lower squares with internal lateral
    g.sew_dart(2, squares[4][5], squares[12][4] )
    g.sew_dart(2, squares[5][3], squares[13][4] )
    g.sew_dart(2, squares[6][1], squares[14][4] )
    g.sew_dart(2, squares[7][7], squares[15][4] )
 
    # sew external lateral squares between each other
    g.sew_dart(2, squares[8][7], squares[9][2] )
    g.sew_dart(2, squares[9][7], squares[10][2] )
    g.sew_dart(2, squares[10][7], squares[11][2] )
    g.sew_dart(2, squares[11][7], squares[8][2] )
 
    # sew internal lateral squares between each other
    g.sew_dart(2, squares[12][2], squares[13][7] )
    g.sew_dart(2, squares[13][2], squares[14][7] )
    g.sew_dart(2, squares[14][2], squares[15][7] )
    g.sew_dart(2, squares[15][2], squares[12][7] )
 
    pos = { 
            (0,0) : [xsize,  ysize,  zsize] ,
            (1,2) : [xsize,  -ysize, zsize] ,
            (2,4) : [-xsize, -ysize, zsize] ,
            (3,6) : [-xsize, ysize,  zsize] ,
 
            (0,5) : [xsize*internalratio,  ysize*internalratio,  zsize] ,
            (1,7) : [xsize*internalratio,  -ysize*internalratio, zsize] ,
            (2,1) : [-xsize*internalratio, -ysize*internalratio, zsize] ,
            (3,3) : [-xsize*internalratio, ysize*internalratio,  zsize] ,
 
            (4,1) : [xsize,  ysize,  -zsize] ,
            (5,7) : [xsize,  -ysize, -zsize] ,
            (6,5) : [-xsize, -ysize, -zsize] ,
            (7,3) : [-xsize, ysize,  -zsize] ,
 
            (4,4) : [xsize*internalratio,  ysize*internalratio,  -zsize] ,
            (5,2) : [xsize*internalratio,  -ysize*internalratio, -zsize] ,
            (6,0) : [-xsize*internalratio, -ysize*internalratio, -zsize] ,
            (7,6) : [-xsize*internalratio, ysize*internalratio,  -zsize] ,
          }
 
    for darti, position in pos.items():
        sqid, dartid = darti
        dart = squares[sqid][dartid]
        g.set_position(dart, position)
 
    return g


if __name__ == '__main__':
    mymap = cube()
    catmullclark(mymap)