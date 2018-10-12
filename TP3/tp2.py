# coding: utf8

try:
    import gmap_solution as gmap 
except:
    import gmap
reload(gmap)


def question1():
    mgmap = gmap.topo_square()
    mgmap.print_alphas(False)
    return mgmap

def question2():
    mgmap = question1()
    expected = [4,4,1]
    for i in range(3):
        ielements = mgmap.elements(i)
        print "Element de degree ",i,":", ielements
        assert len(ielements) == expected[i]

def question3():
    mgmap = gmap.square()
    mgmap.print_alphas()
    return mgmap

def question4a():   
    mgmap = gmap.cube()
    mgmap.print_alphas()
    return mgmap

def question4b(): 
    mgmap = gmap.holeshape()
    mgmap.print_alphas()
    return mgmap

def question5a():
    mgmap = question4a()
    mgmap.display()

def question5b():
    mgmap = question4b()
    mgmap.display()

def question5c():
    mgmap = question4a()
    mgmap.dart_display()

def question5d():
    mgmap = question4b()
    mgmap.dart_display()

def question6():
    mgmap = gmap.cube()
    print 'Euler characteristic :', mgmap.eulercharacteristic()

def question7():
    mgmap = gmap.cube()
    mgmap = mgmap.dual()
    mgmap.print_alphas()
    mgmap.display()


if __name__ == '__main__':
    import sys
    # On recupere l'id de la question. Par defaut c'est 1
    question = 1
    if len(sys.argv) > 1:
        # sinon on prend le chiffre passe en argument
        question = sys.argv[1]
    # dans le namespace courant on recupere la fonction questioni
    test = globals()['question'+str(question)]
    # on appelle la fonction
    test()