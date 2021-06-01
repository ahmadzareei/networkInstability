module VoronoiPy
export Random3DNetwork
using PyCall

function __init__()
    py"""
    from scipy.spatial import Voronoi
    import networkx as nx
    import numpy as np

    def check_outside(point,sizes):
        for i in range(3):
            if point[i]<0 or point[i]>(sizes[i]-1):
                return True
        return False

    def random_3d_graph_voro(Nx,Ny,Nz,*,sigma = 0.1,randomSeed=0):
        np.random.seed(randomSeed)
        deviation = np.zeros((Nx*Ny*Nz,3))
        deviation[Nx*Ny:-Nx*Ny] = np.random.normal(0, sigma, (Nx*Ny*(Nz-2),3))

        points = []
        for z in range(Nz):
            for y in range(Ny):
                for x in range(Nx):
                    points.append([x,y,z])

        points += deviation
        vor = Voronoi(points)

        edges = set()
        for ridge in vor.regions:
            for i in range(len(ridge)):
                if ridge[i-1] != -1 and ridge[i] != -1:
                    edges.add((ridge[i-1],ridge[i]))

        posArray = vor.vertices
        topNodes = []
        botNodes = []
        removed_nodes =[]
        G = nx.Graph(list(edges))
        node_mapping = dict()
        j = 0
        for i in range(len(G.nodes())):
            node = i
            if check_outside(posArray[node],[Nx,Ny,Nz]):
                removed_nodes.append(node)
            else:
                node_mapping[node] = j
                j += 1
        G.remove_nodes_from(removed_nodes)
        G = nx.relabel_nodes(G,node_mapping)
        posArray=np.delete(posArray,np.array(removed_nodes),axis=0)

        removed_nodes =[]
        node_mapping = dict()
        j = 0
        for i in range(len(G.nodes())):
            node = i
            if len(G[node])==0:
                removed_nodes.append(node)
            else:
                node_mapping[node] = j
                j += 1
        if len(removed_nodes)>0:
            G.remove_nodes_from(removed_nodes)
            G = nx.relabel_nodes(G,node_mapping)
            posArray=np.delete(posArray,np.array(removed_nodes),axis=0)

        for node in G.nodes():
            if abs(posArray[node][1]-Ny/2+0.5)< Ny/4:
                if abs(posArray[node][2]-Nz/2+0.5)< Nz/4:
                    if posArray[node][0]<1:
                        botNodes.append(node+1)
                    elif posArray[node][0]>(Nx-2):
                        topNodes.append(node+1)
        sizes = [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
        if len(sizes)>1:
            return random_3d_graph_voro(Nx,Ny,Nz,sigma = sigma,randomSeed=randomSeed*2)
        else:
            return nx.adjacency_matrix(G,nodelist=np.arange(len(G))).todense(),posArray,np.array(topNodes),np.array(botNodes)
    """
end

Random3DNetwork(Nx,Ny,Nz,sigma,rnd_seed) = py"random_3d_graph_voro"(Nx,Ny,Nz,sigma=sigma,randomSeed=rnd_seed)

end
