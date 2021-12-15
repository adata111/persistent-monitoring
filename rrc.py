'''
MIT License
Copyright (c) 2019 Fanjin Zeng
This work is licensed under the terms of the MIT license, see <https://opensource.org/licenses/MIT>.  
'''

import numpy as np
from random import random
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from collections import deque
import itertools
import copy

class Line():
    ''' Define line '''
    def __init__(self, p0, p1):
        self.p = np.array(p0)
        self.dirn = np.array(p1) - np.array(p0)
        self.dist = np.linalg.norm(self.dirn)
        self.dirn /= self.dist # normalize

    def path(self, t):
        return self.p + t * self.dirn


def Intersection(line, center, radius):
    ''' Check line-sphere (circle) intersection '''
    a = np.dot(line.dirn, line.dirn)
    b = 2 * np.dot(line.dirn, line.p - center)
    c = np.dot(line.p - center, line.p - center) - radius * radius

    discriminant = b * b - 4 * a * c
    if discriminant < 0:
        return False

    t1 = (-b + np.sqrt(discriminant)) / (2 * a);
    t2 = (-b - np.sqrt(discriminant)) / (2 * a);

    if (t1 < 0 and t2 < 0) or (t1 > line.dist and t2 > line.dist):
        return False

    return True



def distance(x, y):
    return np.linalg.norm(np.array(x) - np.array(y))


def isInObstacle(vex, obstacles, radius):
    for obs in obstacles:
        if distance(obs, vex) < radius:
            return True
    return False


def isThruObstacle(line, obstacles, radius):
    for obs in obstacles:
        if Intersection(line, obs, radius):
            return True
    return False


def nearest(G, vex, obstacles, radius):
    Nvex = None
    Nidx = None
    minDist = float("inf")

    for idx, v in enumerate(G.vertices):
        line = Line(v, vex)
        if isThruObstacle(line, obstacles, radius):
            continue

        dist = distance(v, vex)
        if dist < minDist:
            minDist = dist
            Nidx = idx
            Nvex = v

    return Nvex, Nidx


def newVertex(randvex, nearvex, stepSize):
    dirn = np.array(randvex) - np.array(nearvex)
    length = np.linalg.norm(dirn)
    dirn = (dirn / length) * min (stepSize, length)

    newvex = (nearvex[0]+dirn[0], nearvex[1]+dirn[1])
    return newvex


def window(startpos, endpos):
    ''' Define seach window - 2 times of start to end rectangle'''
    width = endpos[0] - startpos[0]
    height = endpos[1] - startpos[1]
    winx = startpos[0] - (width / 2.)
    winy = startpos[1] - (height / 2.)
    return winx, winy, width, height


def isInWindow(pos, winx, winy, width, height):
    ''' Restrict new vertex insides search window'''
    if winx < pos[0] < winx+width and \
        winy < pos[1] < winy+height:
        return True
    else:
        return False

def isNotVisited(x, path):
    for node in path:
        if node==x:
            return False
    return True

def findCycle(G, v1, v2, x):
    ''' return cycle containing v1, v2 and x '''
    nodes = list(G.neighbors.keys())
    # print("blehble",v1)
    # for neighbor in G.neighbors[v1]:
    #     newCost = dist[curNode] + cost
    #     if newCost < dist[neighbor]:
    #         dist[neighbor] = newCost
    #         prev[neighbor] = curNode

    cycles = []

    # create a queue which stores the paths
    q = deque()
 
    # path to store the current path
    path = []
    path.append(v1)
    q.append(path)
    loop = 0
    while (len(q)>0):
        # print(q)
        path = q.popleft()
        # print(q)
        # print(path)
        last = path[len(path) - 1]
 
        # if last vertex is the desired destination
        # then print the path
        if (last == v2):
            temp = path
            temp.append(x)
            cycles.append(temp)     
 
        # traverse to all the nodes connected to 
        # current vertex and push new path to queue
        # print(G.neighbors)
        # print(last)
        for neighbour, cost in G.neighbors[last]:
            if (isNotVisited(neighbour, path)):
                newpath = copy.deepcopy(path)
                newpath.append(neighbour)
                # print(neighbour)
                # print(newpath)
                q.append(newpath)
                # print(q)

        # print(len(q))
        loop+=1
        # if loop>3:
            # break

    return cycles
            
        
    


class Graph:
    ''' Define graph '''
    def __init__(self, startpos, endpos):
        self.startpos = startpos
        self.endpos = endpos

        self.vertices = [startpos]
        self.edges = []
        self.success = False

        self.vex2idx = {startpos:0}
        self.neighbors = {0:[]}
        self.distances = {0:0.}

        self.sx = endpos[0] - startpos[0]
        self.sy = endpos[1] - startpos[1]

    def add_vex(self, pos):
        try:
            idx = self.vex2idx[pos]
        except:
            idx = len(self.vertices)
            self.vertices.append(pos)
            self.vex2idx[pos] = idx
            self.neighbors[idx] = []
        return idx

    def add_edge(self, idx1, idx2, cost):
        self.edges.append((idx1, idx2))
        self.neighbors[idx1].append((idx2, cost))
        self.neighbors[idx2].append((idx1, cost))


    def randomPosition(self):
        # SampleFree
        rx = random()
        ry = random()

        posx = self.startpos[0] - (self.sx / 2.) + rx * self.sx * 2
        posy = self.startpos[1] - (self.sy / 2.) + ry * self.sy * 2
        return posx, posy


def RRT(startpos, endpos, obstacles, n_iter, radius, stepSize):
    ''' RRT algorithm '''
    G = Graph(startpos, endpos)

    for _ in range(n_iter):
        randvex = G.randomPosition()
        if isInObstacle(randvex, obstacles, radius):
            continue

        nearvex, nearidx = nearest(G, randvex, obstacles, radius)
        if nearvex is None:
            continue

        newvex = newVertex(randvex, nearvex, stepSize)

        newidx = G.add_vex(newvex)
        dist = distance(newvex, nearvex)
        G.add_edge(newidx, nearidx, dist)

        dist = distance(newvex, G.endpos)
        if dist < 2 * radius:
            endidx = G.add_vex(G.endpos)
            G.add_edge(newidx, endidx, dist)
            G.success = True
            #print('success')
            # break
    return G

def near(G, x, radius):
    x_near = []
    for vex in G.vertices:
        if vex == x:
            continue

        dist = distance(vex, x)
        if dist > radius:
            continue

        line = Line(vex, x)
        if isThruObstacle(line, obstacles, radius):
            continue
        x_near.append(vex)
    return x_near

def findbestcycle(cycles):
    ''' returns cycle with minimum cost and vertex[0] of the cycle '''
    best = None
    maxlen = 0
    for cycle in cycles:
        if len(cycle)>maxlen:
            maxlen = len(cycle)
            best = copy.deepcopy(cycle)

    return best, best[0]

def RRT_star(startpos, endpos, obstacles, n_iter, radius, stepSize):
    ''' RRT star algorithm '''
    G = Graph(startpos, endpos)
    fincycles = []

    for it in range(n_iter):
        randvex = G.randomPosition()
        if isInObstacle(randvex, obstacles, radius):
            continue

        nearvex, nearidx = nearest(G, randvex, obstacles, radius)
        if nearvex is None:
            continue

        newvex = newVertex(randvex, nearvex, stepSize)      # steer

        newidx = G.add_vex(newvex)
        # dist = distance(newvex, nearvex)
        # G.add_edge(newidx, nearidx, dist)
        # G.distances[newidx] = G.distances[nearidx] + dist

        # update nearby vertices distance (if shorter)
        x_near = near(G, newvex, radius)
        # print(x_near)
        if(len(x_near) == 0):
            # print("0 ",it)
            dist = distance(newvex, nearvex)
            G.add_edge(newidx, nearidx, dist)
            G.distances[newidx] = G.distances[nearidx] + dist
        elif(len(x_near) == 1):
            # print("1 ",it)
            dist = distance(newvex, x_near[0])
            G.add_edge(G.vex2idx[x_near[0]], newidx, dist )
            # print(G.vertices)
            # print(G.distances)
            # print(G.vex2idx[x_near[0]])
            G.distances[newidx] = G.distances[G.vex2idx[x_near[0]]] + dist
        else:
            # print("else ",it)
            cycles = []
            for chosenvex in itertools.combinations(x_near, 2):
                vex1 = chosenvex[0]
                vex2 = chosenvex[1]
                if vex1 == newvex or vex2 == newvex:
                    continue

                dist = distance(vex1, newvex)
                if dist > radius:
                    continue

                line = Line(vex1, newvex)
                if isThruObstacle(line, obstacles, radius):
                    continue

                dist = distance(vex2, newvex)
                if dist > radius:
                    continue

                line = Line(vex2, newvex)
                if isThruObstacle(line, obstacles, radius):
                    continue


                cycle = findCycle(G, G.vex2idx[vex1], G.vex2idx[vex2], newidx)
                for mycycle in cycle:
                    cycles.append(mycycle)
            # print(cycles)

            bestcycle, x1 = findbestcycle(cycles)
            # print("best",bestcycle)

            fincycles.append(bestcycle)
            dist1 = distance(newvex, G.vertices[x1])
            G.add_edge(x1, newidx, dist1 )
            G.distances[newidx] = G.distances[x1] + dist1

    finalcycle, _ = findbestcycle(fincycles)  
    final_vex_cycle = []  
    for i in finalcycle:
        final_vex_cycle.append(G.vertices[i])
    return G, final_vex_cycle



def dijkstra(G):
    '''
    Dijkstra algorithm for finding shortest path from start position to end.
    '''
    srcIdx = G.vex2idx[G.startpos]
    dstIdx = G.vex2idx[G.endpos]

    # build dijkstra
    nodes = list(G.neighbors.keys())
    dist = {node: float('inf') for node in nodes}
    prev = {node: None for node in nodes}
    dist[srcIdx] = 0

    while nodes:
        curNode = min(nodes, key=lambda node: dist[node])
        nodes.remove(curNode)
        if dist[curNode] == float('inf'):
            break

        for neighbor, cost in G.neighbors[curNode]:
            newCost = dist[curNode] + cost
            if newCost < dist[neighbor]:
                dist[neighbor] = newCost
                prev[neighbor] = curNode

    # retrieve path
    path = deque()
    curNode = dstIdx
    while prev[curNode] is not None:
        path.appendleft(G.vertices[curNode])
        curNode = prev[curNode]
    path.appendleft(G.vertices[curNode])
    return list(path)



def plot(G, obstacles, radius, path=None):
    '''
    Plot RRT, obstacles and shortest path
    '''
    px = [x for x, y in G.vertices]
    py = [y for x, y in G.vertices]
    fig, ax = plt.subplots()

    for obs in obstacles:
        circle = plt.Circle(obs, radius, color='red')
        ax.add_artist(circle)

    ax.scatter(px, py, c='cyan')
    ax.scatter(G.startpos[0], G.startpos[1], c='black')
    # ax.scatter(G.endpos[0], G.endpos[1], c='black')

    lines = [(G.vertices[edge[0]], G.vertices[edge[1]]) for edge in G.edges]
    lc = mc.LineCollection(lines, colors='green', linewidths=2)
    ax.add_collection(lc)

    if path is not None:
        paths = [(path[i], path[i+1]) for i in range(len(path)-1)]
        lc2 = mc.LineCollection(paths, colors='blue', linewidths=3)
        ax.add_collection(lc2)

    ax.autoscale()
    ax.margins(0.1)
    plt.show()


def pathSearch(startpos, endpos, obstacles, n_iter, radius, stepSize):
    G = RRT_star(startpos, endpos, obstacles, n_iter, radius, stepSize)
    if G.success:
        path = dijkstra(G)
        # plot(G, obstacles, radius, path)
        return path


if __name__ == '__main__':
    startpos = (0., 0.)
    endpos = (5., 5.)
    obstacles = [(1., 1.), (2., 2.)]
    n_iter = 800
    radius = 0.5
    stepSize = 0.7

    G, path = RRT_star(startpos, endpos, obstacles, n_iter, radius, stepSize)

    plot(G, obstacles, radius, path)
    # G = RRT(startpos, endpos, obstacles, n_iter, radius, stepSize)


    # if G.success:
    #     path = dijkstra(G)
    #     print(path)
    #     plot(G, obstacles, radius, path)
    # else:
    #     plot(G, obstacles, radius)

