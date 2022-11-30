#!/usr/bin/python3
from math import floor

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import numpy as np
from TSPClasses import *
import heapq
import itertools
import copy

class TSPSolver:
    def __init__(self, gui_view):
        self._scenario = None


    def setupWithScenario(self, scenario):
        self._scenario = scenario

    ''' <summary>
        This is the entry point for the default solver
        which just finds a valid random tour.  Note this could be used to find your
        initial BSSF.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of solution, 
        time spent to find solution, number of permutations tried during search, the 
        solution found, and three null values for fields not used for this 
        algorithm</returns> 
    '''

    def defaultRandomTour(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time() - start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
            bssf = TSPSolution(route)
            count += 1
            if bssf.cost < np.inf:
                # Found a valid route
                foundTour = True
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    ''' <summary>
        This is the entry point for the greedy solver, which you must implement for 
        the group project (but it is probably a good idea to just do it for the branch-and
        bound project as a way to get your feet wet).  Note this could be used to find your
        initial BSSF.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution, 
        time spent to find best solution, total number of solutions found, the best
        solution found, and three null values for fields not used for this 
        algorithm</returns> 
    '''

    def greedy(self, time_allowance=60.0):
        pass

    ''' <summary>
        This is the entry point for the branch-and-bound algorithm that you will implement
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution, 
        time spent to find best solution, total number solutions found during search (does
        not include the initial BSSF), the best solution found, and three more ints: 
        max queue size, total number of states created, and number of pruned states.</returns> 
    '''

    #  O(n^2 * b^n) space and time
    def branchAndBound(self, time_allowance=60.0):
        # beginning set up variables - O(1) time and space unless noted
        results = {}
        cities = self._scenario.getCities()  # O(n) space
        ncities = len(cities)
        foundTour = False
        start_time = time.time()
        bssf = self.defaultRandomTour(time_allowance).get('soln')  # best solution so far, TODO time complex
        queue = []
        numBSSFupdates = 0
        MatrixState.pruned = 0
        MatrixState.totalStates = 0


        # set up initial matrix - O(n^2), n being num cities, time and space
        initialMatrix = [[-1 for i in range(ncities)] for j in range(ncities)]  # set up blank matrix for filling in - O(n^2 time and space)
        for city in range(ncities):
            city1 = cities[city]
            for secondCity in range(ncities):
                city2 = cities[secondCity]
                initialMatrix[city][secondCity] = city1.costTo(city2)   # allows me to get the cost between cities
        # put initial matrix into a matrix state - O(n^2) size, O(1) complexity
        initialState = MatrixState(initialMatrix)
        initialState.matrix, initialState.lowerBound = self.matrixReduce(initialMatrix, 0)
        initialState.visitedCities.append(0)
        heapq.heappush(queue, initialState)
        maxQueueSize = 1


        while queue and time.time() - start_time < time_allowance:  # while it's not empty, pop off smallest item O(b^n) space and time - where b is the average num of expanded states  and n is the number of states TODO
            state = heapq.heappop(queue) # O(1)
            branchStates = self.searchStates(state, bssf) #  O(n^3) time, O(n^2) space
            print("running")
            for branchState in branchStates:   # O(bn) space where b is num branches, n is num cities - go through all the possible children of the current lowest node, worst case finds a solution and goes through all cities
                # statement to see if it is a solution - long enough, has all items
                # if it is a solution do this:
                if len(branchState.visitedCities) == len(branchState.matrix[0]):  # O(1) - check. if all the cities are in there, it's a full tour TODO
                    if branchState.lowerBound < bssf.cost:
                        route = []
                        for index in branchState.visitedCities: #O(n) time and space - turn the indicies back to city elements for the solution
                            route.append(cities[index])
                        bssf = TSPSolution(route)
                        bssf.cost = branchState.lowerBound
                        numBSSFupdates += 1
                        for _ in queue:
                            MatrixState.pruneState()
                        foundTour = True
                        break   # because I organized my queue by lower bound, if it gets to this point then its the solution with the lowest possible bound, everything else will be higher
                else:  # O(1) - checking a few items, pushing on one item on queue
                    if branchState.lowerBound < bssf.cost:
                        heapq.heappush(queue, branchState)
                        if len(queue) > maxQueueSize:
                            maxQueueSize += 1

        #assignments - O(1) const time
        end_time = time.time()
        results['cost'] = floor(bssf.cost) if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = numBSSFupdates
        results['soln'] = bssf
        results['max'] = maxQueueSize
        results['total'] = MatrixState.totalStates
        results['pruned'] = MatrixState.pruned
        return results

    ''' <summary>
        This is the entry point for the algorithm you'll write for your group project.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution, 
        time spent to find best solution, total number of solutions found during search, the 
        best solution found.  You may use the other three field however you like.
        algorithm</returns> 
    '''

    def fancy(self, time_allowance=60.0):
        pass


    def matrixReduce(self, matrix, lb):   # O(n^2) time and space
        reducedMatrix = copy.deepcopy(matrix)  # O(n^2) makes a new nxn matrix
        for i in range(len(reducedMatrix[0])):  # O(n^2) time, goes through every index of nxn matrix worst case
            minOfRow = min(reducedMatrix[i])
            if minOfRow != np.inf:
                lb += minOfRow
                for j in range(len(reducedMatrix[0])):
                    reducedMatrix[i][j] -= minOfRow

        for i in range(len(reducedMatrix[0])):  # O(n^2) time, goes through every index of nxn matrix worst case
            column = [j[i] for j in reducedMatrix]
            minOfCol = min(column)
            if minOfCol != np.inf:
                lb += minOfCol
                for j in range(len(reducedMatrix[0])):
                    reducedMatrix[j][i] -= minOfCol
        return reducedMatrix, lb


    #  O(n^3) time, O(n^2) space
    def searchStates(self, matrixState, bssf):
        # if matrixState.matrix[0] == cityGoingTo:
        #     if len(matrixState.matrix[0]) == len(matrixState.visitedCities): # then we have a complete path/tour
        #         return matrixState
        #     else:
        #         return None  # incomplete tour something went wrong
        newStates = []
        for i in range(1, len(matrixState.matrix[0])): # O(n^3) -time, n^2 functions repeated n times , O(n^2) space bc lots of matricies
            lastCity = matrixState.visitedCities[-1]
            item = matrixState.matrix[lastCity][i]
            if item == np.inf:
                continue

            # deep copy other matrix - O(n^2) time and space
            newState = MatrixState(copy.deepcopy(matrixState.matrix))
            newState.lowerBound = matrixState.lowerBound
            newState.visitedCities = copy.deepcopy(matrixState.visitedCities)


            # "visit" new city, then flip all the rows to inf
            newState.lowerBound += matrixState.matrix[lastCity][i]  #O(1) - time and space just a visit
            newState.matrix[i][lastCity] = np.inf  # flip its inverse - O(1)
            for k in range(len(matrixState.matrix[0])):  #O(n) time flip all the rows/columns with it in it - 2 x n rows
                newState.matrix[lastCity][k] = np.inf
                newState.matrix[k][i] = np.inf
                # for j in range(len(matrixState.matrix[0])):
                #     if k == lastCity or j == i:
                #         newState.matrix[k][j] = np.inf


            newState.matrix, newState.lowerBound = self.matrixReduce(newState.matrix, newState.lowerBound)  # O(n^2)
            if newState.lowerBound > bssf.cost: # O(1) comparison, and prune count
                MatrixState.pruneState() # prune the state, aka don't even go down the branch
                continue
            else:     # O(1) - add things to lists
                newState.visitedCities.append(i)
                newStates.append(newState)

        return newStates

''' <summary>
    this is where I store all the matrix states for the branch and bound
'''

class MatrixState:
    pruned = 0
    totalStates = 0
    def __init__(self, matrix=None):
        self.matrix = matrix
        self.lowerBound = math.inf
        self.visitedCities = []  #(1,2,4)
        MatrixState.totalStates += 1

    @staticmethod
    def pruneState():
        MatrixState.pruned += 1

    def __lt__(self, nxt):
        return self.lowerBound < nxt.lowerBound
