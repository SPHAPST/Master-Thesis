'''
MIT License

Copyright (c) 2023 Sophia Apostolidou

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

#Import packages
from qgis.core import *
from qgis.utils import iface
import math


# Supply path to where is qgis installed
qgs = QgsApplication([], False)
qgs.initQgis()


# Set the project path
path = 'D:\My Documents\Egna_projekt\Textsättning\Testdata\London\London20211006_Clip.qgz'
project = QgsProject.instance()
project.read(path)
#canvas = iface.mapCanvas()

# Default size of a symbol (moving window)
defWidth = 50
defHeight = 50


###############
## FUNCTIONS ##
###############


# Returns a list of all the layers in a group of layers

def getLayer(group_name): #i.e "Icons"
    # Get all layers from the a group (i.e Icons)
    name = group_name
    root = project.layerTreeRoot()
    group = root.findGroup(name)
    #print(group)
    #Initialize a list to store all the layers in the group
    g_layers = []

    # Iterate through all child layers in a group
    for child in group.children():
        # Append the current child layer being iterated through to the list
        g_layers.append(child.layer())
        # Print the current child layer being iterated through
        #print('Child', child)

    return g_layers



# Returns a list of all the features that belong to the layers of a group

def getFeatures(group_name): #i.e "Icons"
    # Get all layers from the a group (i.e Icons)
    name = group_name
    root = project.layerTreeRoot()
    group = root.findGroup(name)
    # Initialize a list to store all the features in the layers
    g_features = []

    # Iterate through all child layers in a group
    for child in group.children():
        # Get the current child layer being iterated through
        layer = child.layer()
        # Append all the features in the layer to the list
        g_features.append([f for f in layer.getFeatures()])

    # Flatten the list of features
    g_features = [i for sub in g_features for i in sub]
    return g_features

# Call the function
#f = getFeatures('Icons')

'''
Retrieve the size of the symbol of the input layer 
and compare it with the default width and height values passed as arguments.
If the size of the symbol is the same as the default width and height, it returns the size of the symbol
otherwise, it returns the default width and height values
as an array of doubles with 2 elements [width, height] as the final size
'''

# Get the size of the feature in a layer

def getSize(theFeature, defWidth, defHeight): #defWidth = default width, theFeature: layer's name
    #initialize an array containing the default width and height values 
    size = [defWidth, defHeight] 
    layer = project.mapLayersByName(theFeature)[0] 
    renderer = layer.renderer() 
    symbol = renderer.symbol() 
    size = symbol.size() 
    #Exctract the width and height values and store them in a list 
    #But first check if size is a float or a QSizeF object and handle it accordingly 
    if isinstance(size, float): 
        size_list = [size, size] 
    else: 
        size_list = [size.width(), size.height()] 
    #Comparison 
    if size_list[0] == defWidth and size_list[1] == defHeight: 
        return size_list 
    else: 
        return [defWidth, defHeight] 
   

# Define function to get a list(tuple) of the bounding box of a feature in a layer 

def getBoundingBox(icon_feature, theFeature): # theFeature = layer's name, icon_feature = the feature of the icon layer
    coords=[]
    # Get the geometry of the feature 
    geometry = icon_feature.geometry() 
    #print('The geometry of the feature:',geometry)
    # Get the x and y coordinates of the point 
    x, y = geometry.asPoint() 
    # Get the size of the marker symbol 
    size = getSize(theFeature,defWidth,defHeight) 
    width, height = size[0], size[1] 
    '''
    # Create a bounding box using the point coordinates and the symbol size 
    bounding_box = QgsRectangle(x - width/2, y - height/2, x + width/2, y + height/2)
    # Get the individual coordinates of the bounding box and store them in a list 
    coords.append([bounding_box.xMinimum(), bounding_box.yMinimum(), bounding_box.xMaximum(), bounding_box.yMaximum()]) 
    '''
        # Calculate the bounding box differently if the feature belongs to the bus stops or cycle stations layer
    if theFeature == "L_ LMF_Bus_Stops_P_new" or theFeature == "L_ LMF_CH_Stations_P_new" or theFeature == "L_ LMF_Bus_Stops_P" or theFeature == "L_ LMF_CH_Stations_P":
        # For specified layers, calculate bounding box based on point touching upper segment of symbol
        bounding_box = QgsRectangle(x - width/2, y - height , x + width/2, y)
    else:
        # For other layers, calculate bounding box based on center point of symbol
        bounding_box = QgsRectangle(x - width/2, y - height/2, x + width/2, y + height/2)
    # Get the individual coordinates of the bounding box and store them in a list 
    coords.append([bounding_box.xMinimum(), bounding_box.yMinimum(), bounding_box.xMaximum(), bounding_box.yMaximum()]) 

    return bounding_box


# Class that holds variables needed in different methods, 
# but are not specific to a particular method, they can be passed around easily

class Parameters:
    def __init__(self):
        self.orgX = 0.0
        self.orgY = 0.0
        self.width = 0.0
        self.height = 0.0

        self.iconSizeGround = 0.0

        self.iconSize = 0.0

        self.n = 0
        self.step = 0.0
        self.stepGround = 0.0
        self.numStep = 0

        self.bbMinX = 0.0
        self.bbMinY = 0.0
        self.bbMaxX = 0.0
        self.bbMaxY = 0.0

        self.numMatrix = [] #Not fixed size

        self.minX = 0
        self.minY = 0

        self.datasets = None


##############
## MATRICES ##
##############


# addToNumMatrixPoint 
def addToNumMatrixPoint(theFeature, weight, p): #theFeature: point feature, p: object of class Parameters
    u = 0
    v = 0  

    # Retrive the geometry  
    thePoint = theFeature.geometry()
    # Get the coordinates of point
    theCoordinate = thePoint.asPoint()

    # Exctract to a list with only x and y values
    #coords = [theCoordinate.x(), theCoordinate.y()]

    # Now, go through all the coordinates
    # Check whether the x and y coordinates of the point lie within the bounding box
    if theCoordinate.x() > p.bbMinX and theCoordinate.x() < p.bbMaxX and theCoordinate.y() > p.bbMinY and theCoordinate.y() < p.bbMaxY:
        # Calculate the indices of the matrix to place the point in
        u = int(math.floor(((theCoordinate.x() - p.bbMinX) / (p.bbMaxX - p.bbMinX)) / p.step))
        v = int(math.floor(((theCoordinate.y() - p.bbMinY) / (p.bbMaxY - p.bbMinY)) / p.step))
        # Increment the value at the calculated indices by the weight
        p.numMatrix[u][v] = p.numMatrix[u][v] + weight

        # Ensure that the indices 'u' and 'v' are within the bounds of the numMatrix
        assert 0 <= u < len(p.numMatrix), f"Index 'u' out of bounds: {u}"
        assert 0 <= v < len(p.numMatrix[0]), f"Index 'v' out of bounds: {v}"
            
        print('The new value of the matrix at the indices',u,v,'is:',p.numMatrix[u][v])

# distance
    # theCoordinates is a list of lists containing the x and y coordinates 
    # Adjust the function to access the x and y coordinates using list indexing
def distance(c1, c2):
    return round(math.sqrt(math.pow(c1[0] - c2[0], 2) + math.pow(c1[1] - c2[1], 2)))


# Calculate the coordinates of a new point that lies at a certain distance (distLimit) 
# from an existing point (c1) on a line between two points (c1 and c2). c2: next point

def computeNewCoordinate(c1, c2, distLimit):
    # Compute x and y value for the new point
    newX = c1[0] + (c2[0]  - c1[0] ) * distLimit / distance(c1,c2)
    newY = c1[1]  + (c2[1] - c1[1]) * distLimit / distance(c1,c2)
    # Creates a new QgsPointXY object with the new x, y values
    theNewCoordinate = QgsPointXY(newX, newY)
    return theNewCoordinate


# addToNumMatrixLineString 
'''
The goal is to add vertices (or points) to a line feature at regular intervals (determined by the distLimit parameter) 
to improve label placement along the line. 
The new points are added to a "coordinateVector" list, which is used later to create a new object with the updated line geometry. 
This allows for more evenly spaced labels along the line feature, as opposed to just using the original vertices of the line.
Labels are placed at more consistent distances along the line, rather than clustering at certain areas or leaving large gaps at others
'''
def addToNumMatrixLineString(theFeature, weight, featureType, distLimit, p): #  theFeature:  QGIS feature object, 
                                                                            # featureType: string representing type of feature being passed in (e.g. "Road")
                                                                            # distLimit: desired distance between the current point and the new point

    u = 0
    v = 0
    currentCoordinate = None
    nextCoordinate = None
    numPoints = 0
    newCoordinateArray = None

    # Get a list of points that construct the line, with their coords 
    # line's coordinates are first extracted as a list of lists containing QgsPointXY objects  
    theLineString = theFeature.geometry().asMultiPolyline() 

    # Extract only the coords of the points. A tuple is created 
    theCoordinates = [] 
    for line in theLineString: 
        for point in line: 
            theCoordinates.append([point.x(), point.y()]) 

    #Create a vector for the new coordinates
    coordinateVector = []
    # Now, go through all the coordinates
    currentCoordinate = theCoordinates[0]
    # Start with adding the first coordinate to the vector
    coordinateVector.append(currentCoordinate)
    # Get the number of points
    numPoints = len(theCoordinates)

    # Iterate through all coords in the line
    for i in range(1, numPoints):
        # Get the next coordinate
        nextCoordinate = theCoordinates[i]
        # While the distance between the current coordinate and the next coordinate 
        # is longer than the desired distance limit
        while distance(currentCoordinate, nextCoordinate) > distLimit:
            # Calculate the coordinates of a new point that is at the desired distance from the current point
            newCoordinate = computeNewCoordinate(currentCoordinate, nextCoordinate, distLimit)
            # Add the new to the list of coords
            coordinateVector.append(newCoordinate)
            # Update the current coords
            currentCoordinate = newCoordinate
        # Append the next coordinate to the list of coords
        coordinateVector.append(nextCoordinate)
        # Update the current coord
        currentCoordinate = nextCoordinate

    # Create a new object that has the new line geometry
    # Currently the line is not allowed to have any holes
    # Start by converting the coordinate vector to a coordinate array
    size = len(coordinateVector)

    # Construct a new object (newCoordinateArray) with the new coordinates
    newCoordinateArray = [coordinateVector[i] for i in range(size)]

    # Its number of coordinates
    numPoints = len(newCoordinateArray)
    print(f"newCoordinateArray: {newCoordinateArray}")

    # Now, go through all the coordinates
    # For each coordinate in newCoordinateArray, check if it is within the bounding box of the numMatrix
    # If so, update the corresponding cell with the weight value
    for i in range(numPoints):
        # instead of [i].x(), [i].y() use [i][0], [i][1]
        if newCoordinateArray[i][0] > p.bbMinX and newCoordinateArray[i][0] < p.bbMaxX and newCoordinateArray[i][1] > p.bbMinY and newCoordinateArray[i][1] < p.bbMaxY:
            u = int(math.floor(((newCoordinateArray[i][0] - p.bbMinX) / (p.bbMaxX - p.bbMinX)) / p.step))
            v = int(math.floor(((newCoordinateArray[i][1] - p.bbMinY) / (p.bbMaxY - p.bbMinY)) / p.step))
            p.numMatrix[u][v] = p.numMatrix[u][v] + weight
            print(f"Updated numMatrix at ({u}, {v}): {p.numMatrix[u][v]}")
        
            # Ensure that the indices 'u' and 'v' are within the bounds of the numMatrix
            assert 0 <= u < len(p.numMatrix), f"Index 'u' out of bounds: {u}"
            assert 0 <= v < len(p.numMatrix[0]), f"Index 'v' out of bounds: {v}"
        print('The new numMatrix for line is: ', p.numMatrix)


# addToNumMatrixPolygon
def addToNumMatrixPolygon(theFeature, weight, p): # theFeature: QgsFeature object, polygon feature
    u = 0
    v = 0

    theGeometry = theFeature.geometry()

    # Check if the feature's geometry is a multipart geometry or not
    if theGeometry.isMultipart(): 
        # Access the first element of the returned multi-polygon list 
        thePolygon = theGeometry.asMultiPolygon()[0] 
    else: 
        thePolygon = theGeometry.asPolygon()

    # Get all the coordinates in the exterior ring
    # Access the exterior ring of the polygon, which is the first element of the thePolygon list 
    theGeometry = thePolygon[0]  
    # A list of tuples 
    theCoordinates = [(point.x(), point.y()) for point in theGeometry]  

    # Convert coords to a list of QgsPointXY objects
    # Go through each tuple in theCoordinates and add the x and y value list 
    #theCoordinates = [x for point in Coordinates for x in [point[0], point[1]]]    
    
    # Its number of coordinates
    numPoints = len(theCoordinates)
    
    # Now, go through all the coordinates and add the coordinates to the numMatrix
    # Only add coordinates that are within the original bounding box
    for i in range(numPoints):
        # instead of [i].x(), [i].y() use [i][0], [i][1]
        if theCoordinates[i][0] > p.bbMinX and theCoordinates[i][0] < p.bbMaxX and theCoordinates[i][1] > p.bbMinY and theCoordinates[i][1] < p.bbMaxY:
            u = int(math.floor(((theCoordinates[i][0] - p.bbMinX) / (p.bbMaxX - p.bbMinX)) / p.step))
            v = int(math.floor(((theCoordinates[i][1] - p.bbMinY) / (p.bbMaxY - p.bbMinY)) / p.step))
            p.numMatrix[u][v] = p.numMatrix[u][v] + weight

            # Ensure that the indices 'u' and 'v' are within the bounds of the numMatrix
            assert 0 <= u < len(p.numMatrix), f"Index 'u' out of bounds: {u}"
            assert 0 <= v < len(p.numMatrix[0]), f"Index 'v' out of bounds: {v}"


#Handles interference with other POI symbols.
#It handles a point like it's a rectangle.
def addSymbolToMatrixPolygon(theFeature, layer_name, weight, p): # theFeature: QgsFeature object
       
    d = getSize(layer_name, defWidth, defHeight) #to access the name of the layer as string
    #print('The size of the symbol is: ' + str(d))

    width = d[0]
    height = d[1]

    #for feature in theFeature.getFeatures():
    # Access feature's geometry and convert it to point
    pnt = theFeature.geometry().asPoint()
    x = pnt.x()
    y = pnt.y()
    #print('The coordinates of the point are: ' + str(x) + ', ' + str(y))

    # Create polygon to represent the icon
    coords = [
        QgsPointXY(x - (width / 2), y - (height / 2)),
        QgsPointXY(x + (width / 2), y - (height / 2)),
        QgsPointXY(x + (width / 2), y + (height / 2)),
        QgsPointXY(x - (width / 2), y + (height / 2)),
        QgsPointXY(x - (width / 2), y - (height / 2))
    ]

    # Create a polygon using the list of coords 
    poly = QgsGeometry.fromPolygonXY([coords]) 
    #print('The polygon is: ' + str(poly))

    # Create a new QgsFeature 
    featP = QgsFeature(theFeature.fields())
    featP.setGeometry(poly)
    #print('The feature is: ' + str(featP))

    # Insert a polygon feature instead of a point one
    addToNumMatrixPolygon(featP, weight, p)
    

# Class to saves the results of the disturbance value of the icon when sliding accross the space
# Access and change (through the attributes) the disturbance values easily

class IconPlacementValue:
    def __init__(self, xDir, yDir, value = 0):
        self.xDir = xDir
        self.yDir = yDir
        self.value = value
        
    def setValue(self, value):
        self.value = value
        
    def getValue(self):
        return self.value

# localiseWindow
# Used to determine the optimal placement of an icon within a search window

'''
HOW IT WORKS?
Try to find the location where the placement of the icon will cause the least disturbance to the background features, 
as measured by the sum of the values in the area covered by the icon. 
The minX and minY variables store the x and y coordinates of this location, respectively. 
Essentially, the code is moving the icon around the background and computing the sum of values under the icon for each location. 
It is then choosing the location where the sum of values is the smallest, 
as this location is assumed to cause the least disturbance to the background.
'''

def localiseWindow(p): #p = Parameters()
    # Compute the midpoint of the search window (number of steps)
    midPoint = int(math.floor(p.numStep / 2))

    # Compute half side of the search window 
    # This distance is the same in both direction 
    # (even though the size of the steps in ground coordinates might differ)
    halfSide = int(math.floor(p.numStep / 2)) # What is the difference with the midPoint?

    # Compute half side of the icon 
    halfIcon = int(math.floor(p.n / 2))

    # Get Total number of placements within the search window
    numPlacements = p.numStep * p.numStep 

    # Create list of IconPlacementValue objects with the same length of numPlacements
    # These objects contain the "disturbance cost" for each icon's placement within the search window
    thePlacements = [None] * numPlacements 
    
    # i stores the number of the placement in a "spiral" pattern
    # sum: help variable for computing the total disturbance value
    i = 0
    sum = 0

    # Compute disturbance value for the window's midpoint
    # Assign the midpoint to the first element
    thePlacements.insert(i, IconPlacementValue(midPoint, midPoint))
    
    # Iterates through a range of values for x and y directions (u and v) within the search window.
    # Sum the disturbance values for all cells within the icon's placement

    for u in range(thePlacements[i].xDir - halfIcon, thePlacements[i].xDir + halfIcon + 1):
        for v in range(thePlacements[i].yDir - halfIcon, thePlacements[i].yDir + halfIcon + 1):
            sum = sum + p.numMatrix[u][v]

    # Assigns the sum value as the disturbance value for the current placement
    thePlacements[i].setValue(sum)
 
    # Define the midpoint disturbance value as the smallest
    minPlacement = i

    # Iterate through the search window in a spiral pattern, starting from the top right corner and spiraling inward
    for k in range(1, halfSide-halfIcon+1):
        i += 1
        # Compute the next point in the pattern (x,y values of previous plus 1)
        thePlacements[i] = IconPlacementValue(thePlacements[i-1].xDir, thePlacements[i-1].yDir+1)
        sum = 0

        # Compute the disturbance value for the current placement
        for u in range(thePlacements[i].xDir-halfIcon, thePlacements[i].xDir+halfIcon+1):
            sum += p.numMatrix[u][thePlacements[i].yDir+halfIcon]
            sum -= p.numMatrix[u][thePlacements[i].yDir-(halfIcon+1)]
        thePlacements[i].setValue(thePlacements[i-1].getValue()+sum)

        # Check if the disturbance value of the current placement is less than the minimum disturbance
        if thePlacements[i].getValue() < thePlacements[minPlacement].getValue():
            minPlacement = i
            
        # Repeat the process for the remaining placements in the search window spiral
        # Do the same in different direction (x+1 and y-1)
        for j in range(1, 2*k-1+1):
            i += 1
            thePlacements[i] = IconPlacementValue(thePlacements[i-1].xDir+1, thePlacements[i-1].yDir)
            sum = 0
            # Loop through y direction
            for v in range(thePlacements[i].yDir-halfIcon, thePlacements[i].yDir+halfIcon+1):
                # add the value at x+halfIcon and y
                sum += p.numMatrix[thePlacements[i].xDir+halfIcon][v]
                # substract the value at x-(halfIcon+1) and y
                sum -= p.numMatrix[thePlacements[i].xDir-(halfIcon+1)][v]
            # Set value of current placement 
            thePlacements[i].setValue(thePlacements[i-1].getValue()+sum)
            # Check if current has the minimum value
            if thePlacements[i].getValue() < thePlacements[minPlacement].getValue():
                minPlacement = i

        # Iterate over the next set of placements.
        # This set is the next column to the left of the previous set
        for j in range(1, 2*k+1):
            i += 1
            # Generate the next placement by moving up one step from the previous placement
            thePlacements[i] = IconPlacementValue(thePlacements[i-1].xDir, thePlacements[i-1].yDir-1)
            sum = 0

            # Iterate over the range of the icon and sum the values in the numMatrix
            for u in range(thePlacements[i].xDir-halfIcon, thePlacements[i].xDir+halfIcon+1):
                # Sum the values in the numMatrix within the range of the icon
                sum += p.numMatrix[u][thePlacements[i].yDir-halfIcon]
                # Subtract the values in the numMatrix that were added in the previous set
                sum -= p.numMatrix[u][thePlacements[i].yDir+(halfIcon+1)]
            # Set the value of the current placement to the sum of the values in the numMatrix
            thePlacements[i].setValue(thePlacements[i-1].getValue()+sum)
            # Check if the current placement is the minimum disturbance value so far
            if thePlacements[i].getValue() < thePlacements[minPlacement].getValue():
                minPlacement = i    
   
        # (x1-1 and y)
        for j in range(1, 2*k+1):
            i += 1
            thePlacements[i] = IconPlacementValue(thePlacements[i-1].xDir-1, thePlacements[i-1].yDir)
            sum = 0
            for v in range(thePlacements[i].yDir-halfIcon, thePlacements[i].yDir+halfIcon+1):
                sum += p.numMatrix[thePlacements[i].xDir-halfIcon][v]
                sum -= p.numMatrix[thePlacements[i].xDir+halfIcon+1][v]
            thePlacements[i].setValue(thePlacements[i-1].getValue()+sum)
            if thePlacements[i].getValue() < thePlacements[minPlacement].getValue():
                minPlacement = i

        for j in range(1, 2*k+1):
            i += 1
            thePlacements[i] = IconPlacementValue(thePlacements[i-1].xDir, thePlacements[i-1].yDir+1)
            sum = 0
            for u in range(thePlacements[i].xDir-halfIcon, thePlacements[i].xDir+halfIcon+1):
                sum += p.numMatrix[u][thePlacements[i].yDir+halfIcon]
                sum -= p.numMatrix[u][thePlacements[i].yDir-halfIcon-1]
            thePlacements[i].setValue(thePlacements[i-1].getValue()+sum)
            if thePlacements[i].getValue() < thePlacements[minPlacement].getValue():
                minPlacement = i 
            
    
    '''
    The minimum value of thePlacements[minPlacement] is being stored in p.minX and p.minY 
    because the for loop iterates over all possible placements of the icon 
    and compares their values, with the lowest value being stored as minPlacement. 
    The x and y coordinates of this minimum placement are then stored in p.minX and p.minY respectively.    
    '''
    
    # Determine the final location of the icon as minX,minY
    p.minX = thePlacements[minPlacement].xDir
    p.minY = thePlacements[minPlacement].yDir

# addTextWindow
# Calculate the center point of the window that is placed on top of the icon location

def addTextWindow(p):
        #Create the text window geometry object (polygon)
        windowCoordinates = [[None] * 2 for i in range(4)]  #an array of 'Coordinate' objects is created (create an array of None values and then assign values to the elements of the array as needed)
                                        #windowCoordinates = [Coordinate() for i in range(5)]
                                        # If I replace None with an empty list [], the elements in the windowCoordinates list will all be references to the same empty list object       
        
        # Searching with a larger icon than displayed
        # Controls the spacing between the icon and the window
        # Determines the direction in which the spiral moves and hence the final location of the icon
        SPACING = 0 
        # Set the x, y coords of the 4 cornerns
        # The window will be placed at the minX and minY location
        windowCoordinates[0][0] = p.bbMinX - p.iconSizeGround/2 + p.minX * p.stepGround + SPACING
        windowCoordinates[0][1] = p.bbMinY - p.iconSizeGround/2 + p.minY * p.stepGround + SPACING
        windowCoordinates[1][0] = p.bbMinX - p.iconSizeGround/2 + p.minX * p.stepGround + (p.iconSizeGround-(2*SPACING))
        windowCoordinates[1][1] = p.bbMinY - p.iconSizeGround/2 + p.minY * p.stepGround + SPACING
        windowCoordinates[2][0] = p.bbMinX - p.iconSizeGround/2 + p.minX * p.stepGround + (p.iconSizeGround-(2*SPACING))
        windowCoordinates[2][1] = p.bbMinY - p.iconSizeGround/2 + p.minY * p.stepGround + (p.iconSizeGround-(2*SPACING))
        windowCoordinates[3][0] = p.bbMinX - p.iconSizeGround/2 + p.minX * p.stepGround + SPACING
        windowCoordinates[3][1] = p.bbMinY - p.iconSizeGround/2 + p.minY * p.stepGround + (p.iconSizeGround-(2*SPACING))

        # Return only center coordinate (center point):
        tbl = [[None] * 2 for i in range(2)]  #[0,0]
        tbl[0] = ( windowCoordinates[0][0] + windowCoordinates[2][0] ) / 2
        tbl[1] = ( windowCoordinates[0][1] + windowCoordinates[2][1] ) / 2
        return tbl 


##################
## THE FUNCTION ##
##################


def placeIcon(layers, icon_feature ,iconlayer_name, searchDistanceGround): #layers = group of layers i.e `Background areas`, icon_feature: the feature in the icon's layer ,iconlayer_name: name of icon's layer (string, i.e `L_ LMF_Building_Entrances_P`)
    # Get the bounding box of the icon feature
    e = getBoundingBox(icon_feature, iconlayer_name)
   
    p = Parameters()
    
    background_layers = getLayer(layers) # A list of all background layers   
    
    icon_layer = project.mapLayersByName(iconlayer_name)
    p.datasets = background_layers + icon_layer # A list with both the background layers and the icon layer            
        
    # Get size of the icon or use defaults
    d = getSize(iconlayer_name, defWidth, defHeight)
    p.width = d[0]
    p.height = d[1]

    p.iconSizeGround = p.width
    p.iconSize = p.iconSizeGround / searchDistanceGround

    # n: number of permitted positions in each direction, must be odd numbers (center point of the grid can be accurately determined)
    # n = (2 * searchDistance - iconSideLength) / resolution
    p.n = 9
    p.step = p.iconSize / p.n
    p.stepGround = p.iconSizeGround / p.n

    # numStep: number of grid cells in each direction
    # should be odd so that the grid has a central cell and equal number of cells on both sides of the central cell. 
    # This helps ensure that the grid is symmetrical 
    # and the central cell can be easily located for computing the weights of cartographic points in each grid cell.
    p.numStep = int(math.floor(1 / p.step))

    # Check if the numStep is odd 
    if p.numStep % 2 == 0:
        p.numStep += 1

    #Initialise the matrix for the grid
    #as a list of lists (2D matrix), with dimensions (p.numStep+1) by (p.numStep+1) 
    '''Each cell is initially set to 0. 
    The +1 in the dimensions is to account for the extra row and column when dividing the bounding box into equal parts based on numStep.
    When the icon features are added to the numMatrix, 
    their values will be added to the corresponding cells, 
    and the localiseWindow function will use this matrix to determine the best placement for the icons'''

    p.numMatrix = [[0 for j in range(p.numStep+1)] for i in range(p.numStep+1)]

    # Iterate through the tuple that contains the bounding box coordinates for all the features in the icon's layer
    # and get the minimum x and y coordinates
    # Botom left corner
    #for i in range(len(e)):
    p.orgX = e.xMinimum() #center: (e[0][0] + e[0][2]) / 2, top-right:e[0][2], botom-left: e[0][0]
    p.orgY = e.yMinimum() #center: (e[0][1] + e[0][3]) / 2, top-right: e[0][3], botom-left: e[0][1]

    # Calculate the bounding box for the grid where the icon will be placed    
    p.bbMinX = p.orgX - p.numStep * p.stepGround / 2
    p.bbMinY = p.orgY - p.numStep * p.stepGround / 2
    p.bbMaxX = p.orgX + p.numStep * p.stepGround / 2
    p.bbMaxY = p.orgY + p.numStep * p.stepGround / 2


    #Loop through all layers in the list of layers and get number of features in each layer
    for featureDataset in p.datasets:
        #count=0   
        
        # Get the name of the layer that the feature belongs to
        layer_name = featureDataset.name()            

        # Access the renderer of the layer (instead of the feature's because 'QgsFeature' has no 'renderer' attribute)
        # and then use the symbolForFeature method on it passing the feature as argument to get feature's symbol
        renderer = featureDataset.renderer()

        # Iterate through all the features in the layer and count them
        for theFeature in featureDataset.getFeatures():
            #print('The feature:',theFeature)
            #count += 1        
            #nrofFeatures = count
            #print(featureDataset.name(), 'has', nrofFeatures, 'features')

            # Create a render context object and pass it as the second argument
            # The symbolForFeature() expects 2 arguments
            context = QgsRenderContext()   
            
            # Do not compare with itself! 
            if theFeature != icon_feature:
                #print('The feature is:', theFeature, 'and belongs to the layer:', layer_name)
       
                # Assign a weight to each background layer
                if layer_name in ['L_ LMF_CH_Pavement_A', 'Pavement']:
                    weight = 2

                elif layer_name in ['Green', 'Railway', 'Water', 'L_ LMF_LONDON_BACKGROUND_A']:
                    weight = 1

                elif layer_name in ['Road','BuiltUp', 'Landmark'] : 
                    weight = 3
                
                elif layer_name in ['NewIcons', 'L_ LMF_Building_Entrances_P', 'L_ LMF_Station_Entrance_P', 'L_ LMF_CP_Entrances_P', 'L_ LMF_Subways_P']:
                    weight = 8 #4

                #else:
                    #print('The layer', layer_name, 'is not in the list of layers to be assigned a weight')


                # Get the geometry of the feature
                geometry = theFeature.geometry()
                #print(geometry.asPoint())
                type = geometry.type() # it returns a number because QgsWkbTypes is an enumeration that assigns an integer value to each of the different geometry types, such as Point: 0; Line: 1; Polygon: 2; 
            
                
                # Check the geometry type of the feature
                if type == 0: #QgsWkbTypes.PointGeometry:
                    #print('The feature is a point')
                    #print('The feature belongs to the layer:', layer_name)
                    # Get the symbol that is associated with a point feature
                    style = renderer.symbolForFeature(theFeature, context)
                    #print('The style is:', style)
                    # Check if it's QgsMarkerSymbol 
                    #if isinstance(style, QgsMarkerSymbol): #Not sure if it corresponds to if style is None
                        #print('The style is a marker symbol')
                    if style is None:
                        addToNumMatrixPoint(theFeature, weight, p)
                        #print ('added to point matrix')
                    else:
                    #This is for handling interference with other POI symbols.
                        addSymbolToMatrixPolygon(theFeature, layer_name, weight, p) 
                        #print('added symbol to polygon matrix') 
                        #print('for symbol:', theFeature, layer_name, weight, p)
                
            
                elif type == 1: #QgsWkbTypes.LineGeometry:
                        addToNumMatrixLineString(theFeature, weight, any, 1000*p.iconSize, p) 
                        #print('added to line matrix')
                        #print('The feature is a line')

                elif type == 2: #QgsWkbTypes.PolygonGeometry:
                        addToNumMatrixPolygon(theFeature, weight, p)
                        #print('added to polygon matrix')
                        #print('The feature is a polygon')
                        #print ('for polygon:', theFeature, layer_name, weight, p)
                    
                else:
                    raise ValueError(str(theFeature) + "has an unknown geometry!") 
                            
            #else:
                #print('The feature is:', theFeature, 'and belongs to the layer:', layer_name, 'so it is the same as the icon feature')
                   
    # Find the optimal location for the text window
    localiseWindow(p)   

    # Add the text (searching) window to the feature data set by creating a new point geometry at that location
    # the specific x and y coordinates of the new point are returned which is used to set the geometry of the icon, 
    # so that the icon will be placed at the same location as the text window on the map.
    pos = addTextWindow(p)

    # Get the disturbance value for the current position
    disturbance_value = p.numMatrix[p.minX][p.minY]

    # Create a new geometry point object   
    newPoint = QgsGeometry.fromPointXY(QgsPointXY(pos[0], pos[1])) 

    #print('The new point is:', newPoint.asPoint())

    # Set the geometry of the icon object to the newly created point object
    icon_feature.setGeometry(newPoint)

    # Return the new x,y coordinate.
    return pos, disturbance_value


############################################################################
## IMPLIMENTATION OF THE ALGORITHM FOR ALL THE FEATURES IN THE ICON LAYER ##
############################################################################
#RUN FOR ALL ICON LAYERS#

lrs = getLayer('Icons')
layer_names = [layer.name() for layer in lrs]

# Store each icon placed after each iteration in the NewIcon layer
b_layers = getLayer('Background areas')
names = [layer.name() for layer in b_layers]
# The empty layer where the new icons will be stored
seq_layer = project.mapLayersByName(names[13])[0]



for icon_layer_name in layer_names:
    icon_layer = project.mapLayersByName(icon_layer_name)
    ifeatures = []
    new_position = []
    isymbol = icon_layer[0].renderer().symbol()
    # Create layer where the rectangle corresponding to the icon is stored
    epsg = icon_layer[0].crs().postgisSrid()
    #print('The epsg is:', epsg)
    uri = "point?crs=epsg:" + str(epsg) + "&field=id:integer&index=yes" 
    new_layer_name = icon_layer_name + "_new"  # create the new layer name
    output_labels = QgsVectorLayer(uri,new_layer_name,'memory') #'memory'
    #output_labels = QgsVectorLayer(uri, new_layer_name, 'ogr')
    
    #print('The output layer is:', output_labels)
    # For the output layer, set the symbol to be the same as the one for the icon
    output_labels.renderer().setSymbol(isymbol)
    
    prov = output_labels.dataProvider()
    
    # Create a list to store the placed features
    placed_features = []

    # Loop through all the features in the icon layer
    for l in icon_layer:
        ifeatures.append( [f for f in l.getFeatures()]) 

    for i in ifeatures:
        # Clear the placed_features list at the beginning of each inner loop
        for feature in i:
            icon_feature = feature
            location, disturbance = placeIcon('Background areas', icon_feature, icon_layer_name, 40) #80
            new_position.append(location)
            prov.addFeature(icon_feature)            
            placed_features.append(icon_feature)
            #print(f"Disturbance value for feature {icon_feature.id()}: {disturbance}")

    seq_layer.dataProvider().addFeatures(placed_features)
    print('The number of features in the sequence layer is:', seq_layer.featureCount())

    # Refresh the map canvas
    #iface.mapCanvas().refresh()
    project.addMapLayer(output_labels)
   
    print('The output layer for', icon_layer_name,'has', output_labels.featureCount(), 'features') 
    #print('Total number of features in all the icon layers is:', sum(output_labels.featureCount()))
    print('The new position of each feature in the icon layer is:', new_position)



'''
#RUN FOR ONE ICON LAYER#
# Get the features from the icon layer to compare them later with all the background features
ifeatures = [] # A list of all the features in the icon layer
icon_layer = project.mapLayersByName('L_ LMF_Airports_P')

new_position = [] # The new position of all the icon features

# Test
# Get the symbol of the icon object to the symbol of the icon layer
isymbol = icon_layer[0].renderer().symbol()
#print ('The symbol is:', isymbol)


# Create layer where the rectangle corresponding to the icon is stored
epsg = icon_layer[0].crs().postgisSrid()
#print('The epsg is:', epsg)
uri = "point?crs=epsg:" + str(epsg) + "&field=id:integer&index=yes"
# The new layer that contains the icon position
output_labels = QgsVectorLayer(uri,'New_placement','memory')
# For the output layer, set the symbol to be the same as the one for the icon
output_labels.renderer().setSymbol(isymbol)
prov = output_labels.dataProvider()

# Loop through all the features in the icon layer
for l in icon_layer:
    ifeatures.append( [f for f in l.getFeatures()]) 
    for i in ifeatures:
        for feature in i:
            icon_feature = feature

            location = placeIcon('Background areas', icon_feature,'L_ LMF_Airports_P', 100)
            new_position.append(location)
            prov.addFeature(icon_feature)
            

# Add the layer with the output icon labels to the map
project.addMapLayer(output_labels)

#print('The output layer is:', output_labels.name())
print('The output layer has', output_labels.featureCount(), 'features')    
#print('It has been added to the map')

# Create a canvas object
#canvas = iface.mapCanvas()

# Set the extent of the canvas to the extent of the icon layer
#canvas.setExtent(output_labels.extent())

# Refresh the canvas to display the icons
#canvas.refresh()
    
print('The new position of each feature in the icon layer is:', new_position)
#print('The number of new features for the icon layer',icon_layer[0].name(),'is:', len(new_position))   
'''
