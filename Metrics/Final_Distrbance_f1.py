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

#### ALL THE FUNCTIONS WORK FOR THE GROUP OF LAYERS ###
#### f1: For optimization ###

# import necessary QGIS modules
from qgis.core import *

# Supply path to where is qgis installed
qgs = QgsApplication([], False)
qgs.initQgis()

# Set the project path
path = 'D:\My Documents\Egna_projekt\Textsättning\Testdata\London\London20211006_Clip.qgz'
project = QgsProject.instance()
project.read(path)

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

def getSize(theFeature): #defWidth = default width, theFeature: layer's name
    #initialize an array containing the default width and height values 
    size = [] 
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
    
    return size_list

def getBoundingBox(point, theFeature): # point = QgsPointXY(x, y), theFeature = layer's name 
    coords=[]
    #geometry = icon_feature.geometry() 
    #x, y = geometry.asPoint()
    x, y = point.x(), point.y() 
    size = getSize(theFeature) 
    width, height = size[0], size[1] 

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

# Load all the line segments
def load_lines():
    line_layer = project.mapLayersByName('clip_layer_lines')[0]
    line_features = []
    for f in line_layer.getFeatures():
        line_features.append(f)        
    return line_features

# Initialize a list of line features
line_features = load_lines()

# Load the point layer containing the points and extract points
def load_points():
    point_layer = project.mapLayersByName('clip_layer_points')[0]
    point_features = []
    for f in point_layer.getFeatures():
        point_features.append(f)
    return point_features

# Initialize a list of point features
point_features = load_points()

# Return the x and y positions of the icons in a group of layers as a tuple of x,y coordinates
def getpositions(layergroup): # layergroup = getLayer("Icons")
    coords = []
    for layer in layergroup:
        features = layer.getFeatures()
        for feature in features:
            geom = feature.geometry()
            x_coord = geom.asPoint().x()
            y_coord= geom.asPoint().y()
            coords.append((x_coord, y_coord))
    return coords

def disturbance_factor(bp, tl, w1=0.01, w2=0.0037):
    return round(w1 * bp + w2 * tl, 2)

def readability_lines(x_coords, y_coords, i_layer):
    line_layer = project.mapLayersByName('clip_layer_lines')[0]
    spatialIndexStructure_lines = QgsSpatialIndex()
    spatialIndexStructure_lines.addFeatures(line_features)
    readability_lines_factor = {}

    # Iterate over the x and y coordinates
    for x, y in zip(x_coords, y_coords):
        point = QgsPointXY(x, y)
        bbox = getBoundingBox(point, i_layer.name())
        bbox_geom = QgsGeometry.fromRect(bbox)  # Convert bbox to QgsGeometry

        # Get the features that intersect with the bounding box
        ids = spatialIndexStructure_lines.intersects(bbox)
        overlapping_line_features = []
        total_length = 0

        for id in ids:
            line_feature = line_layer.getFeature(id)
            if line_feature.geometry().intersects(bbox):
                overlapping_line_features.append(line_feature)
                intersection = line_feature.geometry().intersection(bbox_geom)
                total_length += intersection.length()
        
        icon_id = None

        # Get the id of the icon feature
        for feature in i_layer.getFeatures():
            if feature.geometry().asPoint().x() == x and feature.geometry().asPoint().y() == y:
                icon_id = feature.id()
                break

        readability_lines_factor[icon_id] = total_length

    return readability_lines_factor

def readability_points(x_coords, y_coords, i_layer):
    point_layer = project.mapLayersByName('clip_layer_points')[0]
    spatialIndexStructure_points = QgsSpatialIndex()
    spatialIndexStructure_points.addFeatures(point_features)
    readability_points_factor = {}

    for x, y in zip(x_coords, y_coords):
        point = QgsPointXY(x, y)
        bbox = getBoundingBox(point, i_layer.name())
        bbox_geom = QgsGeometry.fromRect(bbox)

        # Get the features that intersect with the bounding box
        ids = spatialIndexStructure_points.intersects(bbox)
        overlapping_point_features = []

        for id in ids:
            point_feature = point_layer.getFeature(id)
            if point_feature.geometry().intersects(bbox):
                overlapping_point_features.append(point_feature)
        
        icon_id = None

        # Get the id of the icon feature
        for feature in i_layer.getFeatures():
            if feature.geometry().asPoint().x() == x and feature.geometry().asPoint().y() == y:
                icon_id = feature.id()
                break

        count_points = len(overlapping_point_features)

        readability_points_factor[icon_id] = count_points

    return readability_points_factor

# Return the avg readability factor for each layer
def f1_single(x_coords, y_coords, i_layer):
    line_factors = readability_lines(x_coords, y_coords, i_layer)
    point_factors = readability_points(x_coords, y_coords, i_layer)

    disturbance_factors = {id: disturbance_factor(point_factors.get(id, 0), tl) for id, tl in line_factors.items()}

    # List to store all disturbance factors
    all_disturbance_factors = []

    # Check if disturbance_factors is empty probably because some layers have no features
    if disturbance_factors:
        # Get maximum value among the disturbance factors
        max_value = max(disturbance_factors.values())
        
        if max_value != 0:
            # Normalize disturbance factors to the range 0-1
            normalized_disturbance_factors = {id: round(value / max_value, 2) for id, value in disturbance_factors.items()}

            # Add the normalized disturbance factors to the list
            all_disturbance_factors.extend(normalized_disturbance_factors.values())

    #else:
        #print("Maximum value is zero. Can't normalize.")

    # Calculate and return the average disturbance factor
    return round(sum(all_disturbance_factors) / len(all_disturbance_factors), 2) if all_disturbance_factors else 0

# Return the avg readability factor for the group
def f1(group_name, x_coords, y_coords):
    # Get all layers in the group
    group_layers = getLayer(group_name)

    # Initialize a list to store all disturbance factors
    all_disturbance_factors = []

    # Start index for slicing the coordinate lists
    start_index = 0

    for i_layer in group_layers:
        # Get the number of features in the current layer
        num_features = i_layer.featureCount()

        # Slice the x and y coordinate lists to match the number of features in the current layer
        x_coords_layer = x_coords[start_index : start_index + num_features]
        y_coords_layer = y_coords[start_index : start_index + num_features]

        # Update the start index for the next slice
        start_index += num_features

        # Get the disturbance factors for the current layer using the coordinates
        layer_disturbance_factor = f1_single(x_coords_layer, y_coords_layer, i_layer)

        # Append the disturbance factor to the list of all disturbance factors
        all_disturbance_factors.append(layer_disturbance_factor)

    # Return the average disturbance factor
    return round(sum(all_disturbance_factors) / len(all_disturbance_factors), 2) if all_disturbance_factors else 0

'''
#_______________________________________________________________________________________________________________________

# This will be the X variable in the optimization problem
group_name = '60center'
coords = getpositions(getLayer(group_name))
x_coords = [coord[0] for coord in coords]
y_coords = [coord[1] for coord in coords]

average_f1 = f1(group_name, x_coords, y_coords)
print('Average F1 Disturbance Factor for Group', group_name, ':', average_f1)

#_______________________________________________________________________________________________________________________
'''
