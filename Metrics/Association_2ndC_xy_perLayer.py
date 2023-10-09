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

# Functions take as input x and y coords of features in a layer 
# Final association function that returns the average factor for each layer 


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

def getcoords(layer):
    coords = []
    features = layer.getFeatures()
    for feature in features:
        geom = feature.geometry()
        x_coord = geom.asPoint().x()
        y_coord= geom.asPoint().y()
        coords.append((x_coord, y_coord))
    return coords

'''
def is_bbox_inside_bbox(inner_bbox, outer_bbox):
    # Get the corner points of the inner_bbox
    bottom_left = QgsPointXY(inner_bbox.xMinimum(), inner_bbox.yMinimum())
    bottom_right = QgsPointXY(inner_bbox.xMaximum(), inner_bbox.yMinimum())
    top_left = QgsPointXY(inner_bbox.xMinimum(), inner_bbox.yMaximum())
    top_right = QgsPointXY(inner_bbox.xMaximum(), inner_bbox.yMaximum())

    # Check if all the corner points are inside the outer_bbox
    return (outer_bbox.contains(bottom_left) and outer_bbox.contains(bottom_right) and
            outer_bbox.contains(top_left) and outer_bbox.contains(top_right))
'''
def is_center_of_bbox_inside_bbox(inner_bbox, outer_bbox):
    # Calculate the center of the inner_bbox
    center_x = (inner_bbox.xMinimum() + inner_bbox.xMaximum()) / 2
    center_y = (inner_bbox.yMinimum() + inner_bbox.yMaximum()) / 2
    center_point = QgsPointXY(center_x, center_y)

    # Check if the center point is inside the outer_bbox
    return outer_bbox.contains(center_point)

def create_buffered_bbox_dict(o_layer):
    buffered_bbox_dict = {}

    for feature_o in o_layer.getFeatures():
        feature_o_id = feature_o['FID']

        # Get the point from the feature's geometry
        geom_o = feature_o.geometry()
        point = geom_o.asPoint()
        
        # Calculate the buffered bounding box for the feature
        bbox = getBoundingBox(point, o_layer.name())

        # Set a default value for buffer_distance
        buffer_distance = 20

        buffered_coordinates = {
            'min_x': bbox.xMinimum() - buffer_distance,
            'max_x': bbox.xMaximum() + buffer_distance,
            'min_y': bbox.yMinimum() - buffer_distance,
            'max_y': bbox.yMaximum() + buffer_distance,
        }

        buffered_bbox = QgsRectangle(
            buffered_coordinates['min_x'],
            buffered_coordinates['min_y'],
            buffered_coordinates['max_x'],
            buffered_coordinates['max_y'],
        )
        
        # Add the buffered bounding box to the dictionary with the feature ID as the key
        buffered_bbox_dict[feature_o_id] = buffered_bbox
    
    return buffered_bbox_dict

def getDistance(x_coords, y_coords, i_layer, o_layer, buffered_bbox_dict):
    distance_dict = {}

    # Create a zip object for iterating over the x and y coordinates together
    xy_coords = zip(x_coords, y_coords)

    for (x_coord, y_coord), feature_j in zip(xy_coords, i_layer.getFeatures()):
        feature_j_id = feature_j['id'] 
        # Replace the geometry of feature_j with a point geometry based on the provided coordinates
        feature_j.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x_coord, y_coord)))

        # Find the corresponding feature in o_layer
        for feature_o in o_layer.getFeatures():
            if feature_o['FID'] == feature_j_id:
                # Calculate the distance between the geometries of feature_j and feature_o
                distance = round(feature_j.geometry().distance(feature_o.geometry()),2)

                # Check if the new icon's bounding box is inside the buffered bounding box
                i_bbox = getBoundingBox(QgsPointXY(x_coord, y_coord), i_layer.name())
                is_inside = is_center_of_bbox_inside_bbox(i_bbox, buffered_bbox_dict[feature_j_id])

                distance_dict[feature_j.id()] =  {"distance": distance, "is_inside": is_inside, "x": x_coord, "y": y_coord}
        
    return distance_dict

def association_factor(distance, is_inside, w = 0.5):
    return round( (w * distance + (1 - w) * int(not is_inside)), 2) #Calculate the objective function for association metric for the 2nd category
# is_inside factor to return 1 when a feature is not inside the bounding box and 0 when it is
# a lower association factor will represent a stronger association (features are closer and inside the bounding box)

def normalize_distance(distance, min_distance, max_distance):
    #if min_distance == max_distance:
        #return 0
    return (distance - min_distance) / (max_distance - min_distance)

def f3_2(x_coords, y_coords, i_layer, o_layer):
    count = 0
    total_association = 0

    buffered_bbox_dict = create_buffered_bbox_dict(o_layer)
    distance_dict = getDistance(x_coords, y_coords, i_layer, o_layer, buffered_bbox_dict)

    min_distance = min((data["distance"] for data in distance_dict.values()), default=None)
    max_distance = max((data["distance"] for data in distance_dict.values()), default=None)
    
    # Skip the normalization and association calculation when there are no distances to process (empty layer)
    if min_distance is not None and max_distance is not None:
        for feature_id, data in distance_dict.items():
            distance = data["distance"]
            is_inside = data["is_inside"]

            normalized_distance = normalize_distance(distance, min_distance, max_distance)
            result = association_factor(normalized_distance, is_inside)
            total_association += result
            count += 1

    # If there are no results, return None
    if count == 0:
        return None

    # Otherwise, return the average association factor
    average_association = round(total_association / count, 2)
    return average_association

'''
#_______________________________________________________________________________________________

group_name = '60center'
g_layer = getLayer(group_name) 

#second_cat_i_layers = g_layer[0:2]
o_layers = getLayer('Icons')

# Get result for only the first 2 layers in the group
for i in range(2):
    i_layer = g_layer[i] #second category layers 
    o_layer = o_layers[i]

    coords = getcoords(i_layer)
    x_coords = [coord[0] for coord in coords]
    y_coords = [coord[1] for coord in coords]

    average_association = f3_2(x_coords, y_coords, i_layer, o_layer)
    print('Average association for layer:',i_layer.name(), average_association)

#_______________________________________________________________________________________________

'''
