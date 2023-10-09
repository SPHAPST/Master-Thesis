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
import numpy as np

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

def getDistance3(x_coords, y_coords, i_layer, o_layer):
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

                distance_dict[feature_j.id()] =  {"distance": distance, "x": x_coord, "y": y_coord}
        
    return distance_dict

def getAssociationMetric3(x_coords, y_coords, i_layer, b_layer):
    Association_dict = {}

    # Build the spatial index for the b_layer
    spatialIndexStructure = QgsSpatialIndex()
    for feature in b_layer.getFeatures():
        spatialIndexStructure.insertFeature(feature)

    # Iterate through all coordinates and features in the i_layer together
    for (x, y), feature_j in zip(zip(x_coords, y_coords), i_layer.getFeatures()):

        # Create a point from the coordinates
        point = QgsPointXY(x, y)

        # Get the bounding box of the point
        labelBB = getBoundingBox(point, i_layer.name())
        labelBB_geom = QgsGeometry.fromRect(labelBB)

        # Get the id of the background feature that intersects with the bounding box of the icon feature
        selected_features = spatialIndexStructure.intersects(labelBB)

        # Get the feature object of the background feature that intersects with the bounding box of the icon feature
        for feature_i in b_layer.getFeatures():
            if feature_i.id() in selected_features:

                # Calculate the intersection between labelBB_geom and feature_i's geometry
                intersection_geom = labelBB_geom.intersection(feature_i.geometry())

                # Identify the true intersection geometries
                if intersection_geom.isGeosValid() and intersection_geom.area() > 0:
                    # Calculate the area of the intersection
                    intersection_area = intersection_geom.area()

                    # Calculate the area of labelBB_geom
                    labelBB_area = labelBB.area()

                    # Calculate the percentage of labelBB covered by feature_i
                    coverage_percentage = round((intersection_area / labelBB_area) * 100, 2)

                    # Add the association metric (percentage covered) to a dictionary
                    Association_dict[feature_j.id()] = feature_i.id(), coverage_percentage

    return Association_dict

def association_factor3(overlap, distance, w1=0.1, w2=0.9):
    if overlap is None or distance is None:
        return float('inf')
    return np.round((w1 * overlap + w2 * distance),2)

def normalize_distance(distance, min_distance, max_distance):
    #if min_distance == max_distance:
        #return 0
    return (distance - min_distance) / (max_distance - min_distance)

'''
def f3_3(x_coords, y_coords, i_layer, o_layer, b_layers):
    count = 0
    total_association = 0
    
    distance_dict = getDistance(x_coords, y_coords, i_layer, o_layer)
    
    min_distance = min((data["distance"] for data in distance_dict.values()), default=None)
    max_distance = max((data["distance"] for data in distance_dict.values()), default=None)

    for b_layer in b_layers:
        association_metrics = getAssociationMetric(x_coords, y_coords, i_layer, b_layer)

        # Skip the normalization and association calculation when there are no distances to process (empty layer)
        if min_distance is not None and max_distance is not None and association_metrics:
            for feature_j_id, (feature_i_id, coverage_percentage) in association_metrics.items():
                distance = distance_dict[feature_j_id]["distance"]
                overlap = 100 - coverage_percentage  # Assuming that a smaller overlap value is better

                normalized_distance = normalize_distance(distance, min_distance, max_distance)

                # I calculate overlap as 100 - coverage_percentage. 
                # This implies that overlap could take values in the range [0, 100], not [0, 1]. 
                # To fix this, you should divide overlap by 100

                association = association_factor(overlap/100, normalized_distance)
                
                total_association += association
                count += 1
                
    # If there are no results, return None
    if count == 0:
        return None

    # Otherwise, return the average association factor
    average_association = round(total_association / count, 2)
    return average_association
'''
def f3_3(x_coords, y_coords, i_layer, o_layer, b_layers):
    # Initialize variables for totals
    total_association = 0
    total_count = 0

    for b_layer in b_layers:
        # Calculate the distances
        distance_dict = getDistance3(x_coords, y_coords, i_layer, o_layer)

        min_distance = min((data["distance"] for data in distance_dict.values()), default=None)
        max_distance = max((data["distance"] for data in distance_dict.values()), default=None)

        # Calculate the associations
        Association_dict = getAssociationMetric3(x_coords, y_coords, i_layer, b_layer)

        count = 0
        for feature_id, data in distance_dict.items():
            distance = data["distance"]
            is_inside = feature_id in Association_dict
            if is_inside:
                overlap = Association_dict[feature_id][1] / 100  # Normalize to [0, 1]
            else:
                overlap = 1

            if min_distance is not None and max_distance is not None:
                normalized_distance = normalize_distance(distance, min_distance, max_distance)
            else:
                normalized_distance = 0 if distance == min_distance else 1

            # Calculate the association factor
            result = association_factor3(overlap, normalized_distance)
            total_association += result
            count += 1

        # Add up the counts for each layer
        total_count += count

    # If there are no results, return None
    if total_count == 0:
        return None

    # Otherwise, return the average association factor
    average_association = round(total_association / total_count, 2)
    return average_association

'''
#_______________________________________________________________________________________________

group_name = '60center'

g_layer = getLayer(group_name) 
third_cat_i_layers = g_layer[-8:]


o_layers = getLayer('Icons')
third_cat_o_layers = o_layers[-8:]


b_layers = getLayer('Background areas')
b_layers = b_layers[4:-1] #Exclude the first 4 layers and the last layer from the background layers because in stores the placed icons 


for i in range(-8, 0):  # Adjusted from len(g_layer) to -8
    i_layer = g_layer[i]
    o_layer = o_layers[i] 

    coords = getcoords(i_layer)
    x_coords = [coord[0] for coord in coords]
    y_coords = [coord[1] for coord in coords]

    average_association = f3_3(x_coords, y_coords, i_layer, o_layer, b_layers)
    print(f'Average association for layer {i_layer.name()}: {average_association}')

#_______________________________________________________________________________________________
'''
