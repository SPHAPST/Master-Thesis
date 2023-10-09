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

####### LEGIBILITY METRIC #######
#------ Degree of overlap between:  ------#
#------ 1. icons - text labels      ------#
#------ 2. icons - other icons      ------#

####### FOR OPTIMIZATION #######
#------ Define a final function that takes as input x,y coords of the new placed icons 
#------ and returns the legibility metric as average for all the features ------#


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

def getBoundingBox2(geometry, theFeature): # point = QgsPointXY(x, y), theFeature = layer's name 
    coords=[]
    #geometry = icon_feature.geometry() 
    #x, y = geometry.asPoint()
    if isinstance(geometry, QgsPointXY):
        x, y = geometry.x(), geometry.y()
    else:
        x, y = geometry.asPoint().x(), geometry.asPoint().y()    
    
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

def getpositions(layergroup):
    coords = []
    for layer in layergroup:
        features = layer.getFeatures()
        for feature in features:
            geom = feature.geometry()
            x_coord = geom.asPoint().x()
            y_coord= geom.asPoint().y()
            coords.append((x_coord, y_coord))
    return coords

def getcoords(layer):
    coords = []
    features = layer.getFeatures()
    for feature in features:
        geom = feature.geometry()
        x_coord = geom.asPoint().x()
        y_coord= geom.asPoint().y()
        coords.append((x_coord, y_coord))
    return coords

########### 1. Compute the legibility factor between the icons and the text labels (roads, landmarks) ###########

# Return the average legibility factor for each layer

def legibility_text(x, y, i_layer, t_layer):

    spatialIndexStructure_leg = QgsSpatialIndex()
    spatialIndexStructure_leg.addFeatures(t_layer.getFeatures())

    legibility_factors = {}
    total_intersection = 0
    Sum_areas = 0
    i = 0

    icon_features = i_layer.getFeatures()
    for icon_feature in icon_features:
        icon_id = icon_feature.id()
        x = icon_feature.geometry().asPoint().x()
        y = icon_feature.geometry().asPoint().y()

        labelBB = getBoundingBox(QgsPointXY(x, y), i_layer.name())  # Get the bounding box of the icon feature
        labelBB_geometry = QgsGeometry.fromRect(labelBB)

        label_area = labelBB.area()
        Sum_areas += label_area
        intersection_i = 0

        selected_features = spatialIndexStructure_leg.intersects(labelBB)

        for feat in selected_features:
            feature_j = t_layer.getFeature(feat)
            multi_polygon = QgsGeometry.fromWkt(feature_j.geometry().asWkt())
            polygon_list = multi_polygon.asMultiPolygon()
            border = []
            for polygon in polygon_list[0]:
                for point in polygon:
                    border.append(point)

            outlinegeom = QgsGeometry.fromPolygonXY([border])

            if outlinegeom.intersects(labelBB_geometry):
                intersection_area = labelBB_geometry.intersection(outlinegeom).area()
                total_intersection += intersection_area
                intersection_i += intersection_area

    
        legibility_lm_i = round((intersection_i / label_area), 2) #ranges from 0 to 1, 0 means no intersection (best result)
        legibility_factors[icon_id] = legibility_lm_i #dictionary with the legibility factor for each icon feature

    # Calculate the average legibility factor for the layer
    average_factor = round (sum(legibility_factors.values()) / len(legibility_factors), 2) if legibility_factors else None

    #return legibility_factors
    return average_factor

'''
def legibility_icons(x, y, i_layer, i_layer2):
    legibility_factors = {}
    overlaps = 0

    icon_features = i_layer.getFeatures()
    for icon_feature in icon_features:
        icon_id = icon_feature.id()
        x = icon_feature.geometry().asPoint().x()
        y = icon_feature.geometry().asPoint().y()

        iconBB = getBoundingBox2(QgsPointXY(x, y), i_layer.name())  # Get the bounding box of the icon feature

        icon_area = iconBB.area()
        intersection_i = 0

        for feature_j in i_layer2.getFeatures():
            if feature_j != icon_feature:
                feature_j_name = feature_j.id()
                iconBB_j = getBoundingBox2(feature_j.geometry(), i_layer2.name())

                if iconBB_j.intersects(iconBB):
                    overlaps += 1
                    intersection = iconBB_j.intersect(iconBB)
                    intersection_area = intersection.area()
                
                    intersection_i += intersection_area

                    legibility = round((intersection_area / icon_area), 2)

                    #if 0 <= legibility < 1:
 
                    # Store the legibility factor in the dictionary
                    legibility_factors[icon_id] = legibility


    return legibility_factors
'''
def legibility_icons(x, y, i_layer, i_layer2):
    legibility_sum = 0
    legibility_count = 0

    icon_features = i_layer.getFeatures()
    for icon_feature in icon_features:
        x = icon_feature.geometry().asPoint().x()
        y = icon_feature.geometry().asPoint().y()

        iconBB = getBoundingBox2(QgsPointXY(x, y), i_layer.name())  # Get the bounding box of the icon feature
        icon_area = iconBB.area()

        for feature_j in i_layer2.getFeatures():
            if feature_j != icon_feature:
                iconBB_j = getBoundingBox2(feature_j.geometry(), i_layer2.name())

                if iconBB_j.intersects(iconBB):
                    intersection = iconBB_j.intersect(iconBB)
                    intersection_area = intersection.area()
                    legibility = round((intersection_area / icon_area), 2)

                    legibility_sum += legibility
                    legibility_count += 1

    if legibility_count > 0:
        average_legibility = round(legibility_sum / legibility_count, 2)
    else:
        average_legibility = None

    return average_legibility

'''
#_________________________________________________________________________________________


g_layer = getLayer("60center")

t_layer = project.mapLayersByName('clip_L_ LMF_Road_Names_T')[0]
t_layer2 = project.mapLayersByName('clip_LMF_Landmark_Building_T')[0] 

average_factors = []  # List to store the average legibility factors for each layer

average_factors = {} 

print('For ROAD labels')

for i_layer in g_layer:
    coords = getcoords(i_layer)
    x_coords = [coord[0] for coord in coords]
    y_coords = [coord[1] for coord in coords]
    average_factor = legibility_text(x_coords, y_coords, i_layer, t_layer)
    print(f"Average legibility factor for layer {i_layer.name()}: {average_factor}")

print('For LANDMARK labels')

for i_layer in g_layer:
    coords = getcoords(i_layer)
    x_coords = [coord[0] for coord in coords]
    y_coords = [coord[1] for coord in coords]
    average_factor2 = legibility_text(x_coords, y_coords, i_layer, t_layer2)
    print(f"Average legibility factor for layer {i_layer.name()}: {average_factor2}")

for i_layer in g_layer:
    coords = getcoords(i_layer)
    x_coords = [coord[0] for coord in coords]
    y_coords = [coord[1] for coord in coords]

    legibility_sum = 0
    legibility_count = 0
    for layer in g_layer:
        if layer != i_layer:
            average_legibility = legibility_icons(x_coords, y_coords, i_layer, layer)
            if average_legibility is not None:
                legibility_sum += average_legibility
                legibility_count += 1

    if legibility_count > 0:
        average_factors[i_layer.name()] = legibility_sum / legibility_count
    else:
        average_factors[i_layer.name()] = None

print('For ICONS')

# Print the average legibility factors for each layer
for layer_name, average_factor in average_factors.items():
    print(f"{layer_name}: {average_factor}")


#_________________________________________________________________________________________
'''
