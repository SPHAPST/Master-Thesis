#### ALL THE FUNCTIONS WORK FOR THE GROUP OF LAYERS ###
#### f2t, f2i: For optimization ###


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
# That works for each text layer, so it has to be called twice in the optimization
def f2t(x_coords, y_coords, group_name, t_layer):
    average_factors = []
    group_layers = getLayer(group_name)
    for i_layer in group_layers:
        average_factor = legibility_text(x_coords, y_coords, i_layer, t_layer)
        if average_factor is not None:
            average_factors.append(average_factor)
    return round(sum(average_factors) / len(average_factors),2) if average_factors else None
'''
def f2t(x_coords, y_coords, group_name, t_layers):
    average_factors = []
    group_layers = getLayer(group_name)
    
    for t_layer in t_layers:
        t_layer_factors = []
        for i_layer in group_layers:
            average_factor = legibility_text(x_coords, y_coords, i_layer, t_layer)
            if average_factor is not None:
                t_layer_factors.append(average_factor)
        if t_layer_factors: 
            average_factors.append(sum(t_layer_factors) / len(t_layer_factors))
    
    return round(sum(average_factors) / len(average_factors),2) if average_factors else None

def f2i(x_coords, y_coords, group_name):
    average_factors = {}
    group_layers = getLayer(group_name)

    for i_layer in group_layers:
        legibility_sum = 0
        legibility_count = 0
        for layer in group_layers:
            if layer != i_layer:
                average_legibility = legibility_icons(x_coords, y_coords, i_layer, layer)
                if average_legibility is not None:
                    legibility_sum += average_legibility
                    legibility_count += 1
        if legibility_count > 0:
            average_factors[i_layer.name()] = legibility_sum / legibility_count
        else:
            average_factors[i_layer.name()] = None

    # Filter out None values and calculate the average
    values = [v for v in average_factors.values() if v is not None]
    return round(sum(values) / len(values), 2) if values else None

'''
#_______________________________________________________________________________________________________________________

group_name = '60center'
t_layer = project.mapLayersByName('clip_LMF_Landmark_Building_T')[0]  
t_layer2 = project.mapLayersByName('clip_L_ LMF_Road_Names_T')[0]

coords = getpositions(getLayer(group_name))
x_coords = [coord[0] for coord in coords]
y_coords = [coord[1] for coord in coords]

print('For group:', group_name)  
print(f'Average legibility factor for ROAD text: {f2t(x_coords, y_coords, group_name, t_layer2)}')
print(f'Average legibility factor for LANDMARK text: {f2t(x_coords, y_coords, group_name, t_layer)}')
print(f'Average legibility factor for icons: {f2i(x_coords, y_coords, group_name)}')

#_______________________________________________________________________________________________________________________
'''
'''
#_______________________________________________________________________________________________________________________

group_name = '60center'
t_layers = [project.mapLayersByName('clip_LMF_Landmark_Building_T')[0],
            project.mapLayersByName('clip_L_ LMF_Road_Names_T')[0]]

coords = getpositions(getLayer(group_name))
x_coords = [coord[0] for coord in coords]
y_coords = [coord[1] for coord in coords]

print('For group:', group_name) 
print(f'Average legibility factor for ROAD and LANDMARK text:', {f2t(x_coords, y_coords, group_name, t_layers)})
print(f'Average legibility factor for icons: {f2i(x_coords, y_coords, group_name)}')

#_______________________________________________________________________________________________________________________
'''
