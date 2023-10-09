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
#### f3 compines f3_2 and f3_3 to apply different association functions for different layers in the group ####

from Association_2ndC_xy import *
from Association_3rdC_xy import *

# import necessary QGIS modules
from qgis.core import *

# Supply path to where is qgis installed
qgs = QgsApplication([], False)
qgs.initQgis()

# Set the project path
path = 'D:\My Documents\Egna_projekt\Textsättning\Testdata\London\London20211006_Clip.qgz'
project = QgsProject.instance()
project.read(path)

'''
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
'''
def getpositions(layergroup):
    coords = []
    for layer in layergroup:
        features = layer.getFeatures()
        if features: # Only add coordinates if the layer has features
            for feature in features:
                geom = feature.geometry()
                x_coord = geom.asPoint().x()
                y_coord = geom.asPoint().y()
                coords.append((x_coord, y_coord))
        else: # Append special value if layer has no features
            coords.append(None)
    return coords

'''
def f3_2_G(i_layers, o_layers):
    total_association = 0
    total_count = 0

    for i_layer, o_layer in zip(i_layers, o_layers):
        coords = getcoords(i_layer)
        x_coords = [coord[0] for coord in coords]
        y_coords = [coord[1] for coord in coords]

        average_association = f3_2(x_coords, y_coords, i_layer, o_layer)

        if average_association is not None:
            #print(f'Average association for layer {i_layer.name()}: {average_association}') 
            total_association += average_association
            total_count += 1

    if total_count == 0:
        return None

    return round(total_association / total_count, 2)
'''
def f3_2_G(x_coords, y_coords, g_layers, o_layers, n_features):
    total_association = 0
    count = 0
    start_index = 0

    sec_cat_i_layers = g_layers[:2]  # Use the first 2 layers for the second category
    sec_cat_o_layers = o_layers[:2]  # Use the first 2 icon layers

    for i_layer, o_layer, n in zip(sec_cat_i_layers, sec_cat_o_layers, n_features):
        end_index = start_index + n
        x_coords_layer = x_coords[start_index:end_index]
        y_coords_layer = y_coords[start_index:end_index]
        start_index = end_index

        association = f3_2(x_coords_layer, y_coords_layer, i_layer, o_layer)
        if association is not None:
            total_association += association
            count += 1

    if count == 0:
        return None
    else:
        average_association = round(total_association / count, 2)
        return average_association

'''
def f3_3_G(i_layers, o_layers, b_layers):
    total_association = 0
    total_count = 0

    for i_layer, o_layer in zip(i_layers, o_layers):
        coords = getcoords(i_layer)
        x_coords = [coord[0] for coord in coords]
        y_coords = [coord[1] for coord in coords]

        average_association = f3_3(x_coords, y_coords, i_layer, o_layer, b_layers)

        if average_association is not None:
            #print(f'Average association for layer {i_layer.name()}: {average_association}') 
            total_association += average_association
            total_count += 1

    if total_count == 0:
        return None

    return round(total_association / total_count, 2)
'''
'''
def f3_3_G(x_coords, y_coords, g_layers, o_layers, b_layers, n_features):
    # Get the layers for the third category
    third_cat_i_layers = g_layers[-8:]
    third_cat_o_layers = o_layers[-8:]

    total_association = 0
    total_count = 0
    start_index = 0


    for i_layer, o_layer, n in zip(third_cat_i_layers, third_cat_o_layers, n_features):
        end_index = start_index + n
        x_coords_layer = x_coords[start_index:end_index]
        y_coords_layer = y_coords[start_index:end_index]
        start_index = end_index

        distance_dict = getDistance(x_coords_layer, y_coords_layer, i_layer, o_layer)
        min_distance = min((data["distance"] for data in distance_dict.values()), default=None)
        max_distance = max((data["distance"] for data in distance_dict.values()), default=None)

        for b_layer in b_layers:
            Association_dict = getAssociationMetric(x_coords_layer, y_coords_layer, i_layer, b_layer)
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
                result = association_factor(overlap, normalized_distance)
                total_association += result
                count += 1

            total_count += count

    if total_count == 0:
        return None

    average_association = round(total_association / total_count, 2)
    return average_association
'''

def f3(x_coords, y_coords, g_layers, o_layers, b_layers):
    # Get the number of features for each layer
    n_features = [layer.featureCount() for layer in g_layers]

    # Separate the coordinates into 2nd and 3rd categories
    n_features_2nd = sum(n_features[:2])
    n_features_3rd = sum(n_features[2:])

    x_coords_2nd = x_coords[:n_features_2nd]
    y_coords_2nd = y_coords[:n_features_2nd]

    # Calculate association factors for the 2nd category
    average_association_2nd = f3_2_G(x_coords_2nd, y_coords_2nd, g_layers[:2], o_layers[:2], n_features[:2])
    #print('Average association factor for 2nd category:', average_association_2nd)

    # Split the coordinates for each layer in the 3rd category
    x_coords_3rd = []
    y_coords_3rd = []
    start_index = n_features_2nd
    for n in n_features[2:]:
        end_index = start_index + n
        x_coords_3rd.append(x_coords[start_index:end_index])
        y_coords_3rd.append(y_coords[start_index:end_index])
        start_index = end_index

    # Calculate the average association for each layer in the 3rd category
    total_association_3rd = 0
    count = 0
    for i, (x_coords_layer, y_coords_layer) in enumerate(zip(x_coords_3rd, y_coords_3rd)):
        if x_coords_layer is None or y_coords_layer is None:
            #print(f'Average association factor for layer {i + 3}: None (no features)')
            continue
        association = f3_3(x_coords_layer, y_coords_layer, g_layers[i + 2], o_layers[i + 2], b_layers)
        if association is not None:
            total_association_3rd += association
            count += 1
        #print('Average association factor for layer', i + 3, ':', association)

    # Calculate the average association for the 3rd category
    if count > 0:
        average_association_3rd = round(total_association_3rd / count, 2)
    else:
        average_association_3rd = None
    #print('Average association factor for 3rd category:', average_association_3rd)

    # Return the average association for the group
    if average_association_2nd is not None and average_association_3rd is not None:
        return (average_association_2nd + average_association_3rd) / 2
    elif average_association_2nd is not None:
        return average_association_2nd
    elif average_association_3rd is not None:
        return average_association_3rd
    else:
        return None

'''
#_______________________________________________________________________________________________________________________

group_name = '60center'

# Get layers
g_layers = getLayer(group_name)

o_layers = getLayer('Icons')

b_layers = getLayer('Background areas')[4:-1]

# Prepare coordinates for each layer
coords = getpositions(g_layers)
x_coords = [coord[0] for coord in coords]
y_coords = [coord[1] for coord in coords]

# Get number of features for each layer
n_features = [len(getcoords(i_layer)) for i_layer in g_layers]

# Get result for the first 2 layers in the group
#average_association_2nd = f3_2_G(x_coords, y_coords, g_layers, o_layers, n_features)
#print('Average association for the first 2 layers:', average_association_2nd)

# Get result for the group
average_association = f3(x_coords, y_coords, g_layers, o_layers, b_layers)
print('Average association for the group:', average_association)

#_______________________________________________________________________________________________________________________
'''
