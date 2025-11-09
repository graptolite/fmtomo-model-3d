#!/usr/bin/env python3

'''
FMTOMO Result 3D Modelling | Generate a 3D Wavefront/.obj model of a velocity anomaly isosurface as produced by an FMTOMO run (i.e. with grid files in FMTOMO format).
    Copyright (C) 2025 Yingbo Li

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

# Script to be loaded by Blender for handling the collation of 3D objects and render if requested.

import bpy
import os
import sys
import math
# Ensure the current path is recognised by Blender Python as containing loadable modules.
sys.path.insert(0,"./")
import materials

def toggle_selection(action):
    ''' Toggle the selection of objects.

    Returns: <None>
    '''
    bpy.ops.object.select_all(action=action)
    return

def deselect_all():
    ''' Deselect all object.

    Returns: <None>
    '''
    toggle_selection("DESELECT")
    return

def select_all():
    ''' Select all object.

    Returns: <None>
    '''
    toggle_selection("SELECT")
    return

def delete_all():
    ''' Delete all objects.

    Returns: <None>
    '''
    select_all()
    bpy.ops.object.delete()
    return

def new_collection(name):
    ''' Create a new collection.

    Returns: <bpy.types.Collection>
    '''
    # Define new collection
    coll = bpy.data.collections.new(name)
    # Load collection into the scene.
    bpy.context.scene.collection.children.link(coll)
    return coll

def import_svg(scale,map_svg="map.svg",shadow=False):
    ''' Load an SVG into Blender at a specified map scaling.

    scale   | <float> | scale to apply to the svg after importing into Blender.
    map_svg | <str>   | file to import into Blender.
    shadow  | <bool>  | whether the imported SVG is to cast shadows or not.

    Returns: <bpy.types.Object> | the imported svg as a Blender object.
    '''
    # Import SVG map.
    bpy.ops.import_curve.svg(filepath=map_svg)
    # Join all individual SVG beziers into one bezier
    objs = list(bpy.data.collections[map_svg].all_objects)
    for o in objs:
        o.select_set(True)
    bpy.context.view_layer.objects.active = o
    bpy.ops.object.join()
    combined_map = bpy.context.object
    # Set shadow casting option.
    combined_map.visible_shadow = shadow
    # Set scale.
    combined_map.scale[0] = scale
    combined_map.scale[1] = scale
    # Give area to the lines for render visibility.
    combined_map.data.bevel_depth = 1.5*1e-2/scale
    combined_map.data.twist_mode = "Z_UP"
    return combined_map

def color_material(material_key,RGBA,brightness_on=True):
    ''' Give a material an RGBA color with a "glowing" appearance where possible (for the Blender version and for whether nodes usage is turned on) if brightness_on is requested, otherwise just set the color.

    material_key  | <str>            | key of the material to assign the color to.
    RGBA          | <list> [<float>] | RGBA (in 0-1 range) values to color the material with.
    brightness_on | <bool>           | whether to attempt to change additional options to give the material a "glowing" appearance.

    Returns: None
    '''
    # Identify the material using the key.
    mtl = bpy.data.materials[material_key]
    if brightness_on:
        # Attempt to assign the color in the way requested.
        try:
            # Set the emission strength to 1.
            mtl.node_tree.nodes["Principled BSDF"].inputs[28].default_value = 1
            # Set the emission color.
            mtl.node_tree.nodes["Principled BSDF"].inputs[27].default_value = RGBA
            # Specular off.
            mtl.node_tree.nodes["Principled BSDF"].inputs[13].default_value = 0
            # Roughness off.
            mtl.node_tree.nodes["Principled BSDF"].inputs[2].default_value = 0
            # Basecolor matching.
            mtl.node_tree.nodes["Principled BSDF"].inputs[0].default_value = mtl.node_tree.nodes["Principled BSDF"].inputs[27].default_value
            # If setting these options worked, return early.
            return
        except:
            pass
    # If no brightness requested or setting it failed, set the color in a node-less manner.
    mtl.diffuse_color = RGBA
    # Set roughness to 1 to give a matte appearance.
    mtl.roughness = 1
    return

def render_all(camera,append=""):
    ''' Autorender the 3D model looking from the 4 ordinal directions, and save these renders with the direction of looking in the filename.

    camera | <bpy.types.Camera> | camera to use for rendering.
    append | <str>              | string to append to the end of the image file saving the render.

    Returns: <None>
    '''
    cam.rotation_euler = [math.radians(60),0,0]
    render_conf = {"SW":(-math.radians(45),[-x0,-y0,max_z+2]),
                   "SE":(math.radians(45),[max_ew+x0,-y0,max_z+2]),
                   "NE":(math.radians(135),[max_ew+x0,max_ns+y0,max_z+2]),
                   "NW":(math.radians(225),[-x0,max_ns+y0,max_z+2]),
                   }
    for orientation,(z_rot,loc) in render_conf.items():
        cam.rotation_euler[2] = z_rot
        cam.location = loc
        bpy.context.scene.render.filepath = os.path.join("render-%s%s.png" % (orientation,append))
        bpy.ops.render.render(write_still=True)
    return

# Remove all default items upon Blender load.
delete_all()

# Create new collection for handling optical objects (lighting and cameras).
optics_col = new_collection("optics")
# Add sun.
bpy.ops.object.light_add(type='SUN',radius=1,align='WORLD',location=(0,0,0),rotation=(0,0,0))
sun = bpy.context.scene.objects.get("Sun")
# Increase strength of sun.
sun.data.energy = 3
optics_col.objects.link(sun)

deselect_all()

# Move to the temp dir.
os.chdir("tmp")
# Load parameters.
with open("tmp.txt") as infile:
    scale,max_z,max_depth,max_ew,max_ns,render,is_recovery = map(float,infile.read().split(" "))
render = bool(int(render))
# Declare working dir.
basefolder = os.getcwd()
print(basefolder)

# Import SVG map.
combined_map = import_svg(scale,"map.svg")
combined_map.location[2] = max_z

# Upper bound map.
l_map_mtl = "map-mtl"
# New material.
map_mtl = bpy.data.materials.new(l_map_mtl)
map_mtl.use_nodes = True
combined_map.data.materials.append(map_mtl)
# Assign color to material.
color_material(l_map_mtl,(*materials.upper_map_rgb,1))

# Lower bound map
bpy.ops.object.duplicate()
combined_map_lower = bpy.context.object
combined_map_lower.location[2] = 0
l_map_mtl_lower = "map-mtl-lower"
map_mtl_lower = bpy.data.materials.new(l_map_mtl_lower)
map_mtl_lower.use_nodes = False
# Clear material from duplication of pre-existing object.
combined_map_lower.data.materials.clear()
combined_map_lower.data.materials.append(map_mtl_lower)
color_material(l_map_mtl_lower,(*materials.lower_map_rgb,1))

deselect_all()

# Load as many things as possible.
file_bases = ["vgridstrue.in","vgrids.in","interfaces.in"]
files = []
all_files = os.listdir()
for base in file_bases:
    files.extend([f for f in all_files if base in f and f.endswith(".obj")])
print(files)
for f in files:
    path = os.path.join(basefolder,f)
    if os.path.exists(path):
        bpy.ops.wm.obj_import(filepath=path)
        active_collection = new_collection(f.split(".")[0])
        for o in bpy.context.selected_objects:
            active_collection.objects.link(o)
            if "vgrids.in" in f:
                # Apply smoothing to isosurfaces if of the velocity grid.
                smooth_mod = o.modifiers.new(name="Smooth",type='SMOOTH')
                smooth_mod.factor = 0.5
                smooth_mod.iterations = 2
                bpy.context.view_layer.objects.active = o
                bpy.ops.object.modifier_apply(modifier="Smooth")
                bpy.ops.object.shade_smooth()
                deselect_all()
                if not is_recovery:
                    # No shadows
                    o.visible_shadow = False

# Remove the base Collection.
bpy.data.collections.remove(bpy.data.collections["Collection"])

deselect_all()

# Camera stuff
bpy.ops.object.camera_add(enter_editmode=False,align='VIEW',location=(0,0,0),rotation=(0,0,0),scale=(1,1,1))
cam = bpy.context.scene.objects.get("Camera")
# Different coord convention.
x0,y0 = max_ns/2,max_ew/2
# Position the camera.
cam.rotation_euler = [math.pi/3,0,-math.pi/4]
cam.location = [-x0,-y0,max_z+2]
optics_col.objects.link(cam)
bpy.context.scene.render.resolution_x = 1024
bpy.context.scene.render.resolution_y = 1024

# Orthographic camera.
bpy.context.object.data.type = 'ORTHO'
bpy.context.object.data.ortho_scale = 20
bpy.context.object.data.shift_y = 0.15

bpy.context.scene.camera = cam
background = materials.background
if materials.background == None:
    bpy.context.scene.render.film_transparent = True
else:
    # Set black background
    bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[0].default_value = background

if render:
    # If rendering is requested:
    if not is_recovery:
        # Render the recovered velocity structures.
        render_all(cam)
    else:
        colls = [c.name for c in bpy.data.collections]
        # Get list of collections containing recovered velocity structure.
        coll_recv = [c for c in colls if "vgrid" in c and "true" not in c]
        # Get list of collections containing true velocity structure.
        coll_true = [c for c in colls if "true" in c]
        # Hide all recovered velocity structures.
        for c in coll_recv:
            bpy.data.collections[c].hide_render = True
        # Render the true velocity structures.
        render_all("-true")
        # Show all recovered velocity structures and hide all true velocity structures.
        for c in coll_recv:
            bpy.data.collections[c].hide_render = False
        for c in coll_true:
            bpy.data.collections[c].hide_render = True
        # Render the recovered velocity structures.
        render_all(cam,"-recovered")

deselect_all()
# Move back out of the temp dir.
os.chdir("../")
