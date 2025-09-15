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
import numpy as np
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
# Declare working dir.
basefolder = os.getcwd()
print(basefolder)

# Load parameters.
with open("tmp.txt") as infile:
    scale,max_z,max_ew,max_ns,render,is_recovery = map(float,infile.read().split(" "))
render = bool(int(render))

# Import SVG map.
map_svg = "map.svg"
bpy.ops.import_curve.svg(filepath=map_svg)
# Join all individual SVG beziers into one bezier
objs = list(bpy.data.collections[map_svg].all_objects)
for o in objs:
    o.select_set(True)
bpy.context.view_layer.objects.active = o
bpy.ops.object.join()
combined_map = bpy.context.object
# Remove shadows cast by maps.
combined_map.visible_shadow = False
# Scaling and other geometry options.
combined_map.scale[0] = scale
combined_map.scale[1] = scale
combined_map.location[2] = max_z
combined_map.data.bevel_depth = 1e-2 / scale
combined_map.data.twist_mode = "Z_UP"

# Upper bound map.
l_map_mtl = "map-mtl"
map_mtl = bpy.data.materials.new(l_map_mtl)
map_mtl.use_nodes = False
combined_map.data.materials.append(map_mtl)
map_mtl.diffuse_color = (*materials.upper_map_rgb,1)
map_mtl.roughness = 1

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
map_mtl_lower.diffuse_color = (*materials.lower_map_rgb,1)
map_mtl_lower.roughness = 1

deselect_all()

# Load as many things as possible.
file_bases = ["vgridstrue.in","vgrids.in"]
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
cam.rotation_euler = [np.pi/3,0,-np.pi/4]
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

def render_all(append=""):
    ''' Autorender the 3D model looking from the 4 ordinal directions, and saves these renders with the direction of looking in the filename.

    append | <str> | string to append to the end of the image file saving the render.

    Returns: <None>
    '''
    cam.rotation_euler = [np.radians(60),0,0]
    render_conf = {"SW":(-np.radians(45),[-x0,-y0,max_z+2]),
                   "SE":(np.radians(45),[max_ew+x0,-y0,max_z+2]),
                   "NE":(np.radians(135),[max_ew+x0,max_ns+y0,max_z+2]),
                   "NW":(np.radians(225),[-x0,max_ns+y0,max_z+2]),
                   }
    for orientation,(z_rot,loc) in render_conf.items():
        cam.rotation_euler[2] = z_rot
        cam.location = loc
        bpy.context.scene.render.filepath = os.path.join("render-%s%s.png" % (orientation,append))
        bpy.ops.render.render(write_still=True)
    return

if render:
    # If rendering is requested:
    if not is_recovery:
        # Render the recovered velocity structures.
        render_all()
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
        render_all("-recovered")

deselect_all()
# Move back out of the temp dir.
os.chdir("../")
