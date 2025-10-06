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

# Material/color configuration for the Blender models.

# RGB color for the map placed at the top of the 3D model domain. Each channel's value ranges from 0 to 1.
upper_map_rgb = (1,1,1) # White
# RGB color for the map placed at the bottom of the 3D model domain. Each channel's value ranges from 0 to 1.
lower_map_rgb = (1,0,0) # Red
# RGBA color for the Blender model's background. Each channel's value ranges from 0 to 1.
background = (0.03, 0.03, 0.03, 1) # Nearly black

try:
    # Predefined materials for use by scripts/code importing this file via system (not via Blender).
    from processing_functions import Material

    # Materials defined by <name> : dict(rgb=<RGB>,alpha=<alpha [optional]>)
    materials = {"blue":dict(rgb=[0,0,255]),
                 "lightblue_translucent":dict(rgb=[30,30,255],alpha=0.5),
                 "red":dict(rgb=[255,0,0]),
                 "slab_green":dict(rgb=[3,132,79]),
                 }

    materials = {k:Material(k,**spec) for k,spec in materials.items()}
except ModuleNotFoundError:
    # Due to Blender python pathing, the import may not work. Since only the rgb(a)s are needed, it's fine to just proceed.
    pass
