#!/usr/bin/env bash

# FMTOMO Result 3D Modelling | Generate a 3D Wavefront/.obj model of a velocity anomaly isosurface as produced by an FMTOMO run (i.e. with grid files in FMTOMO format).
#     Copyright (C) 2025 Yingbo Li
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

bounds=$1
width=$2

gmt begin map pdf
gmt coast -W1p,black -A100 -JM${width}c -R$bounds -Bf --MAP_FRAME_TYPE=plain --MAP_TICK_LENGTH=0
gmt end

inkscape --without-gui --file=map.pdf --export-plain-svg=map.svg
