bl_info = {
    "name": "NRA2 Corona Exporter",
    "author": "LD",
    "version": (1, 0),
    "blender": (4, 5, 4),
    "location": "File > Export > NRA2 (.nra2)",
    "description": "Export scene, materials, camera, and geometry to NRA2 format",
    "category": "Export",
}

import bpy
import bmesh
import struct
import os
from mathutils import Vector, Quaternion

def write_geo(obj, basepath):
    mesh = obj.to_mesh()
    mesh.calc_loop_triangles()
    bm = bmesh.new()
    bm.from_mesh(mesh)
    bm.verts.ensure_lookup_table()

    vertices = [v.co[:] for v in bm.verts]
    triangles = [tri.vertices for tri in mesh.loop_triangles]
    geo_filename = os.path.join(basepath, obj.name + ".geo")

    with open(geo_filename, "wb") as f:
        f.write(struct.pack("<I", 0xc01337))  # GEO_MAGIC
        f.write(struct.pack("<I", 2))         # GEO_VERSION
        f.write(struct.pack("<Q", len(triangles)))
        f.write(struct.pack("<Q", 32))
        f.write(struct.pack("<Q", 32 + len(triangles) * 12))
        for tri in triangles:
            f.write(struct.pack("<III", *tri))
        for v in vertices:
            f.write(struct.pack("<fff", *v))

def write_material(mat, index):
    lines = []
    brdf_mat = index
    num_mats = brdf_mat + 1
    

    # Base color
    if mat.use_nodes:
        for node in mat.node_tree.nodes:
            if node.type == 'BSDF_PRINCIPLED':
                c = node.inputs['Base Color'].default_value
                lines.append(f"diffuse # {brdf_mat}")
                lines.append(f"color d {c[0]:.5f} {c[1]:.5f} {c[2]:.5f} {brdf_mat} # {num_mats}")
                '''
                base_color = node.inputs['Base Color'].default_value
                metallic = node.inputs['Metallic'].default_value
                roughness = node.inputs['Roughness'].default_value
                emission = node.inputs['Emission'].default_value
                emission_strength = node.inputs['Emission Strength'].default_value
                specular = node.inputs['Specular'].default_value
                '''
                num_mats += 1
                break

    
    # Emission (process after base color)
    emission_found = False
    if mat.use_nodes:
        for node in mat.node_tree.nodes:
            if node.type == 'EMISSION':
                e = node.inputs["Color"].default_value
                strength = node.inputs["Strength"].default_value
                lines.append(f"color e {e[0]*strength:.5f} {e[1]*strength:.5f} {e[2]*strength:.5f} {brdf_mat} # {num_mats}")
                num_mats += 1
                emission_found = True

    # Combine components
    if emission_found:
        lines.append(f"mult 2 {brdf_mat+1} {brdf_mat+2} {brdf_mat} # {num_mats} {mat.name}")
    else:
        lines.append(f"mult 1 {brdf_mat+1} {brdf_mat} # {num_mats} {mat.name}")

    return lines

def write_camera(filepath, cam_obj):
    cam_data = cam_obj.data
    cam_matrix = cam_obj.matrix_world
    loc = cam_matrix.to_translation()
    rot = cam_matrix.to_quaternion()

    # Duplicate for shutter-close (no motion blur)
    pos_t1 = loc
    rot_t1 = rot

    # Film and lens settings
    film_width = cam_data.sensor_width
    film_height = cam_data.sensor_height
    crop_factor = 1.0  # assume full frame
    aperture = int(cam_data.dof.aperture_fstop) if cam_data.dof.use_dof else 8
    exposure = 100  # dummy value
    focal_length = cam_data.lens
    iso = 100.0

    # Focus
    focus = cam_data.dof.focus_distance if cam_data.dof.use_dof else 10.0
    focus_sensor_offset = 0.0
    speed = 0.0

    cam_filename = os.path.splitext(filepath)[0] + "0.cam"
    with open(cam_filename, "wb") as f:
        f.write(b'CCAM')                          # magic
        f.write(struct.pack("<i", 1))             # version
        f.write(struct.pack("<3f", loc.x * 10, loc.y * 10, loc.z * 10))         # pos
        f.write(struct.pack("<3f", pos_t1.x, pos_t1.y, pos_t1.z))# pos_t1
        f.write(struct.pack("<4f", -rot.w, rot.x, rot.y, rot.z))  # orient
        f.write(struct.pack("<4f", -rot_t1.w, rot_t1.x, rot_t1.y, rot_t1.z))  # orient_t1
        f.write(struct.pack("<f", speed))
        f.write(struct.pack("<f", focus_sensor_offset))
        f.write(struct.pack("<f", focus * 10))
        f.write(struct.pack("<f", film_width))
        f.write(struct.pack("<f", film_height))
        f.write(struct.pack("<f", crop_factor))
        f.write(struct.pack("<i", aperture))
        f.write(struct.pack("<i", exposure))
        f.write(struct.pack("<f", focal_length/100))
        f.write(struct.pack("<f", iso))

    print(f"Exported full camera to {cam_filename}")

def write_nra2(filepath):
    if not filepath.lower().endswith('.nra2'):
        filepath += '.nra2'

    scene = bpy.context.scene
    objects = [obj for obj in scene.objects if obj.visible_get() and obj.type == 'MESH']
    materials = []
    mat_lines = []
    geo_lines = []

    for obj in objects:
        for mat in obj.data.materials:
            if mat.name not in materials:
                materials.append(mat.name)
                mat_index = len(mat_lines)
                mat_lines += write_material(mat, mat_index)

    for obj in objects:
        for mat in obj.data.materials:
            mat_index = materials.index(mat.name)
            geo_lines.append(f"{mat_index + 2} {obj.name}")
            write_geo(obj, os.path.dirname(filepath))

    with open(filepath, "w") as f:
        f.write("black\n")
        f.write(f"{len(mat_lines)}\n")
        f.write("\n".join(mat_lines) + "\n")
        f.write(f"{len(geo_lines)}\n")
        f.write("\n".join(geo_lines) + "\n")

    cam = scene.camera
    if cam:
        write_camera(filepath, cam)
    else:
        print("No camera found in scene!")

    print(f"Exported {filepath} with {len(materials)} materials and {len(objects)} geometry entries.")

class ExportNRA2(bpy.types.Operator):
    bl_idname = "export_scene.nra2_corona"
    bl_label = "Export NRA2 (Corona)"
    filename_ext = ".nra2"
    filepath: bpy.props.StringProperty(subtype="FILE_PATH")

    def execute(self, context):
        write_nra2(self.filepath)
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

def menu_func_export(self, context):
    self.layout.operator(ExportNRA2.bl_idname, text="NRA2 (Corona) Exporter")

def register():
    bpy.utils.register_class(ExportNRA2)
    bpy.types.TOPBAR_MT_file_export.append(menu_func_export)

def unregister():
    bpy.utils.unregister_class(ExportNRA2)
    bpy.types.TOPBAR_MT_file_export.remove(menu_func_export)

if __name__ == "__main__":
    register()
