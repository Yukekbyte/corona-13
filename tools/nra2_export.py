bl_info = {
    "name": "NRA2 Raytracer Exporter",
    "author": "LD",
    "version": (1, 3),
    "blender": (4, 1, 0),
    "location": "File > Export > NRA2 (.nra2)",
    "description": "Native binary .geo and .cam export for NRA2",
    "category": "Export",
}

import bpy
import struct
import os
from bpy_extras.io_utils import ExportHelper

# --- Bitwise Encoders ---

def encode_uv(u, v):
    # geo_encode_uv uses half-floats packed into 32-bit uint
    return struct.unpack('<I', struct.pack('<ee', u, v))[0]

def to_uint(f):
    return struct.unpack('<I', struct.pack('<f', f))[0]

def encode_normal(x, y, z):
    # Matches the octahedral encoding in geo.h
    invL1Norm = 1.0 / (abs(x) + abs(y) + abs(z) + 1e-8)
    enc0, enc1 = x * invL1Norm, y * invL1Norm
    if z < 0.0:
        enc0_old, enc1_old = enc0, enc1
        enc0 = (1.0 - abs(enc1_old)) * (-1.0 if enc0_old < 0.0 else 1.0)
        enc1 = (1.0 - abs(enc0_old)) * (-1.0 if enc1_old < 0.0 else 1.0)
    
    # Quantization to match fixed-point logic in geo.h
    enci0 = int(((enc0 + 1.0) * 0.5) * 65535) & 0xFFFF
    enci1 = int(((enc1 + 1.0) * 0.5) * 65535) & 0xFFFF
    return (enci1 << 16) | enci0

# --- Binary .geo Writer ---

def write_geo(obj, basepath):
    depsgraph = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(depsgraph)
    mesh = obj_eval.to_mesh()
    mesh.calc_loop_triangles()
    
    # World matrix to bake scale/rotation into the binary data
    matrix = obj.matrix_world
    uv_layer = mesh.uv_layers.active.data if mesh.uv_layers else None

    prims_data = []      # list of primid_t structs
    vtxidx_data = []     # list of prims_vtxidx_t structs
    unique_verts = []    # list of (co, n_enc)
    vert_map = {}        # {(co, n_enc): index}

    vtxidx_global_count = 0

    for tri in mesh.loop_triangles:
        # primid_t: shapeid(4), vi(4), mb(1), vcnt(1), pad(2)
        prims_data.append(struct.pack("<IIBBxx", 0, vtxidx_global_count, 0, 3))
        
        for loop_idx in tri.loops:
            loop = mesh.loops[loop_idx]
            co = matrix @ mesh.vertices[loop.vertex_index].co
            no = (matrix.to_3x3() @ loop.normal).normalized()
            
            n_enc = encode_normal(no.x, no.y, no.z)
            u, v = uv_layer[loop_idx].uv if uv_layer else (0.0, 0.0)
            uv_enc = encode_uv(u, v)

            v_key = (tuple(co), n_enc)
            if v_key not in vert_map:
                vert_map[v_key] = len(unique_verts)
                unique_verts.append(v_key)
            
            # prims_vtxidx_t: v(4), uv(4)
            vtxidx_data.append(struct.pack("<II", vert_map[v_key], uv_enc))
            vtxidx_global_count += 1

    # File Structure Offsets
    num_prims = len(prims_data)
    header_size = 32
    prims_size = num_prims * 12
    vtxidx_size = len(vtxidx_data) * 8
    
    vtxidx_offset = header_size + prims_size
    # Vertices must be 16-byte aligned as per geo.h/obj2geo logic
    vtx_offset = (vtxidx_offset + vtxidx_size + 15) & ~15
    padding_needed = vtx_offset - (vtxidx_offset + vtxidx_size)

    geo_path = os.path.join(basepath, f"{obj.name}.geo")
    with open(geo_path, "wb") as f:
        # prims_header_t: magic, version, num_prims, vtxidx_offset, vertex_offset
        f.write(struct.pack("<IIQQQ", 0xc01337, 2, num_prims, vtxidx_offset, vtx_offset))
        f.write(b"".join(prims_data))
        f.write(b"".join(vtxidx_data))
        f.write(b"\x00" * padding_needed)
        for (co, n_enc) in unique_verts:
            # prims_vtx_t: v[3](12), n(4)
            f.write(struct.pack("<fffI", co[0], co[1], co[2], n_enc))

    obj_eval.to_mesh_clear()

# --- Camera Writer ---

def write_camera(filepath, cam_obj):
    scene = bpy.context.scene
    loc = cam_obj.matrix_world.translation
    # Convert Blender orientation (Z-up) to Raytracer expectation if needed
    # Using raw world matrix translation and quat from legacy script
    rot = cam_obj.matrix_world.to_quaternion()
    
    cam_path = os.path.splitext(filepath)[0] + "0.cam"
    with open(cam_path, "wb") as f:
        f.write(struct.pack("<I", 0)) # flags
        f.write(struct.pack("<fff", loc.x, loc.y, loc.z)) # position
        f.write(struct.pack("<ffff", -rot.w, rot.x, rot.y, rot.z)) # orientation
        f.write(struct.pack("<f", 1.0)) # zoom/scale
        f.write(struct.pack("<ii", scene.render.resolution_x, scene.render.resolution_y))
        f.write(b'\x00' * 36) # Empty float arrays (3x3 matrix etc)
        f.write(struct.pack("<ff", 0.0, 0.0)) # CA
        f.write(struct.pack("<f", 100.0)) # ISO
        f.write(b'\x00' * 48) # Padding
        f.write(struct.pack("<f", (cam_obj.data.dof.focus_distance if cam_obj.data.dof.use_dof else 10.0)))
        f.write(struct.pack("<fff", 1.0, 1.0, cam_obj.data.sensor_width))
        f.write(struct.pack("<iffi", 8, cam_obj.data.lens / 1000.0, 5.0, 14))

class ExportNRA2(bpy.types.Operator, ExportHelper):
    bl_idname = "export_scene.nra2_final"
    bl_label = "Export NRA2 Scene"
    filename_ext = ".nra2"

    def execute(self, context):
        filepath = self.filepath
        basepath = os.path.dirname(filepath)
        objs = [o for o in context.scene.objects if o.type == 'MESH' and o.visible_get()]
        
        node_idx = 0
        mat_map = {}
        mat_lines = []
        
        for obj in objs:
            m = obj.data.materials[0] if obj.data.materials else None
            name = m.name if m else "Default"
            if name not in mat_map:
                col = m.node_tree.nodes.get("Principled BSDF").inputs[0].default_value[:3] if (m and m.use_nodes) else (0.7,0.7,0.7)
                mat_lines.append(f"diffuse # {node_idx}")
                mat_lines.append(f"color d {col[0]:.4f} {col[1]:.4f} {col[2]:.4f} # {node_idx+1}")
                mat_lines.append(f"mult 1 {node_idx+1} 0 # {node_idx+2}")
                mat_map[name] = node_idx + 2
                node_idx += 3

        obj_lines = []
        for obj in objs:
            write_geo(obj, basepath)
            midx = mat_map.get(obj.data.materials[0].name if obj.data.materials else "Default", 0)
            obj_lines.append(f"{midx} {obj.name}")

        with open(filepath, "w") as f:
            f.write(f"black\n{node_idx}\n" + "\n".join(mat_lines) + f"\n\n{len(obj_lines)}\n" + "\n".join(obj_lines))

        cam = next((o for o in context.scene.objects if o.type == 'CAMERA'), None)
        if cam: write_camera(filepath, cam)
        
        return {'FINISHED'}

def register():
    bpy.utils.register_class(ExportNRA2)
    bpy.types.TOPBAR_MT_file_export.append(lambda self, context: self.layout.operator(ExportNRA2.bl_idname, text="NRA2 (.nra2)"))

if __name__ == "__main__":
    register()