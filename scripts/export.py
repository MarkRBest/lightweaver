import os
import bpy
import math

def mult(vec, matrix) :
    x = vec[0] * matrix[0][0] + vec[1] * matrix[0][1] + vec[2] * matrix[0][2]
    y = vec[0] * matrix[1][0] + vec[1] * matrix[1][1] + vec[2] * matrix[1][2]
    z = vec[0] * matrix[2][0] + vec[1] * matrix[2][1] + vec[2] * matrix[2][2]
    return((x,y,z))

def writeLight(f, item) :
    #vert = rst(item.rotation_euler, [1.0, 1.0, 1.0], [0.0, 0.0, 0.0], item.location)
    vert = item.matrix_world.to_translation()
    f.write("Light ")
    f.write("x: " + str(vert[0]) +" ")
    f.write("y: " + str(vert[1]) +" ")
    f.write("z: " + str(vert[2]) +" ")
    f.write("\n")
    f.flush()

def writeCamera(f, item) :
    #vert = rst(item.rotation_euler, [1.0, 1.0, 1.0], [0.0, 0.0, 0.0], item.location)
    pos = item.matrix_world.to_translation()
    dir = item.matrix_world.to_euler()
    f.write("Camera ")
    f.write("angle_x: "+ str(item.data.angle_x) +" ")
    f.write("angle_y: "+ str(item.data.angle_y) +" ")
    #f.write("angle_z: "+ str(item.data.angle_z) +" ")
    f.write("x: " + str(pos[0]) +" ")
    f.write("y: " + str(pos[1]) +" ")
    f.write("z: " + str(pos[2]) +" ")
    f.write("\n")
    f.flush()

def rst(rot, scale, trans, vert) :
    # rotate vertex
    cr = (math.cos(rot[1]),math.cos(rot[2]),math.cos(-rot[0]))
    sr = (math.sin(rot[1]),math.sin(rot[2]),math.sin(-rot[0]))
    
    x_rot = ((cr[0],sr[0],0),(-sr[0],cr[0],0),(0,0,1))
    y_rot = ((cr[1],0,-sr[1]),(0,1,0),(sr[1],0,cr[1]))
    z_rot = ((1,0,0),(0,cr[2],sr[2]),(0,-sr[2],cr[2]))
    
    #vert = mult(vert, x_rot)
    #vert = mult(vert, y_rot)
    #vert = mult(vert, z_rot)
    # scale and translate
    vx = (vert[0] * scale[0]) + trans[0]
    vy = (vert[1] * scale[1]) + trans[1]
    vz = (vert[2] * scale[2]) + trans[2]
    return (vx, vy, vz)

def writeVertex(f, uv, normal, vert):
    # write vertex info
    f.write("vert: ")
    f.write("vx: " + str(vert[0]) +" ")
    f.write("vy: " + str(vert[1]) +" ")
    f.write("vz: " + str(vert[2]) +" ")
    # write normals
    f.write("nx: " + str(normal[0]) +" ")
    f.write("ny: " + str(normal[1]) +" ")
    f.write("nz: " + str(normal[2]) +" ")
    # write uv coords
    f.write("tu: " + str(uv[0]) +" ")
    f.write("tv: " + str(uv[1]) +" ")
    f.write("\n")
    f.flush()

def output(f) :
    for item in bpy.data.objects:

        print(item.name, item.type)
        
        if item.type == 'LAMP':
            print (item.location)
            r=item.rotation_axis_angle
            print("rot axis: ", r[0]," ",r[1]," ",r[2]," ",r[3])
            print("rot: ", item.rotation_euler)
            writeLight(f, item)
        if item.type == 'CAMERA':
            #print (item.location)
            #print (item.data.angle_x)
            #print (item.data.angle_y)
            writeCamera(f, item)
        if item.type == 'MESH':
            count = 0
            #print(item.name)
            print("location: ", item.location)
            print("scale: ", item.scale)
            r=item.rotation_axis_angle
            print("rot axis: ", r[0]," ",r[1]," ",r[2]," ",r[3])
            print("rot: ", item.rotation_euler)
            loc = item.location
            scale = item.scale
            rot = item.rotation_euler
            f.write("Mesh "+ item.name + "\n")
            for face in item.data.polygons:
                verts_in_face = face.vertices[:]
                n = face.normal
                UV_ENABLED = False
                uv = (0,0);
                if len(item.data.uv_layers) > 0:
                    uv_layer = item.data.uv_layers[0]
                    UV_ENABLED = True
                #    print(face.loop_indices[i], " " , i);                    
                #    print(uv)
                #i+=1
                if len(verts_in_face) == 3: # write triangle
                    f.write("Tri: " + str(count) +"\n")
                    if UV_ENABLED :
                        uv = uv_layer.data[face.loop_indices[0]].uv;
                    #vert = rst(rot,scale,loc,item.data.vertices[verts_in_face[0]].co)
                    
                    vert = item.matrix_world * item.data.vertices[verts_in_face[0]].co
                    writeVertex(f, uv, n, vert)
                    if UV_ENABLED :
                        uv = uv_layer.data[face.loop_indices[1]].uv;
                    #vert = rst(rot,scale,loc,item.data.vertices[verts_in_face[1]].co)
                    vert = item.matrix_world * item.data.vertices[verts_in_face[1]].co
                    writeVertex(f, uv, n, vert)
                    if UV_ENABLED :
                        uv = uv_layer.data[face.loop_indices[2]].uv;
                    #vert = rst(rot,scale,loc,item.data.vertices[verts_in_face[2]].co)
                    vert = item.matrix_world * item.data.vertices[verts_in_face[2]].co
                    writeVertex(f, uv, n, vert)
                    count = count + 1;
                elif len(verts_in_face) == 4:  # convert square to two triangles
                    # first tri
                    f.write("Tri: " + str(count) +"\n")
                    if UV_ENABLED :
                        uv = uv_layer.data[face.loop_indices[0]].uv;
                    #vert = rst(rot,scale,loc,item.data.vertices[verts_in_face[0]].co)
                    vert = item.matrix_world * item.data.vertices[verts_in_face[0]].co
                    writeVertex(f, uv, n, vert)
                    if UV_ENABLED :
                        uv = uv_layer.data[face.loop_indices[1]].uv;
#                     vert = rst(rot,scale,loc,item.data.vertices[verts_in_face[1]].co)
                    vert = item.matrix_world * item.data.vertices[verts_in_face[1]].co
                    writeVertex(f, uv, n, vert)
                    if UV_ENABLED :
                        uv = uv_layer.data[face.loop_indices[2]].uv;
                    #vert = rst(rot,scale,loc, item.data.vertices[verts_in_face[2]].co)
                    vert = item.matrix_world * item.data.vertices[verts_in_face[2]].co
                    writeVertex(f, uv, n, vert)
                    count = count + 1;
                    # second tri
                    f.write("Tri: " + str(count) +"\n")
                    if UV_ENABLED :
                        uv = uv_layer.data[face.loop_indices[0]].uv;
#                     vert = rst(rot,scale,loc,item.data.vertices[verts_in_face[0]].co)
                    vert = item.matrix_world * item.data.vertices[verts_in_face[0]].co
                    writeVertex(f, uv, n, vert)
                    if UV_ENABLED :
                        uv = uv_layer.data[face.loop_indices[3]].uv;
#                     vert = rst(rot,scale,loc,item.data.vertices[verts_in_face[3]].co)
                    vert = item.matrix_world * item.data.vertices[verts_in_face[3]].co
                    writeVertex(f, uv, n, vert)
                    if UV_ENABLED :
                        uv = uv_layer.data[face.loop_indices[2]].uv;
#                     vert = rst(rot,scale,loc,item.data.vertices[verts_in_face[2]].co)
                    vert = item.matrix_world * item.data.vertices[verts_in_face[2]].co
                    writeVertex(f, uv, n, vert)
                    count = count + 1;
            f.write("MeshEnd\n")


def export():
    f = "/home/markbest/lightweaver/data/stairs.txt"
    f = open(f, "w+")
    f.truncate()
    f.write("Sceen\n");
    output(f)
    f.write("\n")
    f.flush()
    f.close()
    print ("sceen export complete")

    
export()