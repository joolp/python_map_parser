import numpy as np
import os
import sys
from PIL import Image

EPS = 1e-9

TEXPATH = 'textures/'
TEXEXT  = '.jpg'

MAPPATH = 'maps/'
MAPEXT  = '.map' 

def parse_file(name):
    entities = []

    f = open(name)
    lines = f.readlines()
    level = 0

    entity = {'prop':{}, 'geo':[]}
    brush  = []

    for line in lines:
        if line.startswith('//'):
            continue
        if line.startswith('{'):
            level += 1
            if level == 1:
                entity = {'prop':{}, 'geo':[]}
            if level == 2:
                brush  = []
            continue
        if line.startswith('}'):
            level -= 1
            if level == 0:
                entities.append(entity)
            if level == 1:
                entity['geo'].append(brush)
            continue

        if level == 1:
            tokens = line.split('\"')
            # print(tokens)
            entity['prop'][tokens[1]] = tokens[3]

        if level == 2:
            tokens = line.split()
            # print(tokens)
            plane = {}
            if (len(tokens) > 30):
                # print('Map format is Valve')
                global mapformat
                mapformat = 'valve'

                plane['x1'] = tokens[1]
                plane['y1'] = tokens[2]
                plane['z1'] = tokens[3]

                plane['x2'] = tokens[6]
                plane['y2'] = tokens[7]
                plane['z2'] = tokens[8]

                plane['x3'] = tokens[11]
                plane['y3'] = tokens[12]
                plane['z3'] = tokens[13]

                plane['tn'] = tokens[15]

                plane['ux'] = tokens[17]
                plane['uy'] = tokens[18]
                plane['uz'] = tokens[19]
                plane['uo'] = tokens[20]

                plane['vx'] = tokens[23]
                plane['vy'] = tokens[24]
                plane['vz'] = tokens[25]
                plane['vo'] = tokens[26]

                plane['ro'] = tokens[28]
                plane['sx'] = tokens[29]
                plane['sy'] = tokens[30]
            elif(len(tokens) < 30):
                # print('Map format is Standard')
                mapformat = 'standard'

                plane['x1'] = tokens[1]
                plane['y1'] = tokens[2]
                plane['z1'] = tokens[3]

                plane['x2'] = tokens[6]
                plane['y2'] = tokens[7]
                plane['z2'] = tokens[8]

                plane['x3'] = tokens[11]
                plane['y3'] = tokens[12]
                plane['z3'] = tokens[13]

                plane['tn'] = tokens[15]

                plane['ux'] = tokens[00]
                plane['uy'] = tokens[00]
                plane['uz'] = tokens[00]
                plane['uo'] = tokens[16]

                plane['vx'] = tokens[00]
                plane['vy'] = tokens[00]
                plane['vz'] = tokens[00]
                plane['vo'] = tokens[17]

                plane['ro'] = tokens[18]
                plane['sx'] = tokens[19]
                plane['sy'] = tokens[20]

            brush.append(plane)
    return entities

def calculate_uv(plane, position):

    file = TEXPATH + plane['tn'] + TEXEXT
    tex = None

    if (os.path.exists(file)):
        tex = Image.open(file)
    else:
        return [0.0, 0.0]

    nu = nv = np.array([])

    if mapformat == 'valve':
        nu = np.array([float(plane['ux']),float(plane['uy']),float(plane['uz'])])
        nv = np.array([float(plane['vx']),float(plane['vy']),float(plane['vz'])])
    elif mapformat == 'standard':
        a = np.array([float(plane['x1']),float(plane['y1']),float(plane['z1'])])
        b = np.array([float(plane['x2']),float(plane['y2']),float(plane['z2'])])
        c = np.array([float(plane['x3']),float(plane['y3']),float(plane['z3'])])
        n = np.cross(a-c,b-c)
        n /= np.linalg.norm(n)

        nv = np.array([0.0, 0.0, 1.0])
        nu = np.cross(nv, n)

    uo = float(plane['uo'])
    vo = float(plane['vo'])
    sx = float(plane['sx'])
    sy = float(plane['sy'])
    w  = float(tex.size[0])
    h  = float(tex.size[1])

    tu = (np.dot(position,nu)/sx + uo)/w
    tv = (np.dot(position,nv)/sy + vo)/h

    return [tu, tv]

def sort_vertices(poly):
    tmppoly = poly
    n = len(tmppoly)

    center = np.array([0.,0.,0.])
    for i in tmppoly:
        center += i[0]
    center /= float(n)

    for i in range(n-2):
        a = tmppoly[i][0] - center
        a /= np.linalg.norm(a)

        nm = np.cross(a, tmppoly[i][1])

        smallestAngle = -1.0
        smallest = -1

        for j in range(i+1, n):
            if np.dot(poly[j][0] - center, nm) < 0:
                b = poly[j][0] - center
                b /= np.linalg.norm(b)
                angle = np.dot(a,b)

                if angle > smallestAngle:
                    smallestAngle = angle
                    smallest = j

        tmp = tmppoly[i+1]
        tmppoly[i+1]=tmppoly[smallest]
        tmppoly[smallest] = tmp

    return tmppoly

def build_geo(brush):
    geos     = []
    vertices = []
    indices  = []
    eqs      = []
    poly     = []

    for plane in brush:
        a = np.array([float(plane['x1']),float(plane['y1']),float(plane['z1'])])
        b = np.array([float(plane['x2']),float(plane['y2']),float(plane['z2'])])
        c = np.array([float(plane['x3']),float(plane['y3']),float(plane['z3'])])
        n = np.cross(a-c,b-c)
        n /= np.linalg.norm(n)
        d = -np.dot(n,a)
        eqs.append([n,d])
        poly.append([])

    n = len(eqs)
    for i in range(n-2):
        for j in range(i+1,n-1):
            for k in range(j+1,n):

                if i==j or j==k or i==k:
                    continue
                den = np.dot(eqs[i][0], np.cross(eqs[j][0], eqs[k][0]))
                if (-EPS < den and den < EPS):
                    continue

                poi = ( -eqs[i][1]*np.cross(eqs[j][0],eqs[k][0]) \
                        -eqs[j][1]*np.cross(eqs[k][0],eqs[i][0]) \
                        -eqs[k][1]*np.cross(eqs[i][0],eqs[j][0]) ) / den

                ispoin = True
                for m in eqs:
                    if np.dot(poi, m[0]) + m[1] < 0:
                        ispoin = False
                        break

                if not ispoin:
                    continue

                poly[i].append([poi, eqs[i][0], calculate_uv(brush[i], poi)])
                poly[j].append([poi, eqs[j][0], calculate_uv(brush[j], poi)])
                poly[k].append([poi, eqs[k][0], calculate_uv(brush[k], poi)])

    for p in poly:
        p = sort_vertices(p)

    for i in range(n):
        geos.append([brush[i]['tn'], poly[i]])

    return geos

def write_to_poly_list(name, geos):
    o = open(name, 'w')

    for i in range(len(geos)):
        o.write(f'{i}:\n')
        for t in geos[i]:
            o.write(f'\t{t[0]}:\n')
            for f in t[1]:
                o.write(f'\t\t{f[0][0]} {f[0][1]} {f[0][2]} {f[1][0]} {f[1][1]} {f[1][2]} {f[2][0]} {f[2][1]}\n')

def write_vertices(name, geos):
    o = open(name, 'w')
    for i in range(len(geos)):
        for t in geos[i]:
            for f in t[1]:
                o.write(f'({f[0][0]}, {f[0][1]}, {f[0][2]})\n')

def write_wavefront(name, geos):
    o = open(name, 'w')

    def issame(v1, v2):
        if abs(v1[0]-v2[0])<EPS and abs(v1[1]-v2[1])<EPS and abs(v1[2]-v2[2])<EPS:
            return True
        return False

    def issame2(v1, v2):
        if abs(v1[0]-v2[0])<EPS and abs(v1[1]-v2[1])<EPS:
            return True
        return False

    vertex_list = []
    normal_list = []
    texcoord_list = []

    index_list = []
    index_normal_list = []
    index_texcoord_list = []

    for i in range(len(geos)):
        for t in geos[i]:
            for f in t[1]:
                taken = False
                for vv in vertex_list:
                    if issame(vv, f[0]):
                        taken = True
                        break
                if not taken:
                    vertex_list.append(f[0])

    for i in range(len(geos)):
        for t in geos[i]:
            for f in t[1]:
                taken = False
                for vv in normal_list:
                    if issame(vv, f[1]):
                        taken = True
                        break
                if not taken:
                    normal_list.append(f[1])

    def fract(v):
        return np.array([v[0] - int(v[0]), v[1] - int(v[1])])


    for i in range(len(geos)):
        for t in geos[i]:
            for f in t[1]:
                taken = False
                for vv in texcoord_list:
                    if issame2(vv, f[2]):
                        taken = True
                        break
                if not taken:
                    texcoord_list.append(fract(f[2]))


    for i in range(len(geos)):
        for t in geos[i]:
            # nn = len(t[1])
            pvi = []
            pni = []
            pti = []
            for f in t[1]:
                for i in range(len(vertex_list)):
                    if issame(f[0], vertex_list[i]):
                        pvi.append(i+1)
                # print(pvi)
                if len(pvi) == 4:
                    pvi.append(pvi[3])
                    pvi.append(pvi[3])
                    pvi[3] = pvi[0]
                    pvi[4] = pvi[2]

                for i in range(len(normal_list)):
                    if issame(f[1], normal_list[i]):
                        pni.append(i+1)

                if len(pni) == 4:
                    pni.append(pni[3])
                    pni.append(pni[3])
                    pni[3] = pni[0]
                    pni[4] = pni[2]

                for i in range(len(texcoord_list)):
                    if issame2(f[2], texcoord_list[i]):
                        pti.append(i+1)

                if len(pti) == 4:
                    pti.append(pti[3])
                    pti.append(pti[3])
                    pti[3] = pti[0]
                    pti[4] = pti[2]

                for pp in pvi:
                    index_list.append(pp)

                for pn in pni:
                    index_normal_list.append(pn)

                for pt in pti:
                    index_texcoord_list.append(pt)

    #print(index_list)

    for v in vertex_list:
        o.write(f'v {v[0]} {v[1]} {v[2]}\n')

    for v in texcoord_list:
        o.write(f'vt {v[0]} {v[1]}\n')

    for v in normal_list:
        o.write(f'vn {v[0]} {v[1]} {v[2]}\n')

    for i in range(0, len(index_list), 3):
        o.write('f ')
        o.write(f'{index_list[i]}/{index_texcoord_list[i]}/{index_normal_list[i]} {index_list[i+1]}/{index_texcoord_list[i+1]}/{index_normal_list[i+1]} {index_list[i+2]}/{index_texcoord_list[i+2]}/{index_normal_list[i+2]}\n')

def main():
    mapfile = 'E3M4'

    ents = parse_file(MAPPATH + mapfile + MAPEXT)
    geos = []

    for i in ents[0]['geo']:
        geos.append(build_geo(i))

    write_to_poly_list(f'polylist_{mapfile}.txt', geos)

    write_vertices('vt.txt', geos)

    write_wavefront('vt.obj', geos)

main()