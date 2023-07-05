import math
import numpy as np

colorList = [
    "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
    "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
    "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
    "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"
]

def appendToDict(d:dict, key, value):
    if key not in d:
        d[key] = value
    elif type(d[key]) == list:
        d[key].append(value)
    else:
        d[key] = [d[key], value]

def testCenter(C_x, C_y, x1, y1, x2, y2, r):
    testFirstPt = math.isclose((x1 - C_x) ** 2 + (y1 - C_y) ** 2, r ** 2)
    testSecondPt = math.isclose((x2 - C_x) ** 2 + (y2 - C_y) ** 2, r ** 2)

    if testFirstPt and testSecondPt:
        print(f'Passed Test for Center: ({C_x}, {C_y})')
    else:
        print(f'Failed Test for Center: ({C_x}, {C_y})')

def findCenter(x1, y1, x2, y2, r):
    p1 = np.array([x1, y1])
    p2 = np.array([x2, y2])

    # Calculate midpoint between two points
    midPoint = (p1 + p2) / 2

    # Calculate normalized direction vector perpendicular to line segment
    v = np.array([p2[1] - p1[1], p1[0] - p2[0]])
    v = v / np.linalg.norm(v)

    # Calculate distance between midpoint and center
    dist = np.sqrt(r ** 2 - (np.linalg.norm(p1 - p2) / 2) ** 2)

    # Calculate centers of circles
    C1 = midPoint + dist * v
    C2 = midPoint - dist * v

    # Test the centers
    testCenter(C1[0], C1[1], x1, y1, x2, y2, r)
    testCenter(C2[0], C2[1], x1, y1, x2, y2, r)

    return C1, C2

def getIntersection(x1, y1, r1, x2, y2, r2):
    d = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    a = (r1 ** 2 - r2 ** 2 + d ** 2) / (2 * d)
    h = math.sqrt(r1 ** 2 - a ** 2)

    x3 = x1 + a * (x2 - x1) / d
    y3 = y1 + a * (y2 - y1) / d

    i_x1 = x3 + h * (y2 - y1) / d
    i_y1 = y3 - h * (x2 - x1) / d

    i_x2 = x3 - h * (y2 - y1) / d
    i_y2 = y3 + h * (x2 - x1) / d

    return (i_x1, i_y1), (i_x2, i_y2)

def determineWhichIntersection(i_x1, i_y1, i_x2, i_y2, firstLayerBound):
    i1_angle_wrt_org = math.atan2(i_y1, i_x1)
    i1_angle_wrt_org = convertNegRadian(i1_angle_wrt_org)

    i2_angle_wrt_org = math.atan2(i_y2, i_x2)
    i2_angle_wrt_org = convertNegRadian(i2_angle_wrt_org)

    if i1_angle_wrt_org > math.pi:
        if firstLayerBound == 0:
            firstLayerBound = 2 * math.pi
    else:
        if firstLayerBound == 2 * math.pi:
            firstLayerBound = 0
    
    i1_wedge_angle = min(abs(i1_angle_wrt_org - firstLayerBound), abs(firstLayerBound - i1_angle_wrt_org))

    if i2_angle_wrt_org > math.pi:
        if firstLayerBound == 0:
            firstLayerBound = 2 * math.pi
    else:
        if firstLayerBound == 2 * math.pi:
            firstLayerBound = 0
    
    i2_wedge_angle = min(abs(i2_angle_wrt_org - firstLayerBound), abs(firstLayerBound - i2_angle_wrt_org))

    # print(f'I1 Wedge Angle: {i1_wedge_angle}, I2 Wedge Angle: {i2_wedge_angle}')

    if i1_wedge_angle < i2_wedge_angle:
        return i1_angle_wrt_org
    else:
        return i2_angle_wrt_org

# Converts negative radians to positive
def convertNegRadian(angle):
    return angle if angle >= 0 else 2 * math.pi - abs(angle)