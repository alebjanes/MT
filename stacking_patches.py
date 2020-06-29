import gdal
import numpy as np
from osgeo import gdal_array
from sklearn.feature_extraction import image
import itertools, random

base_tif=gdal.Open("/Users/alexandrabjanes/MT/Datos500m/normalized/DistanceToUA.tif")

#SETTINGS
outdir = "/Users/alexandrabjanes/MT/Samples7x7"
indir = "/Users/alexandrabjanes/MT/Datos500m/Normalized/"
year=2018
yearmonth=201812
patchsize=7

tmax = gdal_array.LoadFile(indir +str(year)+ "/" + str(yearmonth) + "tx.tif")
tmin = gdal_array.LoadFile(indir +str(year)+ "/" + str(yearmonth) + "tn.tif")
windspeed = gdal_array.LoadFile(indir +str(year)+ "/" + str(yearmonth) + "ws.tif")
prec = gdal_array.LoadFile(indir +str(year)+ "/" + str(yearmonth) + "rn.tif")
aet = gdal_array.LoadFile(indir +str(year)+ "/" + str(yearmonth) + "aet.tif")
cwd = gdal_array.LoadFile(indir +str(year)+ "/" + str(yearmonth) + "cwd.tif")
NDVI = gdal_array.LoadFile(indir +str(year)+ "/" + str(yearmonth) + "NDVI.tif")
aspect = gdal_array.LoadFile(indir + "Aspect.tif")
slope = gdal_array.LoadFile(indir + "Slope.tif")
elev = gdal_array.LoadFile(indir + "Elevation.tif")
surface = gdal_array.LoadFile(indir + "SurfaceRoughness.tif")
distroads = gdal_array.LoadFile(indir + "DistanceToRoads.tif")
distrivers = gdal_array.LoadFile(indir + "Distancetowater.tif")
distUA = gdal_array.LoadFile(indir + "DistanceToUA.tif")

ignition_class = gdal_array.LoadFile("/Users/alexandrabjanes/MT/IncendiosBiobio/Rasters/clipped/IF" + str(yearmonth) + ".sdat")
#IRD = gdal_array.LoadFile("/Users/alexandrabjanes/MT/Datos500m/IRD/IRD" + str(yearmonth) + ".tif")

Stack=np.dstack((ignition_class,tmax,tmin,windspeed,prec,aet, cwd, NDVI, aspect, slope, elev, surface, distroads, distrivers, distUA))
#Stack=np.dstack((IRD,tmax,tmin,windspeed,prec,aet, cwd, NDVI, aspect, slope, elev, surface, distroads, distrivers, distUA))

print(Stack.shape)

i=0
j=0
n=1
covered = np.zeros((550,480), dtype=int)

def GeoTiff1banda(Name, Array, driver, NDV, DataType, base):
    Array[np.isnan(Array)] = NDV
    driver = gdal.GetDriverByName(driver)
    DataSet = driver.Create(Name, 480, 550, 1, eType=gdal.GDT_Float32)
    Geotr = bas.GetGeoTransform()
    proj = bas.GetProjection()
    DataSet.SetGeoTransform(Geotr)
    DataSet.SetProjection(proj)
    DataSet.GetRasterBand(1).WriteArray(Array)
    DataSet.GetRasterBand(1).SetNoDataValue(NDV)
    DataSet.FlushCache()
    return Name

def select(size, pair_size):
    g =itertools.permutations(range(size),pair_size)
    alist = list(g)
    random.shuffle(alist)
    return alist

#busqueda muestra clase 0
def search(Matriz_inicial, x, y, cover, patchsize, r):
    centroid = int(patchsize/2 - 0.5)
    if r == 0:
        if x<patchsize: a=patchsize
        elif x>550-2*patchsize: a=550-2*patchsize
        else: a=x
        if y<patchsize: b=patchsize
        elif y>480-2*patchsize: b=480-2*patchsize
        else: b=y
        area1 = Matriz_inicial[a-patchsize:a+(2*patchsize), b-patchsize:b+(2*patchsize), :]
        pairs=select(area1.shape[0]-patchsize+1, 2)
        for t in range(len(pairs)):
            i=pairs[t][0]
            j=pairs[t][1]
            randompatch = area1[i:i+patchsize, j:j+patchsize,:]
            gt = area1[i:i+patchsize, j:j+patchsize,0]
            countzero = np.count_nonzero(gt)
            if countzero == 0 and cover[a-patchsize+i+centroid, b-patchsize+j+centroid] == 0:
                patch=np.squeeze(randompatch)
                cover[a-patchsize+i+centroid, b-patchsize+j+centroid] = 2
                return patch, cover
            elif t == (len(pairs)-1):
                r = 1
                break
    if r == 1:
        if x<2*patchsize: a=2*patchsize
        elif x>550-3*patchsize: a=550-3*patchsize
        else: a=x
        if y<2*patchsize: b=2*patchsize
        elif y>480-3*patchsize: b=480-3*patchsize
        else: b=y
        area2 = Matriz_inicial[a-(2*patchsize):a+(3*patchsize), b-(2*patchsize):b+(3*patchsize), :]
        pairs=select(area2.shape[0] - patchsize + 1, 2)
        for t in range(len(pairs)):
            i=pairs[t][0]
            j=pairs[t][1]
            randompatch = area2[i:i+patchsize, j:j+patchsize,:]
            gt = area2[i:i+patchsize, j:j+patchsize,0]
            countzero = np.count_nonzero(gt)
            if countzero == 0 and cover[a-2*patchsize+i+centroid, b-2*patchsize+j+centroid] == 0:
                patch=np.squeeze(randompatch)
                cover[a-2*patchsize+i+centroid, b-2*patchsize+j+centroid] = 2
                return patch, cover
            elif t == (len(pairs)-1):
                r = 2
                break
    if r == 2:
        if x<3*patchsize: a=3*patchsize
        elif x>550-4*patchsize: a=550-4*patchsize
        else: a=x
        if y<3*patchsize: b=3*patchsize
        elif y>480-4*patchsize: b=480-4*patchsize
        else: b=y
        area3 = Matriz_inicial[a-(3*patchsize):a+(4*patchsize), b-(3*patchsize):b+(4*patchsize), :]
        pairs=select(area3.shape[0]-patchsize+1, 2)
        for t in range(len(pairs)):
            i=pairs[t][0]
            j=pairs[t][1]
            randompatch = area3[i:i+patchsize, j:j+patchsize,:]
            gt = area3[i:i+patchsize, j:j+patchsize,0]
            countzero = np.count_nonzero(gt)
            if countzero==0 and cover[a-3*patchsize+i+centroid, b-3*patchsize+j+centroid] == 0:
                patch=np.squeeze(randompatch)
                cover[a-3*patchsize+i+centroid, b-3*patchsize+j+centroid] = 2
                return patch, cover
            elif t == (len(pairs)-1):
                r = 3
                break
    if r == 3:
        if x<4*patchsize: a=4*patchsize
        elif x>550-5*patchsize: a=550-5*patchsize
        else: a=x
        if y<4*patchsize: b=4*patchsize
        elif y>480-5*patchsize: b=480-5*patchsize
        else: b=y
        area4 = Matriz_inicial[a-(4*patchsize):a+(5*patchsize), b-(4*patchsize):b+(5*patchsize), :]
        pairs=select(area4.shape[0]-patchsize+1, 2)
        for t in range(len(pairs)):
            i=pairs[t][0]
            j=pairs[t][1]
            randompatch = area4[i:i+patchsize, j:j+patchsize,:]
            gt = area4[i:i+patchsize, j:j+patchsize,0]
            countzero = np.count_nonzero(gt)
            if countzero==0 and cover[a-4*patchsize+i+centroid, b-4*patchsize+j+centroid] == 0:
                patch=np.squeeze(randompatch)
                cover[a-4*patchsize+i+centroid, b-4*patchsize+j+centroid] = 2
                return patch, cover
            elif t == (len(pairs)-1):
                r=4
                print('r=4')
                break
    if r == 4:
        if x<5*patchsize: a=5*patchsize
        elif x>550-6*patchsize: a=550-6*patchsize
        else: a=x
        if y<5*patchsize: b=5*patchsize
        elif y>480-6*patchsize: b=480-6*patchsize
        else: b=y
        area5 = Matriz_inicial[a-(5*patchsize):a+(6*patchsize), b-(5*patchsize):b+(6*patchsize), :]
        print(area5.shape)
        pairs=select(area5.shape[0]-patchsize+1, 2)
        for t in range(len(pairs)):
            i=pairs[t][0]
            j=pairs[t][1]
            randompatch = area5[i:i+patchsize, j:j+patchsize,:]
            gt = area5[i:i+patchsize, j:j+patchsize,0]
            countzero = np.count_nonzero(gt)
            if countzero==0 and cover[a-5*patchsize+i+centroid, b-5*patchsize+j+centroid] == 0:
                patch=np.squeeze(randompatch)
                #patch0=patch.transpose(2,0,1)
                cover[a-5*patchsize+i+centroid, b-5*patchsize+j+centroid] = 2
                return patch, cover
            elif t == (len(pairs)-1):
                print('r = 5')
                break    

nodata = -99999

for i in range(550-patchsize):
    for j in range(480-patchsize):
        sample = Stack[i:i+patchsize, j:j+patchsize, 0]
        centroid = int(patchsize/2 - 0.5)
        if sample[centroid,centroid] == 1 and not nodata in sample:
                patch = Stack[i:i+patchsize, j:j+patchsize, :]
                covered[i+centroid, j+centroid] = 1
                patch0, covered = search(Stack, i, j, covered, patchsize, 0)
                if n == 1:
                    Dataset = patch
                    Dataset = np.stack((Dataset, patch0), axis=0)
                else:
                    patch = np.expand_dims(patch, axis=0)
                    patch0 = np.expand_dims(patch0, axis=0)
                    Dataset = np.concatenate((Dataset, patch, patch0), axis=0)
                n+=1

print(Dataset.shape)
np.save(outdir + '/sample' + str(yearmonth) + '.npy', Dataset)

GeoTiff1banda('/Users/alexandrabjanes/MT/Cover_3/cover'+ str(yearmonth) +'.tif', covered, driver="GTiff", NDV=-99999, DataType=gdal.GDT_Float32, base=base_tif)
