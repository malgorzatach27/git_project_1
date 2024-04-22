

from math import sin, cos, sqrt, atan, atan2, degrees, radians
import sys
import numpy as np
#from pprint

o = object()

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2
    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 
        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
            
    def plh2xyz(self, phi, lam, h):
        phi = radians(phi)
        lam = radians(lam)
        Rn = self.a /sqrt(1 - self.ecc2 * sin(phi)**2)
        q = Rn * self.ecc2 * sin(phi)
        x = (Rn + h) * cos(phi) * cos(lam)
        y = (Rn + h) * cos(phi) * sin(lam)
        z = (Rn + h) * sin(phi) - q
        return x, y, z
    
    #sama napisalam
    def xyz2neu(self, phi, lam):
        phi = radians(phi)
        lam = radians(lam)
        R = np.array([-sin(lam), -sin(phi)*cos(lam), cos(phi)*cos(lam)], 
                     [cos(lam), -sin(phi)*sin(lam), cos(phi)*sin(lam)],
                     [0, cos(phi), sin(phi)])
    #podac srodek ukladu wsp
    
    
    #def bl2pl2000:
        
    #def bl2pl1992:
    

#z 3 miejscami dziesiętnymi
#transhol? żeby sprawdzic czy dziala poprawnie - czy wspolrzedne dobre wychodzą

if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    # dane XYZ geocentryczne
    # X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    # phi, lam, h = geo.xyz2plh(X, Y, Z)
    # print(phi, lam, h)
    # # phi, lam, h = geo.xyz2plh2(X, Y, Z)
    # # print(phi, lam, h)
    print(sys.argv)
    input_file_path = sys.argv[-1]
    if  '--xyz2plh' in sys.argv and '--plh2xyz' in sys.argv:
        print("mozesz podac tylko jedną flagę")
    elif '--xyz2plh' in sys.argv:
        with open(input_file_path, 'r')as f:
            lines = f.readlines()
            coords_lines = lines[4:]
    
            coords_plh = []
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                x_str, y_str, z_str = coord_line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                phi, lam, h = geo.xyz2plh(x, y, z)
                coords_plh.append([phi, lam, h])
    
        with open('result_xyz2plh.txt', 'w')as f:
            f.write('phi[deg], lam[deg], h[m] \n')
            for coords_list in coords_plh:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')

    elif '--plh2xyz' in sys.argv:        
        with open(input_file_path, 'r')as f:
            lines = f.readlines()
            coard_lines = lines[1:]
    
            coords_xyz = []
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                phi_str, lam_str, h_str = coord_line.split(',')
                phi, lam, h = (float(phi_str), float(lam_str), float(h_str))
                x, y, z = geo.plh2xyz(phi, lam, h)
                coords_xyz.append([x, y, z])
        
            # pprint(coords_xyz)
    
        with open('result_plh2xyz.txt', 'w')as f:
            f.write('x[m], y[m], z[m] \n')
            for coords_list in coords_xyz:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')
        

    
