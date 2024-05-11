from math import sin, cos, sqrt, tan, atan, atan2, degrees, radians
import sys
import numpy as np

o = object()

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        Modele:
        + WGS84
        + GRS80
        """
        if model == "wgs84":
            self.a = 6378137.0 
            self.b = 6356752.31424518
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model nie jest implementowany")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) #WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) 
    
    
    def deg2dms(self, dec_deg):
        """
        Funkcja pomocnicza zamieniająca stopnie dziesiętne 
        na stopnie, minuty i sekundy.

        Parameters
        ----------
        dec_degree : FLOAT
            stopnie dziesiętne

        Returns
        -------
        d : stopnie
        m : minuty
        s : sekundy
        """
        d = int(dec_deg)
        m_float = (dec_deg - d) * 60
        m = int(m_float)
        s = (m_float - m) * 60
        return d, m, f' {s:.5f}\"'


    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.  
        
        Parameters
        ----------
        X, Y, Z : FLOAT - współrzędne w układzie orto-kartezjańskim
             
        Returns
        -------
        phi : FLOAT [stopnie dziesiętne] - szerokość geodezyjna
        lam : FLOAT [stopnie dziesiętne] - długośc geodezyjna
        h   : FLOAT [metry] - wysokość elipsoidalna
            
        output [STR] - opcjonalnie, domyslnie
            dec_degree - stopnie dziesiętne - domyslnie
            dms - stopnie, minuty, sekundy 
        """
        r = sqrt(X**2 + Y**2)   
        phi_prev = atan(Z / (r * (1 - self.ecc2))) 
        phi = 0
        while abs(phi_prev - phi) > 0.000001/206265:    
            phi_prev = phi
            N = self.a / sqrt(1 - self.ecc2 * sin(phi_prev)**2)
            h = r / cos(phi_prev) - N
            phi = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lam = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(phi))**2);
        h = r / cos(phi) - N       
        if output == "dec_degree":
            return degrees(phi), degrees(lam), h 
        elif output == "dms":
            phi = self.deg2dms(degrees(phi))
            lam = self.deg2dms(degrees(lam))
            return f"{phi[0]:02d}:{phi[1]:02d}:{phi[2]:.2f}", f"{lam[0]:02d}:{lam[1]:02d}:{lam[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - format output nie jest zdefiniowany")
            
            
    def plh2xyz(self, phi, lam, h):
        """
        Funkcja transformacji współrzędnych geodezyjnych (phi, lam, h) 
        na współrzędne ortokartezjańskie (x, y, z).

        Parameters
        ----------
        phi : FLOAT [stopnie dziesiętne] - szerokość geodezyjna
        lam : FLOAT [stopnie dziesiętne] - długośc geodezyjna
        h   : FLOAT [metry] - wysokość elipsoidalna

        Returns
        -------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim

        """
        phi = radians(phi)
        lam = radians(lam)
        N = self.a /sqrt(1 - self.ecc2 * sin(phi)**2)
        x = (N + h) * cos(phi) * cos(lam)
        y = (N + h) * cos(phi) * sin(lam)
        z = (N * (1 - self.ecc2) + h) * sin(phi)
        return x, y, z


    def xyz2neu(self, x, y, z, x0, y0, z0):
        """
        Funkcja transformacji współrzędnych geocentrycznych punktu (x, y, z) 
        przy użyciu współrzędnych geocentrycznych srodka ukladu neu (x0, y0, z0)
        na współrzędne topocentryczne tego punktu (e, n, u).
        
        Parameters
        ----------
        x, y, z : FLOAT - współrzędne geocentryczne punktu
        x0, y0, z0 : FLOAT - współrzędne geocentryczne srodka ukladu (neu)
             
        Returns
        -------
        E, N, U : wspolrzedne topocentryczne punktu w ukladzie neu (easting, northing, up)

        """
        phi, lam, _= [radians(coord) for coord in self.xyz2plh(x0, y0, z0)]
        Rneu = np.array([[-sin(lam), -sin(phi)*cos(lam), cos(phi)*cos(lam)], 
                         [ cos(lam), -sin(phi)*sin(lam), cos(phi)*sin(lam)],
                         [        0,           cos(phi),          sin(phi)]])
        xyz_t = np.array([[x - x0],
                          [y - y0],
                          [z - z0]])
        [[E], [N], [U]] = Rneu.T @ xyz_t
        return E, N, U
        
    
    def sigma(self, phi):
        """
        Funkcja pomocnicza przy algorytmie bl2pl2000. 
        Oblicza długosc luku poludnika.

        Parameters
        ----------
        phi : FLOAT [stopnie dziesiętne] - szerokość geodezyjna

        Returns
        -------
        sigma - dlugosc luku poludnika

        """
        A0 = 1 - (self.ecc2/4) - (3*(self.ecc2)**2)/64 -  (5*(self.ecc2)**3)/256
        A2 = 3/8 * (self.ecc2 + (self.ecc2)**2/4 + 15*(self.ecc2)**3/128)
        A4 = 15/256 * ( (self.ecc2)**2 + (3*((self.ecc2)**3))/4 )
        A6 = 35 * (self.ecc2)**3 / 3072
        sigma = self.a * ( A0 * phi - A2 * np.sin(2*phi) + A4 * np.sin(4*phi) - A6 * np.sin(6*phi))
        return(sigma)
    
    
    def bl2pl2000(self, phi, lam):
        """
        Funkcja transformacji współrzędnych geodezyjnych (phi, lam, h) 
        na współrzędne w układzie pl-2000 (x2000, y2000).

        Parameters
        ----------
        phi : FLOAT [stopnie dziesiętne] - szerokość geodezyjna
        lam : FLOAT [stopnie dziesiętne] - długośc geodezyjna

        Returns
        -------
        x2000, y2000 : [FLOAT] - wspolrzedne w ukladzie pl-2000

        """
        m2000 = 0.999923
        phi = radians(phi)
        lam = radians(lam)
        if lam > radians(13.5) and lam < radians(16.5):
            lam0 = radians(15)
            nrstrefy = 5
        if lam > radians(16.5) and lam < radians(19.5):
            lam0 = radians(18)
            nrstrefy = 6
        if lam > radians(19.5) and lam < radians(22.5):
            lam0 = radians(21)
            nrstrefy = 7
        if lam > radians(22.5) and lam < radians(25.5):
            lam0 = radians(24)
            nrstrefy = 8
        b2 = (self.a**2)*(1 - self.ecc2)
        ecc22 = (self.a**2 - b2)/b2
        dlam = lam - lam0
        t = tan(phi)
        eta2 = ecc22 * ((cos(phi))**2)
        N = self.a/np.sqrt(1 - self.ecc2 * sin(phi)**2)
        sigma = self.sigma(phi)
        xgk= sigma + (dlam**2/2) * N * sin(phi)*cos(phi)*((1+(dlam**2/12)*(cos(phi))**2*(5-t**2+9*eta2+4*eta2**2)+(dlam**4/360)*cos(phi)**4*(61-58*t**2+t**4+270*eta2-330*eta2*t**2)))  
        ygk = dlam * N * cos(phi) * (1 + ((dlam**2)/6) * (cos(phi)**2) * (1 - t**2 + eta2) + (dlam**4/120) * (cos(phi)**4) * (5 - 18 * t**2 + t**4 + 14 * eta2 - 58 * eta2 * t**2))
                     
        x2000 = xgk * m2000
        y2000 = ygk * m2000 + nrstrefy * 1000000 + 500000
        return x2000, y2000
    
    def bl2pl1992(self, phi, lam):
        """
        Funkcja transformacji współrzędnych geodezyjnych (phi, lam, h) 
        na współrzędne w układzie pl-1992 (x1992, y1992).

        Parameters
        ----------
        phi : FLOAT [stopnie dziesiętne] - szerokość geodezyjna
        lam : FLOAT [stopnie dziesiętne] - długośc geodezyjna

        Returns
        -------
        x92, y92 :  [FLOAT] - współrzędne w układzie pl-1992

        """
        m92 = 0.9993
        phi = radians(phi)
        lam = radians(lam)
        lam0 = radians(19)
        b2 = (self.a**2)*(1 - self.ecc2)
        ecc22 = (self.a**2 - b2)/b2
        dlam = lam - lam0
        t = tan(phi)
        eta2 = ecc22 * ((cos(phi))**2)
        N = self.a/np.sqrt(1 - self.ecc2 * sin(phi)**2)
        sigma = self.sigma(phi)
        xgk = sigma +(dlam**2/2) * N * sin(phi)*cos(phi)*((1+(dlam**2/12)*(cos(phi))**2*(5-t**2+9*eta2+4*eta2**2)+(dlam**4/360)*cos(phi)**4*(61-58*t**2+t**4+270*eta2-330*eta2*t**2)))  
        ygk = dlam * N * cos(phi) * (1 + ((dlam**2)/6) * (cos(phi)**2) * (1 - t**2 + eta2) + (dlam**4/120) * (cos(phi)**4) * (5 - 18 * t**2 + t**4 + 14 * eta2 - 58 * eta2 * t**2))
                             
        x92 = xgk * m92 - 5300000
        y92 = ygk * m92 + 500000
        return x92, y92 



if __name__ == "__main__":
    
    print(sys.argv)
    
    if len(sys.argv) < 2:
        print("Musisz podac więcej argumentow")
        sys.exit(1)
        
    model = sys.argv[1]
    if model not in ['wgs84', 'grs80']:
        print("Mozesz podac tylko jedna elipsoide: 'wgs84' lub 'grs80'")
        sys.exit(1)
    else:
        if model == 'wgs84':
            print("Wybrano elipsoide WGS84")
        elif model == 'grs80':
            print("Wybrano elipsoide GRS80")
            
    geo = Transformacje(model)
    
    
    if '--header_lines':
        number_of_header_lines = int(sys.argv[4])
    else:
        number_of_header_lines = 0
        
    input_file_path = sys.argv[-1]
    
    if '--xyz2plh' in sys.argv and '--plh2xyz' in sys.argv[2]:
        print('mozesz podac tylko jedna flage!')
    
    elif '--xyz2plh' in sys.argv:
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[number_of_header_lines:]
        
            coords_plh = []
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                x_str,y_str,z_str = coord_line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                phi,lam,h = geo.xyz2plh(x, y, z)
                coords_plh.append([phi, lam, h])
                
        
        with open('result_xyz2plh.txt', 'w') as f:
            f.write('phi[deg], lam[deg], h[m]\n')
            for coords_list in coords_plh:
                line = ','.join([f'{coord:11.7f}' for coord in coords_list])
                f.writelines(line + '\n')
                
    elif '--plh2xyz' in sys.argv:
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[number_of_header_lines:]
            
            coords_xyz = []
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                phi_str, lam_str, h_str = coord_line.split(',')
                phi, lam, h = (float(phi_str), float(lam_str), float(h_str))
                x, y, z = geo.plh2xyz(phi, lam, h)
                coords_xyz.append([x, y, z])
                
        with open('result_plh2xyz.txt', 'w') as f:
            f.write('x[m], y[m], z[m]\n')
            for coords_list in coords_xyz:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')
                
    
    elif '--xyz2neu' in sys.argv:
        if len(sys.argv) < 6:
            print("Nalezy podac wspolrzedne x0, y0, z0.")
            sys.exit(1)
        else:
            x0 = float(input("Podaj wspolrzedna x0: "))
            y0 = float(input("Podaj wspolrzedna y0: "))
            z0 = float(input("Podaj wspolrzedna z0: "))
            
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[number_of_header_lines:]
              
            coords_neu = []
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                x, y, z = coord_line.split(',')
                x, y, z = float(x), float(y), float(z)
                n, e, u = geo.xyz2neu(x, y, z, x0, y0, z0)
                coords_neu.append([n, e, u])
                  
        with open('result_xyz2neu.txt', 'w') as f:
            f.write('easting[m], northing[m], up[m]\n')
            for coords_list in coords_neu:
                line = ','.join([f'{coord:11.3f}' for coord in coords_list])
                f.writelines(line + '\n')
                 
    elif '--bl2pl2000' in sys.argv:
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[number_of_header_lines:]
            
            coords_pl2000 = []
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                phi_str, lam_str = coord_line.split(',')
                phi, lam = (float(phi_str), float(lam_str))
                x2000, y2000 = geo.bl2pl2000(phi, lam)
                coords_pl2000.append([x2000, y2000])
                
        with open('result_bl2pl2000.txt', 'w') as f:
            f.write('x2000[m], y2000[m]\n')
            for coords_list in coords_pl2000:
                line = ','.join([f'{coord:11.3f}' for coord in coords_list])
                f.writelines(line + '\n')
                
    elif '--bl2pl1992' in sys.argv:
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[number_of_header_lines:]
            
            coords_pl1992 = []
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                phi_str, lam_str = coord_line.split(',')
                phi, lam = (float(phi_str), float(lam_str))
                x1992, y1992 = geo.bl2pl1992(phi, lam)
                coords_pl1992.append([x1992, y1992])
                
        with open('result_bl2pl1992.txt', 'w') as f:
            f.write('x1992[m], y1992[m]\n')
            for coords_list in coords_pl1992:
                line = ','.join([f'{coord:11.3f}' for coord in coords_list])
                f.writelines(line + '\n')
                