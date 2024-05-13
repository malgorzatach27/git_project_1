## FUNKCJONALNOŚĆ:

Program ten został napisany w celu transformacji współrzędnych geograficznych wedle potrzeb użytkownika. Oferuje on transformacje:

1) wsp. XYZ do wsp. BLH (--xyz2plh)*
2) wsp. BLH do wsp. XYZ (--plh2xyz)
3) wsp. XYZ do układu NEU (--xyz2neu)
4) wsp. BL do układu PL2000 (--bl2pl2000)
5) wsp. BL do układu PL1992 (--bl2pl1992)

Obsługuje on takie elipsoidy jak:

- GRS80 (grs80)
- WGS84 (wgs84)
- Elipsoida Krasowskiego (krasowski)

*W nawiasach znajdują się nazwy argumentów, które reprezentują dane transformacje. Należy ich użyć podczas uruchamiania programu w Wierszu poleceń.

## WYMAGANIA PROGRAMU:

Aby program zadziałał poprawnie, urządzenie użytkownika powinno zawierać:

- System operacyjny Windows 10, 11
- Python w wersji 3.12.3
- Zaimportowane biblioteki numpy, math, sys

## OPIS DZIAŁANIA PROGRAMU:

Aby poprawnie uruchomić program należy podać argumenty reprezentujące wymagane dane odpowiednio w kolejności:

- model wybranej elipsoidy (grs80/wgs84/krasowski), na której chcemy dokonać transformacji
- rodzaj transformacji (np. --xyz2plh), którą chcemy wykonać
- ilość nagłówków do pominięcia (--header_lines [podaj liczbę]), jeżeli w pliku z danymi wejściowymi znajdują się nagłówki niezawierające interesujących nas współrzędnych
- nazwa pliku z jego rozszerzeniem (lub ścieżka), zawierający współrzędne wejściowe (np. wspolrzedne.txt)

Przed w.w argumentami należy dopisać również "-python geo_v1.py" - w przeciwnym wypadku program nie zadziała. 

W przypadku gdy użytkownik wpisze jedynie polecenie "-python geo_v1.py" wyświetli mu się wzór komendy oraz kilka przykładowych komend zawierających wszystkie niezbędne argumenty sys.argv. 
Użytkownik musi jedynie zmienić dane na swoje potrzeby.
 
## PRZYKŁADOWE WYWOŁANIE PROGRAMU ZA POMOCĄ WIERSZA POLECEŃ WINDOWS:

W celu poprawnego uruchomienia programu należy kolejno:

1) Otworzyć Wiersz poleceń Windows
2) Skopiować ścieżkę do folderu z plikiem ze współrzędnymi, które chcemy transformować
3) Przejść do tego folderu za pomocą komendy "cd [ścieżka do folderu]" w Wierszu poleceń
4) Uruchomić program poprzez wpisanie odpowiednich argumentów np.

python geo_v1.py grs80 --xyz2plh --header_lines 4 wspolrzedne.txt

Jeżeli wszystkie dane i argumenty zostały podane prawidłowo to w folderze z plikiem wejściowym utworzy się nowy plik tekstowy o nazwie "result_[rodzaj transformacji]".
W pliku tym znajdziemy przetransformowane przez program współrzędne.

Należy również pamiętać o poprawnym zapisie danych w pliku wejściowym:

- XYZ -> PLH dla danych z pliku (kolejno X[m], Y[m], Z[m]) otrzymujemy (kolejno phi[st.d]*, lamda[st.d], H[m])
- PLH -> XYZ analogicznie z (phi[st.d],lamda[st.d], H[m]) otrzymujemy (X[m], Y[m], Z[m])
- XYZ -> NEU  dla danych z pliku (kolejno X[m], Y[m], Z[m]) oraz danych wymaganych do wpisania w wierszu poleceń (X0[m], Y0[m], Z0[m]) 
otrzymujemy (E[m], N[m], U[m]) czyli (easting[m], northing[m], up[m])
- PL -> PL2000 dla danych z pliku (phi[st.d],lamda[st.d]) otrzymujemy (X2000[m], Y2000[m])
- PL -> PL1992 dla danych z pliku (phi[st.d],lamda[st.d]) otrzymujemy (X1992[m], Y1992[m])

Wszystkie dane należy oddzielić przecinkiem.

* st.d - stopnie dziesiętne

## PRZYKŁADOWE TRANSFORMACJE WYKONANE PRZEZ PROGRAM:

1) Transformacja XYZ -> PLH (grs80)

```
Dane wejściowe (X[m], Y[m], Z[m]):
3664840.500,1409154.690,5009572.270
3664840.510,1409154.680,5009572.267
3664840.520,1409154.670,5009572.267
```
```
Dane wyjściowe (phi[st.d]*, lamda[st.d], H[m]):
52.0979374, 21.0320720, 85.1698187
52.0979373, 21.0320718, 85.1709807
52.0979373, 21.0320716, 85.1745098
```

2) Transformacja XYZ -> NEU
```
Dane wejściowe (X[m], Y[m], Z[m]):
3664940.500,1409154.690,5009571.170
3664940.510,1409153.580,5009571.167
3664940.520,1409153.570,5009571.167
```
ORAZ 
```
Dane wymagane do wpisania w wierszu poleceń (X0[m], Y0[m], Z0[m]):
50.000
60.000
70.000
```
```
Dane wyjściowe (easting[m], northing[m], up[m]):
-1913367.136,5015125.842,3420468.986
-1913367.150,5015125.839,3420468.985
-1913367.164,5015125.839,3420468.983
```

3) Transformacja PL -> PL2000
```
Dane wejściowe (phi[st.d]*, lamda[st.d]):
51.0879374, 20.0420720
51.0879373, 20.0420718
51.0879373, 20.0420716
```
```
Dane wyjściowe (X2000[m], Y2000[m]):
5661869.269,7432888.693
5661869.258,7432888.679
5661869.258,7432888.665
```

## ZNANE BŁĘDY:

Program może nie zostać poprawnie uruchomiony jeśli:
 
- Współrzędne w pliku wejściowym zostały oddzielone za pomocą czegoś innego niż przecinek
- Została podana błędna ilość nagłówków do pominięcia z użyciem argumentu --header_lines
- Współrzędne zostały podane w nieodpowiedniej kolejności
  
