import pandas as pd
import math
import numpy as np
import random
import matplotlib.pyplot as plt
import graphviz
import imageio

class Mutacja():
    def __init__(self, id, cecha, zmiana, wartosc, pokolenie):
        self.id = id
        self.cecha = cecha
        self.zmiana = zmiana
        self.wartosc = wartosc
        self.pokolenie = pokolenie
        self.kolejne_mutacje = []
        self.kolejne_mutacje = None

    def __str__(self) -> str:
        return f"Mutacja\nid: {self.id}\ncecha: {self.cecha}, pokolenie:{self.pokolenie}, zmiana: {self.zmiana}"

    def poprzednia_mutacja(self, poprzednia_mutacja):
        self.poprzednia_mutacja = poprzednia_mutacja

    def kolejne_mutacje(self, kolejna_mutacja):
        self.kolejne_mutacje.append(kolejna_mutacja)

class Osobnik:
    def __init__(self):
        self.korzystne_mutacje = {0:0, 1:0}

    def __str__(self) -> str:
        return f"Plec: {self.plec}, Genotyp: {self.genotyp}, Mutacje: {self.korzystne_mutacje}"

    def wylosuj_plec(self) -> None:
        p = random.random()
        if p < 0.5:
            self.plec = 0
        else:
            self.plec = 1

    def wylosuj_genotyp(self, liczba_cech: int, srednia: float, wariancja: float) -> None:
        self.genotyp = np.random.normal(loc=srednia, scale=np.sqrt(wariancja), size=liczba_cech)

class Populacja:
    def __init__(self,
                liczba_osobnikow: int,
                liczba_cech: int,
                mi_cechy: float,
                wariancja_cech: float):
        self.lista_osobnikow = self.przygotuj_osobnikow_na_poczatek_symulacji(liczba_osobnikow=liczba_osobnikow,
                                                                              liczba_cech=liczba_cech,
                                                                              srednia=mi_cechy,
                                                                              wariancja=wariancja_cech)
    def __str__(self) -> str:
        return f"Populacja o {len(self.lista_osobnikow)} osobnikach"

    def samce(self):
        return [osobnik for osobnik in self.lista_osobnikow if osobnik.plec == 1]

    def samice(self):
        return [osobnik for osobnik in self.lista_osobnikow if osobnik.plec == 0]

    def przygotuj_osobnikow_na_poczatek_symulacji(self, liczba_osobnikow: int, liczba_cech: int, srednia: float, wariancja: float):
        osobnicy = []
        for i in range(liczba_osobnikow):
            osobnik = Osobnik()
            osobnik.wylosuj_genotyp(liczba_cech, srednia, wariancja)
            if i >= liczba_osobnikow/2:
                osobnik.plec = 0
            else:
                osobnik.plec = 1
            osobnicy.append(osobnik)
        return osobnicy

    def zostaw_potomka(self, samiec: Osobnik, samica: Osobnik, mutacja: bool, srodowisko, nr_pokolenia):

        assert samiec.plec == 1 and samica.plec == 0, "Trzeba podac osobnikow plci meskiej(1) i zenskiej(0)"

        genotyp_ojca = samiec.genotyp
        genotyp_matki = samica.genotyp

        dziecko = Osobnik()
        dziecko.wylosuj_plec()
        genotyp_dziecka = np.zeros(len(genotyp_matki))

        # Dziedziczenie cech odbywa się losowo.
        # Każda cecha od ojca i matki jest wybierana z p = 0.5 i przekazywana dziecku. Cechy są dziedziczone niezależnie

        for i in range(len(genotyp_matki)):
            p = random.random()
            if p < 0.5:
                genotyp_dziecka[i] = genotyp_matki[i]
                dziecko.korzystne_mutacje[i] = samica.korzystne_mutacje[i]
            else:
                genotyp_dziecka[i] = genotyp_ojca[i]
                dziecko.korzystne_mutacje[i] = samiec.korzystne_mutacje[i]
            if dziecko.korzystne_mutacje[i] != 0:
                ostatnia_mutacja = dziecko.korzystne_mutacje[i]
                if ostatnia_mutacja.id in zliczanie_korzystnych_mutacji:
                    if ostatnia_mutacja.zmiana in zliczanie_korzystnych_mutacji[ostatnia_mutacja.id]:
                        if nr_pokolenia in zliczanie_korzystnych_mutacji[ostatnia_mutacja.id][ostatnia_mutacja.zmiana]:
                            zliczanie_korzystnych_mutacji[ostatnia_mutacja.id][ostatnia_mutacja.zmiana][nr_pokolenia] += 1
                        else:
                            zliczanie_korzystnych_mutacji[ostatnia_mutacja.id][ostatnia_mutacja.zmiana][nr_pokolenia] = 1
                    else:
                        zliczanie_korzystnych_mutacji[ostatnia_mutacja.id][ostatnia_mutacja.zmiana] = {nr_pokolenia : 1}
                else:
                    zliczanie_korzystnych_mutacji[ostatnia_mutacja.id] = {ostatnia_mutacja.zmiana : {nr_pokolenia : 1}}
        dziecko.genotyp = genotyp_dziecka
        return dziecko

    def wszystkie_osobniki(self):
        genotypy = []
        for osobnik in self.lista_osobnikow:
            genotypy.append(osobnik.genotyp)
        df_genotypy = pd.DataFrame(genotypy)
        return df_genotypy

class Srodowisko:
    def __init__(self, najlepszy_genotyp: list):
        self.najlepszy_genotyp = najlepszy_genotyp

    def __str__(self) -> str:
        return "Srodowisko o najlepszym genotypie: " + str(self.najlepszy_genotyp)

    def porownaj_osobnika_do_optimum(self, osobnik: Osobnik):
        genotyp_osobnika = osobnik.genotyp
        assert len(genotyp_osobnika) == len(self.najlepszy_genotyp), "Optymalny Genotyp musi być tej samej wielkości co genotyp osobnika"
        odlaglosc_kwadrat = sum([(xi - yi)**2 for xi, yi in zip(genotyp_osobnika, self.najlepszy_genotyp)])
        odleglosc = math.sqrt(odlaglosc_kwadrat)
        return odleglosc

    def prawdopodobienstwo_przezycia_osobnika(self, osobnik: Osobnik, sila_selekcji=1):
        fenotyp = self.porownaj_osobnika_do_optimum(osobnik)
        prawdopodbienstwo = np.exp(-fenotyp/(2*(sila_selekcji)**2))
        return prawdopodbienstwo

def krok_populacyjny_reprodukcja(populacja: Populacja, srodowisko: Srodowisko, nr_pokolenia, srednia_liczba_potomstwa: int = 2.1):
    osobniki_po_reprodukcji = []
    samce = populacja.samce()
    samice = populacja.samice()
    for samiec in samce:
        if len(samice) > 0:
            indeks_samicy = random.randint(0, len(samice)-1)
            samica_dla_samca = samice[indeks_samicy]
            liczba_potomstwa = np.random.poisson(lam=srednia_liczba_potomstwa, size=1).item()
            for i in range(liczba_potomstwa):
                p = nr_pokolenia
                potomek = populacja.zostaw_potomka(samiec, samica_dla_samca, mutacja=False, srodowisko=srodowisko, nr_pokolenia=p)
                osobniki_po_reprodukcji.append(potomek)
            samice.pop(indeks_samicy)
        else:
            break

    populacja.lista_osobnikow = osobniki_po_reprodukcji

def krok_populacyjny_selekcja(populacja: Populacja, srodowisko: Srodowisko, maksymalna_liczebnosc: int):
    # Etap I
    osobniki_przed_selekcja = populacja.lista_osobnikow
    osobniki_po_selekcji = []
    for osobnik in osobniki_przed_selekcja:
        czy_przezyje = (srodowisko.prawdopodobienstwo_przezycia_osobnika(osobnik)> random.random())
        if czy_przezyje:
            osobniki_po_selekcji.append(osobnik)
    #Etap II
    if len(osobniki_po_selekcji) > maksymalna_liczebnosc:
        liczba_osobnikow_do_usuniecia = len(osobniki_po_selekcji) - maksymalna_liczebnosc
        for i in range(liczba_osobnikow_do_usuniecia):
            losowy_indeks = random.randint(0, len(osobniki_po_selekcji)-1)
            osobniki_po_selekcji.pop(losowy_indeks)
        assert len(osobniki_po_selekcji) == maksymalna_liczebnosc

    populacja.lista_osobnikow = osobniki_po_selekcji

def wylosuj_mutowana_ceche(liczba_cech: int):
    return random.randint(0, liczba_cech-1)

def krok_populacyjny_mutacje(populacja: Populacja, srodowisko: Srodowisko, prawdop_losowej_mutacji: float, efekt_mutacji: float, nr_pokolenia: int):
    osobniki = populacja.lista_osobnikow
    for osobnik in osobniki:
        czy_osobnik_mutuje = (random.random() < prawdop_losowej_mutacji)
        if czy_osobnik_mutuje:
            numer_cechy = wylosuj_mutowana_ceche(len(osobnik.genotyp))
            dostosowane_bez_mutacji = srodowisko.porownaj_osobnika_do_optimum(osobnik)
            zmiana_cechy = np.random.normal(loc=0, scale=math.sqrt(efekt_mutacji), size=1)[0]
            stary_genotyp = osobnik.genotyp[numer_cechy]
            osobnik.genotyp[numer_cechy] += zmiana_cechy
            dostosowane_po_mutacji = srodowisko.porownaj_osobnika_do_optimum(osobnik)
            if dostosowane_po_mutacji > dostosowane_bez_mutacji:
                if osobnik.korzystne_mutacje[numer_cechy]==0:
                    id = (numer_cechy, stary_genotyp, zmiana_cechy, nr_pokolenia)
                    mutacja = Mutacja(id, numer_cechy, zmiana_cechy, osobnik.genotyp[numer_cechy], nr_pokolenia)
                    mutacja.poprzednia_mutacja(stary_genotyp)
                    sledzenie_korzystnych_mutacji[id] = [mutacja]
                else:
                    id = osobnik.korzystne_mutacje[numer_cechy].id
                    mutacja = Mutacja(id, numer_cechy, zmiana_cechy, osobnik.genotyp[numer_cechy], nr_pokolenia)
                    mutacja.poprzednia_mutacja(osobnik.korzystne_mutacje[numer_cechy].wartosc)
                    sledzenie_korzystnych_mutacji[id].append(mutacja)
                osobnik.korzystne_mutacje[numer_cechy] = mutacja
            else:
                osobnik.korzystne_mutacje[numer_cechy] = 0



def uderzenie_meteorytu(srodowisko, pokolenie):
    if pokolenie % 50 == 0:
        optymalny_genotyp = srodowisko.najlepszy_genotyp
        for i in range(len(optymalny_genotyp)):
            if random.random() < 0.5:
                optymalny_genotyp[i] = srodowisko.najlepszy_genotyp[i] - 1
            else:
                optymalny_genotyp[i] = srodowisko.najlepszy_genotyp[i] + 1
        srodowisko.najlepszy_genotyp = optymalny_genotyp

# Zakładam że ocieplenie klimatu to stała zmiana optymalnego genotypu w konkretnym kierunku.
# Ten kierunek oczywiście można wybrać.
# Ja przyjmuje że ocieplenie klimatu będzie powodować zawsze wzrost wszystkich cech optymalnego genotypu
def ocieplenie_klimatu(srodowisko: Srodowisko):
    optymalny_genotyp = srodowisko.najlepszy_genotyp
    for i in range(len(optymalny_genotyp)):
        zmiana = np.random.lognormal(sigma=0.01)/100
        optymalny_genotyp[i] += zmiana
    srodowisko.najlepszy_genotyp = optymalny_genotyp

def losowe_zmiany(srodowisko):
    #TODO
    pass

def stale_optimum(srodowisko):
    pass

def krok_populacyjny_zmiana_srodowiska(srodowisko: Srodowisko,
                                       pokolenie: int,
                                       losowo: bool=False,
                                       ocieplenie: bool = False,
                                       meteoryt: bool=False,
                                       stale: bool=False):
    if losowo:
        losowe_zmiany(srodowisko)
    elif stale:
        stale_optimum(srodowisko)
    elif ocieplenie:
        ocieplenie_klimatu(srodowisko)
    elif meteoryt:
        uderzenie_meteorytu(srodowisko, pokolenie)

def zapisz_polozenie_populacji(genotypy: pd.DataFrame, srodowisko: Srodowisko, pokolenie: int):
    x = []
    y = []

    if len(srodowisko.najlepszy_genotyp) == 2:
        x = genotypy.iloc[:, 0]
        y = genotypy.iloc[:, 1]
        najlepszy_genotyp = srodowisko.najlepszy_genotyp
    else: #Nie działa póki co dobrze
        # X = []
        # pca = PCA(n_components=2)
        # if  len(populacja.lista_osobnikow) > 0:
        #     for os in populacja.lista_osobnikow:
        #         X.append(os.genotyp)
        #     X.append(srodowisko.najlepszy_genotyp)
        #     X = pca.fit_transform(X)
        #     x = [i[0] for i in X[:-1]]
        #     y = [i[1] for i in X[:-1]]
        #     najlepszy_genotyp = X[-1]
        pass

    plt.scatter(x, y, marker='.', s=0.7)
    plt.scatter(najlepszy_genotyp[0], najlepszy_genotyp[1], color='red')
    plt.title(f'Epoch: {pokolenie}')
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    plt.savefig(f'frames\epoch-{pokolenie}.png')
    plt.close()

# Nasze parametry

liczba_pokolen = 100

# początkowa liczbeność populacji jak i pojemność środowiska
N = 1000

#liczba cech osobnika
n = 2 # dla n > niz dwa animacja nie jest za dobrze zrobiona :(

#średnia wartość cechy losowanej przy tworzeniu osobnika
mi_cechy = 0

#wariancja wartości cechy losowanej przy tworzeniu osobnika
var = 1

#prawdopodobieństwo wystąpienia losowej mutacji u osobnika
mutation_rate = 0.7

#siła mutacji
wariancja_mutacji = 0.1

# Parametry dla wizualizacji szlaków korzystnych mutacji

#minimalna liczba pokoleń podczad których utrzymała się dana mutacja
min_pokolen = 30

# Wybór optymalneg ogenomu, albo losowo, albo [0...0]
najlepszy_genotyp = np.zeros(n)
# np.random.normal(loc=mi_cechy, scale=np.sqrt(var), size=n)
assert len(najlepszy_genotyp) == n

populacja = Populacja(liczba_osobnikow=N,
                      liczba_cech=n,
                      mi_cechy=mi_cechy,
                      wariancja_cech=var)

srodowisko = Srodowisko(najlepszy_genotyp=najlepszy_genotyp)

srednia_cechy = [[] for i in range(n)]
wariancja_cechy = [[] for i in range(n)]
naj_genotyp_cech = [[] for i in range(n)]

global korzystne_mutacje_w_pokoleniach
global zmutowane
global zliczanie_korzystnych_mutacji
global sledzenie_korzystnych_mutacji
korzystne_mutacje_w_pokoleniach = {i:{} for i in range(n)}
zmutowane = {}
zliczanie_korzystnych_mutacji = dict()
sledzenie_korzystnych_mutacji = {}

# Ogólna pętla powinna wyglądać tak:
numery_pokolen = np.arange(liczba_pokolen)
liczba_osobnikow_per_pokolenie = []

for pokolenie in range(liczba_pokolen):

    krok_populacyjny_mutacje(populacja,
                             srodowisko,
                             prawdop_losowej_mutacji=mutation_rate,
                             efekt_mutacji=wariancja_mutacji,
                             nr_pokolenia = pokolenie)

    krok_populacyjny_selekcja(populacja, srodowisko, maksymalna_liczebnosc=N)

    krok_populacyjny_reprodukcja(populacja, srodowisko, srednia_liczba_potomstwa=3, nr_pokolenia=pokolenie)

    krok_populacyjny_zmiana_srodowiska(srodowisko, pokolenie, ocieplenie=True)

    krok_populacyjny_wszystkie_genotypy = populacja.wszystkie_osobniki()

    for i in range(n):
        naj_genotyp_cech[i].append(srodowisko.najlepszy_genotyp[i])
        if len(krok_populacyjny_wszystkie_genotypy[i]) > 0:
            srednia_cechy[i].append(krok_populacyjny_wszystkie_genotypy[i].mean())
            wariancja_cechy[i].append(krok_populacyjny_wszystkie_genotypy[i].var())
        else:
            srednia_cechy[i].append(0)
            wariancja_cechy[i].append(0)

    liczba_osobnikow_per_pokolenie.append(len(populacja.lista_osobnikow))
    zapisz_polozenie_populacji(krok_populacyjny_wszystkie_genotypy, srodowisko, pokolenie)

with imageio.get_writer('frames\line.gif', mode='i') as writer:
    for i in range(0, liczba_pokolen):
        image = imageio.imread(f'frames\epoch-{i}.png')
        writer.append_data(image)

plt.plot(numery_pokolen, liczba_osobnikow_per_pokolenie, color = 'green')
plt.title('liczbność populacji w kolejnych pokoleniach')
plt.xticks(range(0, 110, 10))
plt.xlabel('liczba pokoleń')
plt.ylabel('ilość osobników w populacji')
plt.savefig(f'liczebność_populacji_N:{N}_p:{liczba_pokolen}')

for nr_cechy in range(n):
    plt.plot(numery_pokolen, srednia_cechy[nr_cechy], color = 'red', label = 'średnia')
    plt.fill_between(numery_pokolen, np.array(srednia_cechy[nr_cechy])-np.array(wariancja_cechy[nr_cechy]),
                    np.array(srednia_cechy[nr_cechy])+np.array(wariancja_cechy[nr_cechy]), color = 'lightblue',
                     label = 'wariancja')
    plt.plot(numery_pokolen, naj_genotyp_cech[nr_cechy], color = 'navy', label = 'najlepszy genotyp')
    plt.title(f'cecha {nr_cechy}')
    plt.legend()
    plt.xticks(range(0, 110, 10))
    plt.xlabel('liczba pokoleń')
    plt.ylabel('wartość cechy')
    plt.savefig(f'ocieplenie_cecha_{nr_cechy}')
    plt.show()

kolory = ['blue', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflowerblue', 'cornsilk', 'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgreen', 'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 'darkslateblue', 'darkslategray', 'darkslategrey', 'darkturquoise', 'darkviolet', 'deeppink', 'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 'firebrick']
x = 0
for id_mutacji in sledzenie_korzystnych_mutacji:
    ps = graphviz.Digraph(f'mutacje:{x}', directory='korzystne_mutacje')
    zliczanie_wykres = dict()
    if id_mutacji in zliczanie_korzystnych_mutacji:
        pierwsze_pokolenie = 0
        ostatnie_pokolenie = 0
        for m in range(len(sledzenie_korzystnych_mutacji[id_mutacji])):
            mutacja = sledzenie_korzystnych_mutacji[id_mutacji][m]
            if mutacja.zmiana in zliczanie_korzystnych_mutacji[id_mutacji]:
                ps.node(str(np.round_(mutacja.poprzednia_mutacja, decimals=4)), fontsize="10")
                ps.node(str(np.round_(mutacja.wartosc, decimals=4)), fontsize="10")
                ps.edge(str(np.round_(mutacja.poprzednia_mutacja, decimals=4)),
                        str(np.round_(mutacja.wartosc, decimals=4)), label=str(mutacja.pokolenie),
                        headlabel=str(np.round_(mutacja.zmiana, decimals=4)), fontsize="10")
                p = mutacja.pokolenie
                pokolenia = set(zliczanie_korzystnych_mutacji[id_mutacji][mutacja.zmiana].keys())
                if min(pokolenia) < pierwsze_pokolenie:
                    pierwsze_pokolenie = min(pokolenia)
                if max(pokolenia) > ostatnie_pokolenie:
                    ostatnie_pokolenie = max(pokolenia)
                while p in pokolenia:
                    ile_osobnikow = zliczanie_korzystnych_mutacji[id_mutacji][mutacja.zmiana][p]
                    if id_mutacji in zliczanie_wykres:
                        if mutacja.zmiana in zliczanie_wykres[id_mutacji]:
                            zliczanie_wykres[id_mutacji][mutacja.zmiana][p] = ile_osobnikow
                        else:
                            zliczanie_wykres[id_mutacji][mutacja.zmiana] = np.zeros(liczba_pokolen)
                            zliczanie_wykres[id_mutacji][mutacja.zmiana][p] = ile_osobnikow
                    else:
                        zliczanie_wykres[id_mutacji] = {mutacja.zmiana : np.zeros(liczba_pokolen)}
                        zliczanie_wykres[id_mutacji][mutacja.zmiana][p] = ile_osobnikow
                    pokolenia.remove(p)
                    p += 1
        if ostatnie_pokolenie - pierwsze_pokolenie > min_pokolen:
            ps.save()
            k = 0
            if len(zliczanie_wykres[id_mutacji]) > 0:
                dane_do_pliku = open(f'mutacje:{x}', 'w')
                dane_do_pliku.write(f'cecha:{id_mutacji[0]} wartosc_pocz: {id_mutacji[1]} pokolenie_pocz: {id_mutacji[3]}\n')
                for zmiana in zliczanie_wykres[id_mutacji]:
                    alpha = abs(np.round_(zmiana, decimals=2))
                    if alpha > 1:
                        alpha = 1
                    if alpha > 0.1:
                        plt.fill_between(numery_pokolen, 0, zliczanie_wykres[id_mutacji][zmiana], color = kolory[k], alpha=alpha)
                        dane_do_pliku.write(f'zmiana: {zmiana}\n{list(zliczanie_wykres[id_mutacji][zmiana])}\n')
                        k += 1
                dane_do_pliku.close()
                plt.title(f'cecha:{id_mutacji[0]} wartosc_pocz: {id_mutacji[1]} pokolenie_pocz: {id_mutacji[3]}')
                plt.xticks(range(0, 110, 10))
                plt.xlabel('liczba pokoleń')
                plt.ylabel('liczba osobników')
                plt.savefig(f'mutacje:{x}')
                plt.close()
    x += 1

