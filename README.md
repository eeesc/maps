# Nejmenší elektorát České republiky - Interaktivní mapa

Interaktivní mapa zobrazující obce s nejmenším počtem voličů v České republice v roce 2021, včetně demografických dat ze sčítání lidu.

## 🗺️ Živá mapa

Otevřete `index.html` ve vašem prohlížeči nebo navštivte: [GitHub Pages odkaz] (po nahrání na GitHub)

## 📊 Data

Mapa zobrazuje 10 obcí s nejmenším elektorátem v ČR, včetně:

- **Volební data (2021)**: Počet voličů a vliv jednoho hlasu
- **Demografická data**: Počet mužů, žen 15+, celkem obyvatel 15+
- **Geografické údaje**: Přesné souřadnice a odkazy na Mapy.cz
- **Administrativní údaje**: Kraj, správní centrum, kód obce

### Zobrazovaná data (původní sloupce E, G, J, L):
- **E**: Počet voličů v roce 2021
- **G**: Počet mužů
- **J**: Počet žen starších 15 let
- **L**: Celkový počet obyvatel starších 15 let

## 🎛️ Interaktivní funkce

### Ovládací prvky:
- **Zbarvit podle**: Vyberte metriku pro barevné rozlišení obcí
- **Velikost podle**: Vyberte metriku pro velikost markerů na mapě

### Dostupné metriky:
- Počet voličů (2021)
- Počet mužů
- Ženy 15+
- Celkem 15+

### Funkce mapy:
- **Kliknutím na marker**: Zobrazí detailní informace o obci
- **Zoom a pan**: Standardní ovládání mapy
- **Odkaz na Mapy.cz**: Přímý odkaz pro každou obec
- **Responzivní design**: Funguje na mobilních zařízeních

## 🛠️ Technické detaily

### Použité technologie:
- **Leaflet.js**: Interaktivní mapová knihovna
- **OpenStreetMap**: Mapové podklady
- **Vanilla JavaScript**: Žádné další závislosti
- **CSS3**: Moderní styling

### Struktura souborů:
```
mapy/
├── index.html          # Hlavní HTML soubor s mapou
├── data.js            # Sloučená data (souřadnice + census)
├── coordinates.js     # Původní souřadnice (záloha)
├── nejmenší_elektorát.xlsx  # Původní Excel data
├── merge_simple.py    # Script pro sloučení dat
└── README.md          # Tento soubor
```

## 🚀 Jak sdílet přes GitHub

### 1. Nahrání na GitHub:
```bash
git add .
git commit -m "Add interactive map of smallest electorates"
git push origin main
```

### 2. Aktivace GitHub Pages:
1. Jděte do Settings vašeho GitHub repository
2. Scroll dolů na "Pages"
3. V "Source" vyberte "Deploy from a branch"
4. Vyberte "main" branch
5. Klikněte "Save"

### 3. Sdílení:
Vaše mapa bude dostupná na: `https://[username].github.io/[repository-name]/`

## 📈 Možná rozšíření

- **Další demografická data**: Věkové skupiny, vzdělání, zaměstnanost
- **Historická data**: Porovnání s předchozími volbami
- **Export dat**: CSV, JSON formáty
- **Filtrování**: Podle kraje, velikosti obce
- **Statistiky**: Grafy a trendy
- **Mobilní aplikace**: Progressive Web App

## 📋 Zdroje dat

- **Volební data**: Český statistický úřad
- **Demografická data**: Sčítání lidu 2021
- **Geografická data**: Extrahováno z Mapy.cz
- **Mapové podklady**: OpenStreetMap

## 🤝 Přispívání

Návrhy na vylepšení jsou vítány! Otevřete issue nebo pošlete pull request.

## 📄 Licence

Tento projekt je k dispozici pod MIT licencí. Data jsou z veřejných zdrojů ČSÚ. 