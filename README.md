# FazPPP

**GPS taşıyıcı faz ölçüleriyle istasyon konumu, tarayıcıda.**
*Station positioning from GPS carrier phase, in the browser.* (English below)

Gözlem, navigasyon, SP3 ve meteoroloji dosyalarını yüklediğinizde FazPPP şunları hesaplar:

- RINEX okuma ve ön işleme: donmuş kayıtlar, alıcı saati 1 ms sıçramaları, GF ve Melbourne–Wübbena ile faz kopmaları, yaylar
- Uydu konumu ve saati: SP3 Lagrange ara değerlemesi, gönderim anı, rölativistik düzeltme, Dünya dönüşü, yükseklik ve azimut
- Troposfer: Saastamoinen + Herring eşleme fonksiyonları
- Alıcı ve uydu anteni (ANTEX PCO/PCV), katı yer gelgiti (IERS 2010), Shapiro gecikmesi, faz sarılması (Wu vd. 1993)
- L3 (iyonosferden bağımsız) faz + P3 kod en küçük kareler çözümü: statik konum, seçilen epokta kinematik konum, alıcı saati, ıslak troposfer, ondalıklı belirsizlikler, DOP
- İyonosfer: fazın koda seviyelenmesi, ince tabaka modeli, alıcı ve uydu DCB'leri, STEC/VTEC
- Klobuchar modeliyle L1 ve L2 tek frekans faz çözümleri
- MP1/MP2 kod çok yolluluğu
- Yayın efemerisiyle P3 kod tek nokta konumu
- İsteğe bağlı referans koordinatla (ör. ITRF) bütün çözümlerin karşılaştırılması

Hesaplar Python (yalnızca numpy) ile yazılmıştır ve [Pyodide](https://pyodide.org) sayesinde tarayıcıda çalışır. Dosyalar hiçbir sunucuya gönderilmez.

## Kullanım

1. Siteyi açın, dosyaları ilgili kutulara sürükleyin.
2. Hesap epoğunu (GPS zamanı) ve isterseniz referans koordinatı girin.
3. **Hesapla**'ya basın. İlk açılışta Python ortamı 10-20 saniyede indirilir, sonraki hesaplar birkaç saniye sürer.

Denemek için **Örnek veriyle dene** düğmesi ISTA istasyonunun 28.08.2014 verisini yükler (`examples/ista_2014_240/`).

| Dosya | Gerekli | Not |
|---|---|---|
| Gözlem (.o / .rnx) | evet | RINEX 2.xx veya 3.xx, GPS; P1/C1, P2, L1, L2 |
| SP3 | evet | IGS kesin yörünge; birden çok dosya seçilebilir |
| Navigasyon (.n) | hayır | Kod çözümü ve Klobuchar (başlıkta ION ALPHA/BETA varsa) |
| Meteoroloji (.m) | hayır | Yoksa standart atmosfer |
| ANTEX (.atx) | hayır | Yoksa gömülü `igs08.atx` (2017 öncesi veriler için) |
| DCB | hayır | CODE P1-P2; yoksa uydu DCB'leri veriden kestirilir |

`.gz` dosyalar kabul edilir. Hatanaka (`.d`, `.crx`) ve `.Z` dosyalarını önce açın.

## Sınırlamalar

- Belirsizlikler ondalıklı bırakılır (tam sayı çözümü yok).
- Uydu saatleri SP3 aralığında bilindiği için faz çözümü yalnızca SP3 epoklarını kullanır (CLK desteği yok).
- Okyanus yüklemesi, iyonosferin yüksek dereceli terimleri ve troposfer gradyanları modellenmez.
- Yalnızca GPS.

ISTA örneğinde statik çözüm ITRF2014 koordinatından (İzmit depremi sonrası deformasyon dahil) 1,2 cm, CSRS-PPP sonucundan yaklaşık 1 cm farklıdır.

## Komut satırı

```bash
pip install numpy
python py/cli.py examples/ista_2014_240/ISTA2400.14o examples/ista_2014_240/igs18074.sp3 \
  --nav examples/ista_2014_240/ISTA2400.14n --met examples/ista_2014_240/ISTA2400.14m \
  --epoch 13:15:00 --out sonuc.json
```

## GitHub Pages'te yayınlama

1. Bu klasörün içeriğini bir GitHub reposunun köküne yükleyin.
2. Repo **Settings → Pages** bölümünde *Source* olarak `Deploy from a branch`, dal olarak `main` ve klasör olarak `/ (root)` seçin.
3. Birkaç dakika sonra site `https://<kullanıcı-adı>.github.io/<repo-adı>/` adresinde açılır.

Yerelde denemek için klasörde `python -m http.server` çalıştırıp `http://localhost:8000` adresini açın (dosyayı doğrudan çift tıklayarak açmak Web Worker nedeniyle çalışmaz).

## Klasör yapısı

```
index.html            arayüz
worker.js             Pyodide'i yükleyip hesabı ayrı iş parçacığında çalıştırır
assets/               stil, arayüz kodu, metinler (TR/EN)
py/gnssfaz/           hesap paketi (readers, preprocess, orbit, corrections, estimate, pipeline ...)
py/web_entry.py       tarayıcı giriş noktası
py/cli.py             komut satırı
data/igs08.atx.gz     gömülü anten kalibrasyonu
examples/             ISTA örnek verisi
```

## Kaynaklar

Hofmann-Wellenhof, Lichtenegger ve Wasle (2008) *GNSS*; Xu (2007) *GPS*; Kouba (2009) *A Guide to Using IGS Products*; Petit ve Luzum (2010) *IERS Conventions*; Wu vd. (1993); Schaer (1999); Montenbruck ve Gill (2000).

Örnek veri ve `igs08.atx`: International GNSS Service (IGS).

---

## English

FazPPP computes a GPS station position from carrier-phase observations entirely in the browser (Python + numpy on Pyodide). Upload RINEX observation, navigation, meteorological and SP3 files; it runs preprocessing (clock jumps, cycle slips, arcs), satellite orbits and clocks, troposphere, antenna/tide/Shapiro/wind-up corrections, an L3 + P3 least-squares solution (static and kinematic at a chosen epoch, receiver clock, wet troposphere, float ambiguities, DOP), ionosphere/TEC with DCB estimation, single-frequency L1/L2 solutions with Klobuchar, MP1/MP2 multipath and a broadcast-ephemeris code solution, and compares them with an optional reference coordinate.

To publish: push the folder to a GitHub repository and enable **Settings → Pages → Deploy from a branch → main / (root)**. To run locally: `python -m http.server` in the folder, then open `http://localhost:8000`.

Built by Salih Tatcı for the Satellite Geodesy course (Geomatics Engineering MSc, Istanbul Technical University). Results are for teaching purposes.

## Lisans / License

MIT, see `LICENSE`.
