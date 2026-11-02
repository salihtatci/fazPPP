# FazPPP

**GPS ve Galileo taşıyıcı faz ölçüleriyle statik ve kinematik konum, tarayıcıda.**
*Static and kinematic positioning from GPS and Galileo carrier phase, in the browser.* (English below)

Gözlem, SP3, saat (CLK), navigasyon ve meteoroloji dosyalarını yüklediğinizde FazPPP şunları hesaplar:

- RINEX okuma ve ön işleme (GPS L1/L2, Galileo E1/E5a): donmuş kayıtlar, alıcı saati 1 ms sıçramaları, GF ve Melbourne–Wübbena ile faz kopmaları (yüksek hızlı veride ek pencere testi), yaylar
- Uydu konumu ve saati: SP3 Lagrange ara değerlemesi, RINEX saat (CLK) dosyasından ya da SP3'ten uydu saati, gönderim anı, rölativistik düzeltme, Dünya dönüşü, yükseklik ve azimut
- Troposfer: Saastamoinen + Herring eşleme fonksiyonları
- Alıcı ve uydu anteni (ANTEX PCO/PCV), katı yer gelgiti (IERS 2010), Shapiro gecikmesi, faz sarılması (Wu vd. 1993)
- L3 (iyonosferden bağımsız) faz + P3 kod çözümü, en küçük kareler ya da **Kalman filtresi + RTS yumuşatıcı** ile: statik konum (ve seçilen epokta ayrı konum) ya da her epokta ayrı konumla **kinematik rota**, alıcı saati, ıslak troposfer, GPS–Galileo sistemler arası saat farkı (ISB), ondalıklı belirsizlikler, DOP. Epok parametreleri epok epok önceden elimine edilir; 1 s'lik veri de çözülebilir
- İyonosfer (GPS): fazın koda seviyelenmesi, ince tabaka modeli, alıcı ve uydu DCB'leri, STEC/VTEC
- Klobuchar modeliyle L1 ve L2 tek frekans faz çözümleri (GPS)
- MP1/MP2 kod çok yolluluğu (GPS)
- Yayın efemerisiyle P3 kod tek nokta konumu (GPS); statik modda alıcı hareket ediyorsa uyarı
- İsteğe bağlı referans koordinatla (ör. ITRF) bütün çözümlerin karşılaştırılması
- Bütün sonuçların (tablolar ve grafikler) tek PDF raporu olarak indirilmesi; ayrıca JSON ve tablo bazında CSV

Hesaplar Python (yalnızca numpy) ile yazılmıştır ve [Pyodide](https://pyodide.org) sayesinde tarayıcıda çalışır. Dosyalar hiçbir sunucuya gönderilmez.

## Kullanım

1. Siteyi açın, dosyaları ilgili kutulara sürükleyin.
2. Çözüm türünü (statik / kinematik) seçin, hesap epoğunu (GPS zamanı) ve isterseniz referans koordinatı girin.
3. **Hesapla**'ya basın. İlk açılışta Python ortamı 10-20 saniyede indirilir, sonraki hesaplar birkaç saniye sürer.

Denemek için **Örnek veriyle dene** düğmesi ISTA istasyonunun 28.08.2014 verisini yükler (`examples/ista_2014_240/`).

| Dosya | Gerekli | Not |
|---|---|---|
| Gözlem (.o / .rnx) | evet | RINEX 2.xx veya 3.xx; GPS (P1/C1, P2, L1, L2) ve isteğe bağlı Galileo (E1 + E5a ya da E5b) |
| SP3 | evet | IGS kesin yörünge; Galileo için çok sistemli (ör. `COD0MGXFIN`) ürün; birden çok dosya seçilebilir |
| Saat (.clk) | hayır | RINEX uydu saatleri (ör. 30 s). Kinematik ve yüksek hızlı veri için gerekli; SP3 ile aynı merkezden |
| Navigasyon (.n) | hayır | Kod çözümü ve Klobuchar (başlıkta ION ALPHA/BETA varsa) |
| Meteoroloji (.m) | hayır | Yoksa standart atmosfer |
| ANTEX (.atx) | hayır | Yoksa tarihe göre gömülü `igs20.atx` (27.11.2022 sonrası; GPS + Galileo, küçültülmüş) ya da `igs08.atx` |
| DCB | hayır | CODE P1-P2; yoksa uydu DCB'leri veriden kestirilir |

`.gz` dosyalar kabul edilir. Hatanaka (`.d`, `.crx`) ve `.Z` dosyalarını önce açın.

## Kalman filtresi

**Çözüm türü → Kinematik, Kalman filtresi** seçildiğinde ölçüler epok epok genişletilmiş Kalman filtresiyle (EKF) işlenir (`py/gnssfaz/estimate.py`, `kalman_solution`).

Durum vektörü: `x = [X, Y, Z, (Vx, Vy, Vz), c·dt, ZWD, ISB, B_1 … B_n]`

1. **Tahmin (zaman güncellemesi):** `x⁻ = Φ x⁺`, `P⁻ = Φ P⁺ Φᵀ + Q`
   - Sabit hız modeli: `Φ = [[I, Δt I], [0, I]]`, beyaz ivme gürültüsü `Q_XX = S Δt³/3`, `Q_XV = S Δt²/2`, `Q_VV = S Δt`, `S = diag(q_h, q_h, q_v)` (NEU'da, ECEF'e döndürülür)
   - Rastgele yürüyüş: `Q_XX = S Δt`; statik: `Q_XX = 0`
   - Alıcı saati her epokta yeniden başlatılır (beyaz gürültü), ZWD rastgele yürüyüş, ISB ve yay belirsizlikleri sabittir. Yay başladığında belirsizlik `L3 − P3` farkından büyük varyansla başlatılır.
2. **Güncelleme (ölçü güncellemesi):** L3 faz ve P3 kod, `v = z − h(x⁻)`, `S = H P⁻ Hᵀ + R`, `K = P⁻ Hᵀ S⁻¹`, `x⁺ = x⁻ + K v`, `P⁺ = (I − K H) P⁻ (I − K H)ᵀ + K R Kᵀ` (Joseph biçimi). Doğrusallaştırma yinelemeli EKF ile yapılır. `R = σ² / sin²E`.
3. **Yenilik testi:** `|v_i| / √S_ii > 5` olan kod ölçüsü atılır; faz ölçüsünde o yayın belirsizliği yeniden başlatılır (faz kopması).
4. **RTS yumuşatıcı (geri geçiş):** `C_k = P⁺_k Φᵀ (P⁻_{k+1})⁻¹`, `x^s_k = x⁺_k + C_k (x^s_{k+1} − x⁻_{k+1})`, `P^s_k = P⁺_k + C_k (P^s_{k+1} − P⁻_{k+1}) C_kᵀ`.

İleri filtre yalnızca o ana kadarki ölçüleri kullandığı için başta kod düzeyinde başlar ve belirsizlikler belirlendikçe yakınsar; "Kinematik rota" sekmesindeki standart sapma grafiği bu yakınsamayı gösterir. Yumuşatılmış çözüm toplu en küçük kareler çözümüne eşdeğerdir: ge.23O verisinde iki çözüm yatayda 1 mm'den, düşeyde (farklı troposfer modeli nedeniyle) 4 mm'den az farklıdır; ISTA'da statik modelli Kalman filtresi statik EKK çözümünü 0,5 mm içinde verir.

Komut satırı:

```bash
python py/cli.py gozlem.23O COD0MGXFIN_..._ORB.SP3.gz --clk COD0MGXFIN_..._CLK.CLK.gz \
  --mode kalman --kf-dyn pv --kf-qh 0.5 --kf-qv 0.01 --rate 1 --out kalman.json
```

## Sınırlamalar

- Belirsizlikler ondalıklı bırakılır (tam sayı çözümü yok).
- Saat (CLK) dosyası verilmezse uydu saatleri SP3 aralığında bilindiği için faz çözümü yalnızca SP3 epoklarını kullanır.
- Okyanus yüklemesi, iyonosferin yüksek dereceli terimleri ve troposfer gradyanları modellenmez.
- GPS ve Galileo kullanılır; GLONASS, BeiDou ve QZSS atlanır. İyonosfer, tek frekans, çok yolluluk ve kod çözümü analizleri yalnızca GPS ile yapılır.
- Kinematik çözümde düzeltmeler (yükseklik açısı, gelgit, anten) yaklaşık konumda hesaplanır; birkaç km'ye kadar hareket için etkisi ihmal edilebilir.

ISTA örneğinde statik çözüm ITRF2014 koordinatından (İzmit depremi sonrası deformasyon dahil) 1,2 cm, CSRS-PPP sonucundan yaklaşık 1 cm farklıdır.

## Komut satırı

```bash
pip install numpy
python py/cli.py examples/ista_2014_240/ISTA2400.14o examples/ista_2014_240/igs18074.sp3 \
  --nav examples/ista_2014_240/ISTA2400.14n --met examples/ista_2014_240/ISTA2400.14m \
  --epoch 13:15:00 --out sonuc.json

# hareketli alıcı: GPS + Galileo, 30 s saat dosyasıyla kinematik
python py/cli.py gozlem.23O COD0MGXFIN_20233070000_01D_05M_ORB.SP3.gz \
  --clk COD0MGXFIN_20233070000_01D_30S_CLK.CLK.gz --mode kinematic --out rota.json
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
data/igs08.atx.gz     gömülü anten kalibrasyonu (2022 öncesi)
data/igs20.atx.gz     gömülü anten kalibrasyonu (GPS + Galileo uyduları ve alıcılar, NOAZI)
examples/             ISTA örnek verisi
```

## Kaynaklar

Hofmann-Wellenhof, Lichtenegger ve Wasle (2008) *GNSS*; Xu (2007) *GPS*; Kouba (2009) *A Guide to Using IGS Products*; Petit ve Luzum (2010) *IERS Conventions*; Wu vd. (1993); Schaer (1999); Montenbruck ve Gill (2000).

Örnek veri, `igs08.atx` ve `igs20.atx`: International GNSS Service (IGS).

---

## English

FazPPP computes static or kinematic positions from GPS and Galileo carrier-phase observations entirely in the browser (Python + numpy on Pyodide). Upload RINEX observation, SP3, clock (CLK), navigation and meteorological files; it runs preprocessing (clock jumps, cycle slips, arcs), satellite orbits and clocks, troposphere, antenna/tide/Shapiro/wind-up corrections, an L3 + P3 solution by least squares or by an extended Kalman filter with RTS smoothing (static, or a kinematic track with a position at every epoch; receiver clock, wet troposphere, GPS–Galileo inter-system bias, float ambiguities, DOP), ionosphere/TEC with DCB estimation, single-frequency L1/L2 solutions with Klobuchar, MP1/MP2 multipath and a broadcast-ephemeris code solution, and compares them with an optional reference coordinate.

To publish: push the folder to a GitHub repository and enable **Settings → Pages → Deploy from a branch → main / (root)**. To run locally: `python -m http.server` in the folder, then open `http://localhost:8000`.

Built by Salih Tatcı for the Satellite Geodesy course (Geomatics Engineering MSc, Istanbul Technical University). Results are for teaching purposes.

## Lisans / License

MIT, see `LICENSE`.
