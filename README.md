# FazPPP

**GPS ve Galileo taşıyıcı faz ölçüleriyle statik ve kinematik konum, tarayıcıda.**
*Static and kinematic positioning from GPS and Galileo carrier phase, in the browser.* (English below)

Gözlem, SP3, saat (CLK), navigasyon ve meteoroloji dosyalarını yüklediğinizde FazPPP şunları hesaplar:

- RINEX okuma (Hatanaka `.crx`/`.d` ve `.gz`/`.Z` doğrudan) ve ön işleme (GPS L1/L2, Galileo E1/E5a): donmuş kayıtlar, alıcı saati 1 ms sıçramaları, alıcının LLI bayrakları, GF ve Melbourne–Wübbena ile faz kopmaları (yüksek hızlı veride ek pencere testi), yaylar
- **Kod ölçülerine Kalman filtresi:** her uydu ve yay için P1/P2'nin faz ile süzülmesi (iyonosfer ayrışmasız), yenilik testi
- Uydu konumu ve saati: SP3 Lagrange ara değerlemesi, RINEX saat (CLK) dosyasından ya da SP3'ten uydu saati, gönderim anı, rölativistik düzeltme, Dünya dönüşü, yükseklik ve azimut
- Troposfer: Saastamoinen + Herring eşleme fonksiyonları
- Alıcı ve uydu anteni (ANTEX PCO/PCV), katı yer gelgiti (IERS 2010), Shapiro gecikmesi, faz sarılması (Wu vd. 1993)
- L3 (iyonosferden bağımsız) faz + P3 kod en küçük kareler çözümü: statik konum (ve seçilen epokta ayrı konum) ya da her epokta ayrı konumla **kinematik rota**, alıcı saati, ıslak troposfer, GPS–Galileo sistemler arası saat farkı (ISB), ondalıklı belirsizlikler, DOP. Epok parametreleri epok epok önceden elimine edilir; 1 s'lik veri de çözülebilir
- İyonosfer (GPS): fazın koda seviyelenmesi, ince tabaka modeli, alıcı ve uydu DCB'leri, STEC/VTEC
- Klobuchar modeliyle L1 ve L2 tek frekans faz çözümleri (GPS)
- MP1/MP2 kod çok yolluluğu (GPS)
- Yayın efemerisiyle P3 kod tek nokta konumu (GPS); statik modda alıcı hareket ediyorsa uyarı
- İsteğe bağlı referans koordinatla (ör. ITRF) bütün çözümlerin karşılaştırılması
- Bütün sonuçların (tablolar ve grafikler) tek PDF raporu olarak indirilmesi; ayrıca JSON ve tablo bazında CSV

Hesaplar Python (yalnızca numpy) ile yazılmıştır ve [Pyodide](https://pyodide.org) sayesinde tarayıcıda çalışır. Dosyalar hiçbir sunucuya gönderilmez.

## Kullanım

1. Siteyi açın, dosyaları ilgili kutulara sürükleyin.
2. Çözüm türünü seçin (varsayılan **Otomatik**: alıcı hareket ediyorsa kinematik, etmiyorsa statik), isterseniz hesap epoğunu (GPS zamanı) ve referans koordinatı girin.
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

`.gz` ve `.Z` sıkıştırılmış ve Hatanaka (`.crx`, `.d`) dosyalar doğrudan yüklenebilir.

## Sağlamlık ve doğruluk

Her türlü gözlem dosyasında düşük hatayla çalışması için:

- **Otomatik çözüm türü:** Kod çözümünde alıcı 10 m'den fazla yer değiştiriyorsa kinematik seçilir; değiştirmiyorsa statik çözülür, statik çözüm ölçülere uymuyorsa (s₀ > 2,5) kinematik denenir ve belirgin biçimde daha iyiyse seçilir.
- **Artığa dayalı ayıklama:** Çözümden sonra faz artıklarında ~2 dk'lık pencerelerle basamak aranır; ön işlemenin kaçırdığı faz kopmalarında (ör. GF'yi 3 mm, MW'yi 2 devir değiştiren 9/7 devirlik kopma) yay bölünür, tek tek büyük faz ve kod artıkları çıkarılır, çözüm yinelenir (en çok 5 tur; her turda yalnız en belirgin basamaklar).
- **Alıcının LLI bayrakları** faz kopması olarak kullanılır (neredeyse her epokta işaretleyen alıcılarda yok sayılır).
- **Kaydırılmış ya da ondalıklı epoklar:** Gözlem epokları işlem aralığına oturmuyorsa her işlem anına en yakın epok seçilir.
- **Kısa dosyalar:** En kısa yay, dosya süresinin %30'una (2–10 dk) indirilir.
- **Ürün denetimi:** CLK ve SP3 uydu saatleri karşılaştırılır (farklı analiz merkezi uyarısı), gözlemlerin SP3/CLK kapsamı dışında kalan kısmı bildirilir; başlıktaki yaklaşık konum kod çözümünden 1 km'den fazla farklıysa kod çözümü kullanılır.

Test sonuçları (3 Kasım 2023, CODE MGEX final yörünge + 30 s saat, GPS + Galileo; doğruluk ölçütü IGS günlük koordinatları, `IGS0OPSSNX_20233070000_01D_01D_CRD`):

| Veri | Çözüm | Hata N / E / U [mm] |
|---|---|---|
| ISTA, 24 sa, 30 s, RINEX 3 | statik | 4 / −2 / −6 |
| BUCU, 24 sa, 30 s | statik | 3 / 0 / −8 |
| MATE, 24 sa, 30 s | statik | 3 / −4 / −15 |
| ISTA, CLK olmadan (5 dk SP3 saati) | statik | 4 / −2 / −4 |
| ISTA, 1,5 sa, 1 s | statik | 0 / −4 / 11 |
| ISTA, epoklar :15/:45 | statik | 0 / 3 / 16 |
| ISTA, yapay 3 faz kopması + kod aykırıları | statik | 4 / −2 / −6 (ayıklamasız 5 / −45 / 62) |
| ISTA, 24 sa, 30 s | kinematik, epok RMS | 7 / 5 / 23 |
| ISTA, 1,5 sa, 1 s | kinematik, epok RMS | 5 / 12 / 19 |
| ISTA 2014, RINEX 2, yalnız GPS | statik (ITRF2014'e göre) | −8 / 4 / 9 |

Kısa oturumlarda ondalıklı belirsizlikler tam yakınsamaz: 20 dk statik ~8 cm, 8 dk ~20 cm. Hareketli alıcı ge.23O (1 s, 68 dk) otomatik modda kinematik çözülür (faz RMS 8 mm).

## Kod ölçülerine Kalman filtresi

Kod ölçüleri (P1, P2) çözüme girmeden önce her uydu ve her faz yayı için ayrı ayrı, tek durumlu bir Kalman filtresiyle süzülür (`py/gnssfaz/codefilter.py`). Durum, o frekanstaki kod uzaklığıdır; kodun iki epok arasındaki değişimi fazdan bilindiği için filtre fazın hassasiyetini koda aktarır.

Ölçü modeli (frekans *i*, γ₁ = 1, γ₂ = (f₁/f₂)²):

```
P_i = ρ + γ_i I + ε_P          Φ_i = ρ − γ_i I + B_i + ε_Φ
ΔI  = Δ(Φ₁ − Φ₂) / (γ₂ − 1)
```

1. **Tahmin:** `x⁻_k = x⁺_(k−1) + ΔΦ_i + 2 γ_i ΔI`, `P⁻_k = P⁺_(k−1) + q`. İki frekansın fazıyla iyonosfer değişimi hesaba katıldığı için klasik Hatch filtresindeki kod–faz ayrışması (divergence) oluşmaz.
2. **Güncelleme:** `v = P_i − x⁻_k`, `K = P⁻_k / (P⁻_k + R)`, `x⁺_k = x⁻_k + K v`, `P⁺_k = (1 − K) P⁻_k`.
3. **Süreç gürültüsü:** `q = R κ² / (1 − κ)`, `κ = Δt / τ`; böylece kararlı durumda kazanç `K = Δt / τ` olur. τ (varsayılan 100 s) büyüdükçe süzme güçlenir; τ → ∞ kümülatif ortalamadır (Hatch filtresi).
4. **Yenilik testi:** `|v| > 5 √(P⁻ + R)` ise kod aykırı sayılır, güncelleme atlanır. Üst üste 3 ret olursa ya da yeni bir faz yayı başlarsa (faz kopması) filtre ham kodla yeniden başlar: `x = P_i`, `P = R`.

Süzülmüş kodlar statik/kinematik faz çözümüne, kod çözümüne, tek frekans çözümlerine ve iyonosfere girer; çok yolluluk (MP1/MP2) ham kodla hesaplanır. "Kalman filtresi" sekmesi uydu başına gürültü azalmasını (iyonosferden arındırılmış kod − faz farkı), kazancın ve σ'nın zamanla değişimini, kod çözümünün ve faz çözümünün ham ve süzülmüş kodla karşılaştırmasını gösterir.

Örnek sonuçlar (τ = 100 s): ISTA (30 s) P1 gürültüsü 0,58 → 0,27 m, P2 0,61 → 0,27 m; kod çözümünün faz çözümünden yatay farkı 1,6 → 1,1 m. ge.23O (1 s, hareketli alıcı) P1 0,23 → 0,07 m, P2 0,27 → 0,10 m. Faz çözümü konumu fazdan belirlendiği için 1 mm'den az değişir; kod artıkları ve s₀ küçülür.

Komut satırında filtreyi kapatmak için `--no-code-kf`, τ'yı değiştirmek için `--code-kf-tau 300`.

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
py/gnssfaz/           hesap paketi (readers, crinex, preprocess, codefilter, orbit, corrections, estimate, pipeline ...)
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

FazPPP computes static or kinematic positions from GPS and Galileo carrier-phase observations entirely in the browser (Python + numpy on Pyodide). Upload RINEX observation, SP3, clock (CLK), navigation and meteorological files; it runs preprocessing (clock jumps, cycle slips, arcs), a per-satellite Kalman filter that smooths the code observations with the carrier phase (divergence-free), satellite orbits and clocks, troposphere, antenna/tide/Shapiro/wind-up corrections, an L3 + P3 least-squares solution (static, or a kinematic track with a position at every epoch; receiver clock, wet troposphere, GPS–Galileo inter-system bias, float ambiguities, DOP), ionosphere/TEC with DCB estimation, single-frequency L1/L2 solutions with Klobuchar, MP1/MP2 multipath and a broadcast-ephemeris code solution, and compares them with an optional reference coordinate.

To publish: push the folder to a GitHub repository and enable **Settings → Pages → Deploy from a branch → main / (root)**. To run locally: `python -m http.server` in the folder, then open `http://localhost:8000`.

Built by Salih Tatcı for the Satellite Geodesy course (Geomatics Engineering MSc, Istanbul Technical University). Results are for teaching purposes.

## Lisans / License

MIT, see `LICENSE`.
