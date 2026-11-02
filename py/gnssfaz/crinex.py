"""Hatanaka sikistirmali (Compact RINEX, CRINEX 1.0 / 3.0) gozlem dosyasini RINEX'e acar.

Algoritma (Hatanaka 2008, crx2rnx):
- Epok satiri bir onceki epok satirina gore metin farki olarak saklanir: bosluk = degismedi,
  '&' = bosluk, diger karakterler yenisi. Satir '>' (CRINEX 3) ya da '&' (CRINEX 1) ile
  basliyorsa tam satirdir.
- Her olcu, uydu ve olcu tipi basina 3. dereceye kadar fark alinmis tam sayi (birim 0,001):
  'm&deger' yeni bir yay baslatir (m: fark derecesi), sonraki epoklarda derece m'ye kadar
  artarak m. farklar gelir. Bos alan olcu yok demektir (sonra yeniden baslatilir).
- LLI/SSI bayraklari uydu basina metin farki olarak saklanir.
"""


def _repair(old, new):
    """Metin farkini cozer."""
    n = max(len(old), len(new))
    out = []
    for j in range(n):
        c = new[j] if j < len(new) else None
        if c is None:
            out.append(old[j])
        elif c == " ":
            out.append(old[j] if j < len(old) else " ")
        elif c == "&":
            out.append(" ")
        else:
            out.append(c)
    return "".join(out)


def _fmt(v):
    """0,001 birimli tam sayi -> RINEX F14.3."""
    neg = v < 0
    a = -v if neg else v
    s = f"{a // 1000}.{a % 1000:03d}"
    if neg:
        s = "-" + s
    return s.rjust(14)[-14:]


class _Arc:
    __slots__ = ("order", "d", "n")

    def __init__(self, order, val):
        self.order = order
        self.d = [val]
        self.n = 0

    def step(self, val):
        self.n += 1
        q = min(self.n, self.order)
        new = [0] * (q + 1)
        new[q] = val
        for i in range(q - 1, -1, -1):
            new[i] = self.d[i] + new[i + 1]
        self.d = new
        return new[0]


def crx2rnx(text):
    L = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    if not L or "COMPACT RINEX" not in L[0].upper():
        return text
    cver = L[0][:9].strip()
    v3 = cver.startswith("3")
    out = []
    i = 2
    ntype = {}
    types2 = 0
    cur = None
    rver = 0.0
    while i < len(L):
        s = L[i]
        out.append(s)
        lab = s[60:80].strip()
        if lab == "RINEX VERSION / TYPE":
            try:
                rver = float(s[:9])
            except ValueError:
                rver = 3.0 if v3 else 2.11
        elif lab == "SYS / # / OBS TYPES":
            if s[0] != " ":
                cur = s[0]
                ntype[cur] = int(s[3:6])
        elif lab == "# / TYPES OF OBSERV":
            if s[:6].strip():
                types2 = int(s[:6])
        i += 1
        if lab == "END OF HEADER":
            break
    sat_col = 41 if v3 else 32
    ep_old = ""
    arcs = {}          # (sat, j) -> _Arc
    flags = {}         # sat -> flag string
    n = len(L)
    while i < n:
        s = L[i]
        i += 1
        if not s and i >= n:
            break
        if (v3 and s.startswith(">")) or (not v3 and s.startswith("&")):
            ep = s if v3 else " " + s[1:]
        else:
            ep = _repair(ep_old, s)
        ep_old = ep
        try:
            flag = int(ep[31] if v3 else ep[28]) if len(ep) > (31 if v3 else 28) and (ep[31] if v3 else ep[28]).strip() else 0
        except ValueError:
            flag = 0
        if flag > 1:
            try:
                nrec = int(ep[32:35] if v3 else ep[29:32])
            except ValueError:
                nrec = 0
            out.append(ep[:sat_col].rstrip())
            out.extend(L[i:i + nrec])
            i += nrec
            ep_old = ""
            continue
        sats_txt = ep[sat_col:]
        sats = [sats_txt[k:k + 3] for k in range(0, len(sats_txt), 3) if sats_txt[k:k + 3].strip()]
        i += 1                                   # alici saat duzeltmesi satiri (kullanilmiyor)
        lines_out = []
        for sv in sats:
            d = L[i] if i < n else ""
            i += 1
            sv = sv.replace(" ", "0") if sv[0] != " " else "G" + sv[1:].replace(" ", "0")
            nt = ntype.get(sv[0], 0) if v3 else types2
            parts = d.split(" ", nt)
            fields = parts[:nt] + [""] * (nt - len(parts[:nt]))
            ftxt = parts[nt] if len(parts) > nt else ""
            fl = _repair(flags.get(sv, ""), ftxt)
            flags[sv] = fl
            vals = []
            for j in range(nt):
                f = fields[j]
                key = (sv, j)
                if not f:
                    arcs.pop(key, None)
                    vals.append(None)
                    continue
                if "&" in f:
                    o, v = f.split("&", 1)
                    a = _Arc(int(o), int(v))
                    arcs[key] = a
                    vals.append(a.d[0])
                else:
                    a = arcs.get(key)
                    if a is None:
                        vals.append(None)
                        continue
                    vals.append(a.step(int(f)))
            cells = []
            for j, v in enumerate(vals):
                lli = fl[2 * j] if 2 * j < len(fl) else " "
                ssi = fl[2 * j + 1] if 2 * j + 1 < len(fl) else " "
                if v is None:
                    cells.append(" " * 14 + lli + ssi if (lli.strip() or ssi.strip()) else " " * 16)
                else:
                    cells.append(_fmt(v) + lli + ssi)
            if v3:
                lines_out.append((sv + "".join(cells)).rstrip())
            else:
                for k in range(0, max(nt, 1), 5):
                    lines_out.append("".join(cells[k:k + 5]).rstrip())
        # olmayan uydularin durumu sifirlanir
        present = set(x.replace(" ", "0") for x in sats)
        for key in [k for k in arcs if k[0] not in present]:
            arcs.pop(key, None)
        if v3:
            out.append(ep[:35].rstrip() if len(ep) >= 35 else ep.rstrip())
        else:
            head = ep[:32]
            out.append(head + "".join(sats[:12]))
            for k in range(12, len(sats), 12):
                out.append(" " * 32 + "".join(sats[k:k + 12]))
        out.extend(lines_out)
    return "\n".join(out) + "\n"
