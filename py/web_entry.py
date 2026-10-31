"""Pyodide giris noktasi: JSON metin al, JSON metin dondur."""
import json
import traceback

from gnssfaz import pipeline


def run_json(files_json, params_json, progress=None, default_atx=None):
    try:
        res = pipeline.run(json.loads(files_json), json.loads(params_json), progress=progress,
                           default_atx=default_atx)
        return json.dumps({"ok": True, "result": res}, allow_nan=False)
    except Exception as e:  # noqa: BLE001
        msg = str(e)
        code, _, text = msg.partition("|") if "|" in msg else ("internal", "", msg)
        return json.dumps({"ok": False, "code": code, "message": text or msg,
                           "trace": traceback.format_exc()[-2000:]})
