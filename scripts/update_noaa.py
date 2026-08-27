import json, math, time
from datetime import datetime, timezone, timedelta
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
BASE = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"
APP = "LouisianaCoastalTemperatureTracker"
STATIONS = {
    "east": {"name":"East Coast Louisiana","station":"Grand Isle, LA","id":"8761724","lat":29.263,"lon":-89.957},
    "central": {"name":"Central Coast Louisiana","station":"Eugene Island, North of, LA","id":"8764314","lat":29.373,"lon":-91.384},
    "west": {"name":"West Coast Louisiana","station":"Calcasieu Pass, LA","id":"8768094","lat":29.768,"lon":-93.343},
}
YEARS = list(range(2020, 2027))
def api(station, product, year):
    end_date = datetime.now(timezone.utc).strftime("%Y%m%d") if year == datetime.now(timezone.utc).year else f"{year}1231"
    params = {
        "begin_date": f"{year}0101",
        "end_date": end_date,
        "station": station,
        "product": product,
        "time_zone": "gmt",
        "units": "english",
        "interval": "h",
        "format": "json",
        "application": APP,
    }
    url = BASE + "?" + urlencode(params)
    req = Request(url, headers={"User-Agent": APP + "/1.0"})
    for attempt in range(5):
        try:
            with urlopen(req, timeout=120) as r:
                obj = json.load(r)
            if "error" in obj:
                raise RuntimeError(obj["error"].get("message", str(obj["error"])))
            return obj.get("data", [])
        except Exception:
            if attempt == 4:
                raise
            time.sleep(5 * (attempt + 1))
def heat_index_f(t, rh):
    if not (math.isfinite(t) and math.isfinite(rh)):
        return None
    if t < 80:
        return t
    if rh < 40:
        return t
    h = (-42.379 + 2.04901523*t + 10.14333127*rh
         - 0.22475541*t*rh - 0.00683783*t*t
         - 0.05481717*rh*rh + 0.00122874*t*t*rh
         + 0.00085282*t*rh*rh - 0.00000199*t*t*rh*rh)
    if rh < 13 and 80 <= t <= 112:
        h -= ((13-rh)/4) * math.sqrt((17-abs(t-95))/17)
    elif rh > 85 and 80 <= t <= 87:
        h += ((rh-85)/10) * ((87-t)/5)
    return h
def parse_time(s):
    # CO-OPS GMT strings are "YYYY-MM-DD HH:MM"
    return datetime.strptime(s, "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)
def week_of_year(dt):
    day = (dt - datetime(dt.year,1,1,tzinfo=timezone.utc)).days
    return min(48, day // 7 + 1)
def mean(vals):
    vals = [x for x in vals if isinstance(x,(int,float)) and math.isfinite(x)]
    return round(sum(vals)/len(vals), 3) if vals else None
def last_completed_week():
    now = datetime.now(timezone.utc)
    start = datetime(now.year,1,1,tzinfo=timezone.utc)
    week = (now-start).days // 7
    # completed 7-day weeks are numbered 1..48
    return min(48, week)
def anomaly_projection(actual, baseline):
    diffs = [actual[i]-baseline[i] for i in range(48)
             if actual[i] is not None and baseline[i] is not None]
    # Weight the most recent eight observed weekly anomalies more heavily.
    recent = []
    for i in range(48):
        if actual[i] is not None and baseline[i] is not None:
            recent.append(actual[i]-baseline[i])
    if not recent:
        return baseline[:]
    recent = recent[-8:]
    weights = list(range(1, len(recent)+1))
    anomaly = sum(v*w for v,w in zip(recent,weights))/sum(weights)
    return [round(x+anomaly,3) if x is not None else None for x in baseline]
def build_station(meta):
    all_rows = []
    for year in YEARS:
        temps = api(meta["id"], "air_temperature", year)
        time.sleep(1)
        hums = api(meta["id"], "humidity", year)
        time.sleep(1)
        humidity = {}
        for x in hums:
            try:
                humidity[x["t"]] = float(x["v"])
            except Exception:
                pass
        for x in temps:
            try:
                t = float(x["v"])
                dt = parse_time(x["t"])
            except Exception:
                continue
            rh = humidity.get(x["t"])
            hi = heat_index_f(t, rh) if rh is not None else None
            all_rows.append((dt,t,hi))
    hist_t=[[] for _ in range(48)]
    hist_h=[[] for _ in range(48)]
    actual_t=[[] for _ in range(48)]
    actual_h=[[] for _ in range(48)]
    cutoff_week = last_completed_week()
    for dt,t,hi in all_rows:
        w=week_of_year(dt)-1
        if not 0 <= w < 48:
            continue
        y=dt.year
        if 2020 <= y <= 2025:
            hist_t[w].append(t)
            if hi is not None: hist_h[w].append(hi)
        elif y == 2026 and w < cutoff_week:
            actual_t[w].append(t)
            if hi is not None: actual_h[w].append(hi)
    baseline_t=[mean(x) for x in hist_t]
    baseline_h=[mean(x) for x in hist_h]
    observed_t=[mean(x) for x in actual_t]
    observed_h=[mean(x) for x in actual_h]
    proj_t=anomaly_projection(observed_t, baseline_t)
    proj_h=anomaly_projection(observed_h, baseline_h)
    return {
        **meta,
        "baseline_temp": baseline_t,
        "baseline_heat_index": baseline_h,
        "actual_2026_temp": observed_t,
        "actual_2026_heat_index": observed_h,
        "projected_2026_temp": [proj_t[i] if observed_t[i] is None else None for i in range(48)],
        "projected_2026_heat_index": [proj_h[i] if observed_h[i] is None else None for i in range(48)],
    }
def main():
    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "years": YEARS,
        "weeks": list(range(1,49)),
        "stations": {k: build_station(v) for k,v in STATIONS.items()},
    }
    Path("data").mkdir(exist_ok=True)
    Path("data/coastal.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("Wrote data/coastal.json")
if __name__ == "__main__":
    main()
