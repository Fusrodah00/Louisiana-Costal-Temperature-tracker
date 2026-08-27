import csv
import gzip
import io
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen
# NOAA/NDBC source:
# Historical standard meteorological files: 2020-2025
# 2026 quality-controlled file + current realtime observations
# NDBC standard met fields include ATMP (air temp) and DEWP (dew point).
# Relative humidity is calculated from ATMP + DEWP, then NWS heat index is calculated.
BASE_HIST = "https://www.ndbc.noaa.gov/view_text_file.php"
BASE_REALTIME = "https://www.ndbc.noaa.gov/data/realtime2"
STATIONS = {
    "east": {
        "name": "East Coast Louisiana",
        "station": "Grand Isle, LA",
        "id": "8761724",
        "ndbc": "gisl1",
        "lat": 29.265,
        "lon": -89.958,
    },
    "central": {
        "name": "Central Coast Louisiana",
        "station": "North of Eugene Island, LA",
        "id": "8764314",
        "ndbc": "einl1",
        "lat": 29.373,
        "lon": -91.384,
    },
    "west": {
        "name": "West Coast Louisiana",
        "station": "Calcasieu Pass, LA",
        "id": "8768094",
        "ndbc": "capl1",
        "lat": 29.768,
        "lon": -93.343,
    },
}
HIST_YEARS = list(range(2020, 2026))
ALL_YEARS = list(range(2020, 2027))
WEEKS = 48
def fetch_bytes(url, timeout=120, attempts=5):
    req = Request(
        url,
        headers={
            "User-Agent": "LouisianaCoastalTemperatureTracker/2.0 "
                          "(NOAA/NDBC data processing)"
        },
    )
    last_error = None
    for attempt in range(attempts):
        try:
            with urlopen(req, timeout=timeout) as response:
                return response.read()
        except Exception as exc:
            last_error = exc
            if attempt < attempts - 1:
                time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"Download failed: {url}\n{last_error}")
def fetch_historical(station, year):
    filename = f"{station}h{year}.txt.gz"
    url = (
        BASE_HIST
        + f"?filename={filename}&dir=data/historical/stdmet/"
    )
    raw = fetch_bytes(url)
    try:
        raw = gzip.decompress(raw)
    except OSError:
        # Some NDBC endpoints can return the uncompressed text.
        pass
    return raw.decode("utf-8", errors="replace")
def fetch_2026_realtime(station):
    url = f"{BASE_REALTIME}/{station}.txt"
    raw = fetch_bytes(url)
    return raw.decode("utf-8", errors="replace")
def parse_ndbc_text(text):
    """
    Parse NDBC standard meteorological text.
    Header normally contains:
    #YY MM DD hh mm ... ATMP WTMP DEWP ...
    We only need timestamp, ATMP, and DEWP.
    """
    lines = text.splitlines()
    header = None
    for line in lines:
        if line.startswith("#YY"):
            header = line.lstrip("#").split()
            break
    if not header:
        return []
    # Handle the second units header by simply parsing rows below the header.
    try:
        iy = header.index("YY")
        im = header.index("MM")
        iday = header.index("DD")
        ih = header.index("hh")
        imin = header.index("mm")
        iatmp = header.index("ATMP")
        idewp = header.index("DEWP")
    except ValueError:
        return []
    rows = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) <= max(iy, im, iday, ih, imin, iatmp, idewp):
            continue
        try:
            yy = int(parts[iy])
            mm = int(parts[im])
            dd = int(parts[iday])
            hh = int(parts[ih])
            minute = int(parts[imin])
            # NDBC uses two-digit years in standard met files.
            year = yy if yy >= 1000 else (2000 + yy if yy < 70 else 1900 + yy)
            temp = float(parts[iatmp])
            dew = float(parts[idewp])
            if not (math.isfinite(temp) and math.isfinite(dew)):
                continue
            if temp <= -90 or dew <= -90:
                continue
            dt = datetime(year, mm, dd, hh, minute, tzinfo=timezone.utc)
            rows.append((dt, temp, dew))
        except (ValueError, TypeError):
            continue
    return rows
def relative_humidity_from_dewpoint(temp_f, dewpoint_f):
    """Magnus approximation, sufficient for heat-index calculation."""
    if not (math.isfinite(temp_f) and math.isfinite(dewpoint_f)):
        return None
    # Convert F to C.
    t = (temp_f - 32.0) * 5.0 / 9.0
    td = (dewpoint_f - 32.0) * 5.0 / 9.0
    a = 17.625
    b = 243.04
    gamma_t = a * t / (b + t)
    gamma_td = a * td / (b + td)
    rh = 100.0 * math.exp(gamma_td - gamma_t)
    return max(0.0, min(100.0, rh))
def heat_index_f(temp_f, rh):
    """
    NWS heat-index equation in Fahrenheit.
    For conditions below the standard heat-index threshold,
    use the actual temperature.
    """
    if not (math.isfinite(temp_f) and math.isfinite(rh)):
        return None
    if temp_f < 80.0:
        return temp_f
    # Rothfusz regression is intended for RH >= 40%.
    if rh < 40.0:
        return temp_f
    hi = (
        -42.379
        + 2.04901523 * temp_f
        + 10.14333127 * rh
        - 0.22475541 * temp_f * rh
        - 0.00683783 * temp_f * temp_f
        - 0.05481717 * rh * rh
        + 0.00122874 * temp_f * temp_f * rh
        + 0.00085282 * temp_f * rh * rh
        - 0.00000199 * temp_f * temp_f * rh * rh
    )
    # NWS low-RH adjustment.
    if rh < 13.0 and 80.0 <= temp_f <= 112.0:
        hi -= (
            ((13.0 - rh) / 4.0)
            * math.sqrt((17.0 - abs(temp_f - 95.0)) / 17.0)
        )
    # NWS high-RH adjustment.
    elif rh > 85.0 and 80.0 <= temp_f <= 87.0:
        hi += ((rh - 85.0) / 10.0) * ((87.0 - temp_f) / 5.0)
    return hi
def week_of_year(dt):
    # Seven-day bins starting Jan 1. Weeks 1-48 are displayed.
    day = (dt - datetime(dt.year, 1, 1, tzinfo=timezone.utc)).days
    return min(WEEKS, day // 7 + 1)
def last_completed_week():
    now = datetime.now(timezone.utc)
    start = datetime(now.year, 1, 1, tzinfo=timezone.utc)
    return min(WEEKS, (now - start).days // 7)
def mean(values):
    values = [
        x for x in values
        if isinstance(x, (int, float)) and math.isfinite(x)
    ]
    return round(sum(values) / len(values), 3) if values else None
def load_station_rows(meta):
    rows = []
    # 2020-2025 quality-controlled historical standard meteorological data.
    for year in HIST_YEARS:
        text = fetch_historical(meta["ndbc"], year)
        rows.extend(parse_ndbc_text(text))
        time.sleep(0.5)
    # 2026 quality-controlled historical file.
    # It is maintained through the most recently completed QC month.
    try:
        text = fetch_historical(meta["ndbc"], 2026)
        rows.extend(parse_ndbc_text(text))
    except Exception as exc:
        print(f"Warning: 2026 historical file unavailable: {exc}")
    # Add current realtime data so 2026 reaches the latest observations.
    try:
        text = fetch_2026_realtime(meta["ndbc"])
        rows.extend(parse_ndbc_text(text))
    except Exception as exc:
        print(f"Warning: realtime data unavailable: {exc}")
    # Deduplicate timestamps. Realtime data can overlap the historical file.
    unique = {}
    for dt, temp, dew in rows:
        unique[dt] = (dt, temp, dew)
    return sorted(unique.values(), key=lambda x: x[0])
def anomaly_projection(actual, baseline):
    """
    Project the remaining 2026 weeks from the 2020-2025 seasonal baseline
    using the observed 2026 anomaly. The latest 8 observed weeks are
    weighted increasingly heavily so the forecast follows the current
    2026 pattern without simply copying the last week's temperature.
    """
    recent = []
    for i in range(WEEKS):
        if actual[i] is not None and baseline[i] is not None:
            recent.append(actual[i] - baseline[i])
    if not recent:
        return baseline[:]
    recent = recent[-8:]
    weights = list(range(1, len(recent) + 1))
    anomaly = sum(v * w for v, w in zip(recent, weights)) / sum(weights)
    return [
        round(x + anomaly, 3) if x is not None else None
        for x in baseline
    ]
def build_station(meta):
    rows = load_station_rows(meta)
    hist_t = [[] for _ in range(WEEKS)]
    hist_hi = [[] for _ in range(WEEKS)]
    actual_t = [[] for _ in range(WEEKS)]
    actual_hi = [[] for _ in range(WEEKS)]
    cutoff = last_completed_week()
    for dt, temp, dew in rows:
        if dt.year not in ALL_YEARS:
            continue
        w = week_of_year(dt) - 1
        if not 0 <= w < WEEKS:
            continue
        hi = heat_index_f(
            temp,
            relative_humidity_from_dewpoint(temp, dew)
        )
        if 2020 <= dt.year <= 2025:
            hist_t[w].append(temp)
            if hi is not None:
                hist_hi[w].append(hi)
        elif dt.year == 2026 and w < cutoff:
            actual_t[w].append(temp)
            if hi is not None:
                actual_hi[w].append(hi)
    baseline_t = [mean(x) for x in hist_t]
    baseline_hi = [mean(x) for x in hist_hi]
    actual_2026_t = [mean(x) for x in actual_t]
    actual_2026_hi = [mean(x) for x in actual_hi]
    projected_t = anomaly_projection(actual_2026_t, baseline_t)
    projected_hi = anomaly_projection(actual_2026_hi, baseline_hi)
    return {
        **meta,
        "baseline_temp": baseline_t,
        "baseline_heat_index": baseline_hi,
        "actual_2026_temp": actual_2026_t,
        "actual_2026_heat_index": actual_2026_hi,
        "projected_2026_temp": [
            projected_t[i] if actual_2026_t[i] is None else None
            for i in range(WEEKS)
        ],
        "projected_2026_heat_index": [
            projected_hi[i] if actual_2026_hi[i] is None else None
            for i in range(WEEKS)
        ],
    }
def main():
    generated = datetime.now(timezone.utc).isoformat()
    cutoff = last_completed_week()
    out = {
        "generated_utc": generated,
        "source": "NOAA/NDBC standard meteorological observations",
        "historical_years": HIST_YEARS,
        "actual_year": 2026,
        "actual_through_completed_week": cutoff,
        "weeks": list(range(1, WEEKS + 1)),
        "stations": {},
    }
    for key, meta in STATIONS.items():
        print(f"Building {meta['name']}...")
        out["stations"][key] = build_station(meta)
    Path("data").mkdir(exist_ok=True)
    Path("data/coastal.json").write_text(
        json.dumps(out, indent=2),
        encoding="utf-8",
    )
    print("Wrote data/coastal.json")


if __name__ == "__main__":
    main()
