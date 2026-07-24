# -*- coding: utf-8 -*-
"""GOES ABI L2 ACMF 운량 라벨 수집기 — CTIO 점 추출 (GitHub Actions 무인 실행용).
날짜 범위 × 야간 시각(UTC 0~11,22,23)을 받아 data/goes_ctio_labels.csv에 append.
file_key 중복 스킵 = 재개형. 예산(분) 소진 시 정상 종료 → cron 릴레이가 이어받음.
출처: NOAA GOES-16/19 ABI L2 ACMF (AWS Open Data / GCS 미러, 익명).
로컬 기상수치모델/구름/scripts/collect_goes_labels.py와 동일 추출 경로 (반경 10/20/50km BCM 원판평균).
"""
import argparse, csv, math, os, re, subprocess, sys, tempfile, time
import datetime as dt
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import requests
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "data", "goes_ctio_labels.csv")
MISSING = os.path.join(HERE, "data", "missing_hours.txt")  # 파일 자체가 없는 (date,hour) 재조회 방지

LAT, LON = -30.169, -70.806            # CTIO
RADII = (10, 20, 50)
NIGHT_HOURS = list(range(0, 12)) + [22, 23]   # UTC — 칠레 겨울 밤 전체 커버 (박명 필터는 분석 단계에서)
SWITCH = dt.date(2025, 4, 4)           # GOES-19 GOES-East 운영 승계
COLS = ["scan_mid_utc", "sat", "cf10", "n10", "cf20", "n20", "cf50", "n50", "file_key"]


def bucket_for(day):
    return "noaa-goes19" if day >= SWITCH else "noaa-goes16"


GCS_MAP = {"noaa-goes16": "gcp-public-data-goes-16", "noaa-goes19": "gcp-public-data-goes-19"}


def list_hour(bucket, day, hour):
    doy = day.timetuple().tm_yday
    prefix = f"ABI-L2-ACMF/{day.year}/{doy:03d}/{hour:02d}/"
    r = requests.get(f"https://{bucket}.s3.amazonaws.com/?list-type=2&prefix={prefix}", timeout=60)
    r.raise_for_status()
    return re.findall(r"<Key>([^<]+)</Key>", r.text)


def pick_key(keys):
    good = []
    for k in sorted(keys):
        m = re.search(r"_s(\d+)_e(\d+)_", k)
        if m and m.group(1)[:11] == m.group(2)[:11]:  # 스캔 지속 0 = 이상 파일
            continue
        good.append(k)
    m6 = [k for k in good if "-M6_" in k]
    pool = m6 or good
    return pool[0] if pool else None


def download(bucket, key, dest):
    # 러너(Azure)에선 AWS가 동일 리전이라 우선, GCS 폴백
    urls = [f"https://{bucket}.s3.amazonaws.com/{key}",
            f"https://storage.googleapis.com/{GCS_MAP[bucket]}/{key}"]
    last = None
    for url in urls:
        try:
            with requests.get(url, stream=True, timeout=300) as resp:
                resp.raise_for_status()
                expected = int(resp.headers.get("Content-Length", 0))
                with open(dest + ".part", "wb") as f:
                    for chunk in resp.iter_content(1 << 20):
                        f.write(chunk)
            if expected and os.path.getsize(dest + ".part") != expected:
                os.remove(dest + ".part")
                raise IOError("incomplete download")
            os.replace(dest + ".part", dest)
            return
        except Exception as e:
            last = e
    raise last


def extract(nc_path):
    ds = xr.open_dataset(nc_path)
    p = ds["goes_imager_projection"]
    req = float(p.semi_major_axis); rpol = float(p.semi_minor_axis)
    H = float(p.perspective_point_height) + req
    lam0 = math.radians(float(p.longitude_of_projection_origin))
    e2 = (req**2 - rpol**2) / req**2
    phi = math.radians(LAT); lam = math.radians(LON)
    phic = math.atan((rpol**2 / req**2) * math.tan(phi))
    rc = rpol / math.sqrt(1 - e2 * math.cos(phic) ** 2)
    sx = H - rc * math.cos(phic) * math.cos(lam - lam0)
    sy = -rc * math.cos(phic) * math.sin(lam - lam0)
    sz = rc * math.sin(phic)
    x_t = math.asin(-sy / math.sqrt(sx**2 + sy**2 + sz**2))
    y_t = math.atan(sz / sx)
    xv = ds["x"].values.astype(np.float64); yv = ds["y"].values.astype(np.float64)
    ix = int(np.abs(xv - x_t).argmin()); iy = int(np.abs(yv - y_t).argmin())
    W = 80
    sl_y = slice(iy - W, iy + W + 1); sl_x = slice(ix - W, ix + W + 1)
    xx, yy = np.meshgrid(xv[sl_x], yv[sl_y])
    a = np.sin(xx) ** 2 + np.cos(xx) ** 2 * (np.cos(yy) ** 2 + (req**2 / rpol**2) * np.sin(yy) ** 2)
    b = -2 * H * np.cos(xx) * np.cos(yy)
    c = H**2 - req**2
    rs = (-b - np.sqrt(np.maximum(b**2 - 4 * a * c, 0))) / (2 * a)
    sx2 = rs * np.cos(xx) * np.cos(yy); sy2 = -rs * np.sin(xx); sz2 = rs * np.cos(xx) * np.sin(yy)
    lat_g = np.degrees(np.arctan((req**2 / rpol**2) * sz2 / np.sqrt((H - sx2) ** 2 + sy2**2)))
    lon_g = np.degrees(lam0 - np.arctan(sy2 / (H - sx2)))
    dphi = np.radians(lat_g - LAT); dlam = np.radians(lon_g - LON)
    hav = np.sin(dphi / 2) ** 2 + np.cos(np.radians(LAT)) * np.cos(np.radians(lat_g)) * np.sin(dlam / 2) ** 2
    dist = 2 * 6371.0 * np.arcsin(np.sqrt(hav))
    bcm = ds["BCM"].isel(y=sl_y, x=sl_x).values.astype(float)
    dqf = ds["DQF"].isel(y=sl_y, x=sl_x).values
    bcm[dqf != 0] = np.nan
    t0 = ds.attrs["time_coverage_start"]; t1 = ds.attrs["time_coverage_end"]
    mid = np.datetime64(t0[:23]) + (np.datetime64(t1[:23]) - np.datetime64(t0[:23])) / 2
    row = {"scan_mid_utc": str(mid)[:19], "sat": ds.attrs.get("platform_ID", "?")}
    for Rkm in RADII:
        m = dist <= Rkm
        vals = bcm[m]
        n = int(np.isfinite(vals).sum())
        row[f"cf{Rkm}"] = round(float(np.nanmean(vals)), 4) if n else ""
        row[f"n{Rkm}"] = n
    ds.close()
    return row


def git_flush(msg):
    subprocess.run(["git", "add", "data"], cwd=HERE, check=False)
    r = subprocess.run(["git", "commit", "-m", msg], cwd=HERE, capture_output=True, text=True)
    if r.returncode != 0:
        return
    for _ in range(5):
        if subprocess.run(["git", "push"], cwd=HERE, check=False).returncode == 0:
            return
        subprocess.run(["git", "pull", "--rebase"], cwd=HERE, check=False)
        time.sleep(5)


def process_one(args):
    bucket, key, tmpdir, idx = args
    local = os.path.join(tmpdir, f"cur_{idx}.nc")
    try:
        download(bucket, key, local)
        row = extract(local)
        row["file_key"] = key
        return row
    except Exception as e:
        return {"file_key": key, "_err": f"{type(e).__name__} {e}"}
    finally:
        if os.path.exists(local):
            os.remove(local)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--budget-min", type=float, default=0, help="0=무제한")
    ap.add_argument("--commit-every-min", type=float, default=20)
    a = ap.parse_args()
    t_start = time.time()
    budget_s = a.budget_min * 60 if a.budget_min else float("inf")

    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    done, done_dh = set(), set()
    if os.path.exists(CSV):
        with open(CSV, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                done.add(r["file_key"])
                m = re.search(r"ABI-L2-ACMF/(\d{4})/(\d{3})/(\d{2})/", r["file_key"])
                if m:
                    d = dt.date(int(m.group(1)), 1, 1) + dt.timedelta(days=int(m.group(2)) - 1)
                    done_dh.add((d, int(m.group(3))))
    missing = set()
    if os.path.exists(MISSING):
        with open(MISSING, encoding="utf-8") as f:
            missing = {tuple(l.strip().split(",")) for l in f if l.strip()}

    # 대상 (date,hour) 목록 — 이미 수집됐거나 빈 시각은 listing 자체를 스킵
    todo = []
    day = dt.date.fromisoformat(a.start)
    end = dt.date.fromisoformat(a.end)
    while day <= end:
        for hour in NIGHT_HOURS:
            if (day, hour) in done_dh: continue
            if (day.isoformat(), str(hour)) in missing: continue
            todo.append((day, hour))
        day += dt.timedelta(days=1)
    print(f"todo {len(todo)} (수집완료 {len(done_dh)}, 빈시각 {len(missing)})", flush=True)

    tmpdir = tempfile.mkdtemp(prefix="goes_")
    n_ok = n_err = 0
    last_commit = time.time()
    header_needed = not os.path.exists(CSV)
    fcsv = open(CSV, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(fcsv, fieldnames=COLS, extrasaction="ignore")
    if header_needed:
        w.writeheader()
    fmiss = open(MISSING, "a", encoding="utf-8")

    # listing은 메인스레드에서 순차로 배치 생성, 다운로드+추출은 병렬
    BATCH = a.workers * 4
    i = 0
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        while i < len(todo):
            if time.time() - t_start > budget_s:
                print("예산 소진 — 정상 종료", flush=True)
                break
            jobs = []
            while i < len(todo) and len(jobs) < BATCH:
                day, hour = todo[i]; i += 1
                bucket = bucket_for(day)
                try:
                    keys = list_hour(bucket, day, hour)
                except Exception as e:
                    print(f"{day} {hour:02d}Z listing 실패: {e}", flush=True)
                    continue
                key = pick_key(keys) if keys else None
                if key is None:
                    fmiss.write(f"{day.isoformat()},{hour}\n"); fmiss.flush()
                    continue
                if key in done: continue
                jobs.append((bucket, key, tmpdir, len(jobs)))
            for fut in as_completed([ex.submit(process_one, j) for j in jobs]):
                row = fut.result()
                if "_err" in row:
                    print(f"실패 {row['file_key']}: {row['_err']}", flush=True)
                    n_err += 1
                    continue
                w.writerow(row); fcsv.flush()
                done.add(row["file_key"]); n_ok += 1
            if n_ok and n_ok % 100 < BATCH:
                rate = n_ok / max(time.time() - t_start, 1) * 60
                print(f"{n_ok} rows ({rate:.0f}/min)", flush=True)
            if time.time() - last_commit > a.commit_every_min * 60:
                fcsv.flush(); fmiss.flush()
                git_flush(f"collect: +{n_ok} rows (진행 커밋)")
                last_commit = time.time()
    fcsv.close(); fmiss.close()
    print(f"종료: 신규 {n_ok}, 실패 {n_err}, {(time.time()-t_start)/60:.1f}분", flush=True)


if __name__ == "__main__":
    main()
