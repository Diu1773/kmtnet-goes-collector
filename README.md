# kmtnet-goes-collector

GOES-16/19 ABI L2 ACMF(Clear Sky Mask) → CTIO 야간 운량 라벨 수집 (GitHub Actions 무인).

- 출처: NOAA GOES ABI L2 ACMF, AWS Open Data `noaa-goes16`/`noaa-goes19` (GCS 미러 폴백), 익명 접근
- 산출: `data/goes_ctio_labels.csv` — 시각별(UTC 0~11·22~23시) CTIO 반경 10/20/50 km BCM 원판평균 운량
- 재개형: file_key 중복 스킵 + 빈 시각 기록(`data/missing_hours.txt`) + 6시간 cron 릴레이(예산 320분/회)
- 관측자료(KMTNet 로그 등)는 이 공개 레포에 올리지 않음 — 조인은 로컬에서
- 자매 레포: [kmtnet-gefs-collector](https://github.com/Diu1773/kmtnet-gefs-collector)
- 본체(분석): 로컬 `기상수치모델/구름/` (private)
